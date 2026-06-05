#!/usr/bin/env python3
"""
Virtual screening pipeline for ferulic acid derivatives targeting COX-2.
Phase: DESIGN / TEST (DMTA). Computes molecular properties, applies drug-likeness
filters, scores synthetic accessibility, benchmarks similarity to a known COX-2
inhibitor (celecoxib), and ranks candidates through a screening funnel.

NOTE: the library below is a REPRESENTATIVE ferulic-acid derivative set used to
demonstrate the pipeline end-to-end. Replace `data/library.csv` with the exact
FA1-FA24 SMILES to screen the real set; the rest of the pipeline is unchanged.
"""
import os, sys, json
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import (Descriptors, QED, Crippen, Lipinski, rdMolDescriptors,
                        AllChem, DataStructs, FilterCatalog)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# SA score (RDKit contrib)
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import os as _os
OUT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
os.makedirs(f"{OUT}/data", exist_ok=True)
os.makedirs(f"{OUT}/figures", exist_ok=True)

# ---------------------------------------------------------------------------
# 1. LIBRARY  - read the real FA1-FA24 set from data/library.csv
#    Expected columns: id, smiles  (optional: name, klass)
# ---------------------------------------------------------------------------
_lib = pd.read_csv(f"{OUT}/data/library.csv")
if "name"  not in _lib.columns: _lib["name"]  = _lib["id"]
if "klass" not in _lib.columns: _lib["klass"] = "derivative"
LIBRARY = list(_lib[["id","name","klass","smiles"]].itertuples(index=False, name=None))

REFERENCE = ("CELECOXIB","Celecoxib (COX-2 ref)","reference",
             "Cc1ccc(cc1)-c1cc(C(F)(F)F)nn1-c1ccc(cc1)S(N)(=O)=O")

# ---------------------------------------------------------------------------
# 2. DESCRIPTORS + FILTERS
# ---------------------------------------------------------------------------
pains = FilterCatalog.FilterCatalogParams(); pains.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
brenk = FilterCatalog.FilterCatalogParams(); brenk.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
pains_fc = FilterCatalog.FilterCatalog(pains)
brenk_fc = FilterCatalog.FilterCatalog(brenk)

def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

ref_mol = Chem.MolFromSmiles(REFERENCE[3])
ref_fp = fp(ref_mol)

def describe(cid, name, cls, smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    mw   = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd  = Lipinski.NumHDonors(m)
    hba  = Lipinski.NumHAcceptors(m)
    tpsa = Descriptors.TPSA(m)
    rotb = Lipinski.NumRotatableBonds(m)
    arom = rdMolDescriptors.CalcNumAromaticRings(m)
    fcsp3= rdMolDescriptors.CalcFractionCSP3(m)
    mr   = Crippen.MolMR(m)
    heavy= m.GetNumHeavyAtoms()
    qed  = QED.qed(m)
    sa   = sascorer.calculateScore(m)
    # Lipinski Ro5
    ro5 = sum([mw>500, logp>5, hbd>5, hba>10])
    # Veber
    veber = (rotb<=10) and (tpsa<=140)
    # Egan
    egan = (logp<=5.88) and (tpsa<=131.6)
    # alerts
    n_pains = len(pains_fc.GetMatches(m))
    n_brenk = len(brenk_fc.GetMatches(m))
    # similarity to celecoxib
    sim = DataStructs.TanimotoSimilarity(fp(m), ref_fp)
    return dict(id=cid, name=name, klass=cls, smiles=Chem.MolToSmiles(m),
                formula=rdMolDescriptors.CalcMolFormula(m),
                MW=round(mw,2), LogP=round(logp,2), HBD=hbd, HBA=hba,
                TPSA=round(tpsa,1), RotB=rotb, AromRings=arom,
                FractionCSP3=round(fcsp3,3), MolarRefr=round(mr,1), HeavyAtoms=heavy,
                QED=round(qed,3), SAscore=round(sa,2),
                Lipinski_violations=ro5, Veber_pass=veber, Egan_pass=egan,
                PAINS_alerts=n_pains, Brenk_alerts=n_brenk,
                Sim_celecoxib=round(sim,3))

rows = []
bad = []
for cid,name,cls,smi in LIBRARY:
    d = describe(cid,name,cls,smi)
    if d is None: bad.append(cid)
    else: rows.append(d)
# reference too
rref = describe(*REFERENCE); rref["is_reference"]=True

df = pd.DataFrame(rows)
df["is_reference"] = False
df["is_lead"] = df["id"]=="FA19"

print(f"Library parsed: {len(df)} derivatives | invalid SMILES: {bad}")

# ---------------------------------------------------------------------------
# 3. SCREENING FUNNEL  (drug-likeness pre-filter, before docking)
# ---------------------------------------------------------------------------
df["pass_ro5"]   = df["Lipinski_violations"]<=1
df["pass_veber"] = df["Veber_pass"]
df["pass_pains"] = df["PAINS_alerts"]==0
df["pass_prefilter"] = df["pass_ro5"] & df["pass_veber"] & df["pass_pains"]

# Composite pre-docking drug-likeness score (transparent, normalised 0-1)
def synth_ease(sa):           # SA 1 (easy) -> 1.0 ; 10 (hard) -> 0.0
    return round(max(0.0, min(1.0, (10.0-sa)/9.0)),3)
df["SynthEase"] = df["SAscore"].apply(synth_ease)
df["DrugLikeScore"] = (0.55*df["QED"] + 0.30*df["SynthEase"] + 0.15*df["pass_prefilter"].astype(int)).round(3)

df = df.sort_values("DrugLikeScore", ascending=False).reset_index(drop=True)
df.insert(0,"rank",df.index+1)

# Save master table (+ reference appended for context)
full = pd.concat([df, pd.DataFrame([rref])], ignore_index=True)
cols = ["rank","id","name","klass","smiles","formula","MW","LogP","HBD","HBA","TPSA","RotB",
        "AromRings","FractionCSP3","MolarRefr","HeavyAtoms","QED","SAscore","SynthEase",
        "Lipinski_violations","Veber_pass","Egan_pass","PAINS_alerts","Brenk_alerts",
        "Sim_celecoxib","pass_prefilter","DrugLikeScore","is_lead","is_reference"]
full = full[[c for c in cols if c in full.columns]]
full.to_csv(f"{OUT}/data/screening_results.csv", index=False)
df[["id","name","klass","smiles"]].to_csv(f"{OUT}/data/library_canonical.csv", index=False)

# ---------------------------------------------------------------------------
# 4. SUMMARY
# ---------------------------------------------------------------------------
n_pass = int(df["pass_prefilter"].sum())
lead = df[df["is_lead"]].iloc[0]
summary = {
    "n_library": int(len(df)),
    "n_pass_prefilter": n_pass,
    "n_fail_prefilter": int(len(df)-n_pass),
    "lead_id": "FA19",
    "lead_rank_by_druglikeness": int(lead["rank"]),
    "lead_QED": float(lead["QED"]),
    "lead_SAscore": float(lead["SAscore"]),
    "lead_Lipinski_violations": int(lead["Lipinski_violations"]),
    "lead_passes_prefilter": bool(lead["pass_prefilter"]),
    "top5": df.head(5)[["rank","id","name","DrugLikeScore","QED","SAscore"]].to_dict("records"),
}
with open(f"{OUT}/data/summary.json","w") as f: json.dump(summary,f,indent=2)

print("\n=== SCREENING FUNNEL ===")
print(f"Library size           : {summary['n_library']}")
print(f"Pass drug-likeness gate: {summary['n_pass_prefilter']}  (Ro5 <=1 viol, Veber, no PAINS)")
print(f"FA19 lead -> rank {summary['lead_rank_by_druglikeness']} by drug-likeness | QED {summary['lead_QED']} | SA {summary['lead_SAscore']} | Ro5 viol {summary['lead_Lipinski_violations']} | passes gate: {summary['lead_passes_prefilter']}")
print("\nTop 5 by drug-likeness score:")
print(df.head(5)[["rank","id","name","DrugLikeScore","QED","SAscore","Sim_celecoxib"]].to_string(index=False))
print("\nSaved: data/screening_results.csv, data/library.csv, data/summary.json")
