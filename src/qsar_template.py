#!/usr/bin/env python3
"""
QSAR extension (supervised ML) — COX-2 bioactivity prediction template.

This is a READY-TO-RUN SCAFFOLD, not a trained model. Real predictions require
labelled bioactivity data (e.g. COX-2 IC50 / Ki from ChEMBL).

Usage
-----
  Real data:  python qsar_template.py --data data/chembl_cox2.csv
              CSV must have columns: smiles, pIC50   (pIC50 = -log10(IC50 in M))

  Demo mode:  python qsar_template.py --demo
              Generates SYNTHETIC, RANDOM labels purely to confirm the pipeline
              executes. Predictions in demo mode are MEANINGLESS by construction.

Pipeline: Morgan fingerprints -> RandomForestRegressor -> 5-fold CV (R2, RMSE)
          -> predict pIC50 for the screened ferulic library.
"""
import os, sys, argparse
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.metrics import r2_score, mean_squared_error

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def featurize(smiles_list):
    X = np.zeros((len(smiles_list), 2048))
    keep = []
    for i, smi in enumerate(smiles_list):
        m = Chem.MolFromSmiles(str(smi))
        if m is None: continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        DataStructs.ConvertToNumpyArray(fp, X[i]); keep.append(i)
    return X[keep], keep

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", help="CSV with columns: smiles, pIC50")
    ap.add_argument("--demo", action="store_true", help="run with synthetic labels (illustrative only)")
    a = ap.parse_args()

    if a.demo or not a.data:
        print("="*70)
        print("DEMO MODE — synthetic random labels. Results are ILLUSTRATIVE ONLY")
        print("and must NOT be reported as real bioactivity predictions.")
        print("="*70)
        lib = pd.read_csv(os.path.join(ROOT, "data", "library.csv"))
        train = lib.copy()
        rng = np.random.default_rng(0)
        train["pIC50"] = rng.normal(6.0, 1.0, len(train)).round(2)   # FAKE labels
    else:
        train = pd.read_csv(a.data)
        assert {"smiles","pIC50"} <= set(train.columns), "CSV needs columns: smiles, pIC50"
        print(f"Loaded {len(train)} labelled COX-2 compounds from {a.data}")

    Xtr, keep = featurize(train["smiles"])
    ytr = train["pIC50"].values[keep]

    model = RandomForestRegressor(n_estimators=400, random_state=0, n_jobs=-1)
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    yhat = cross_val_predict(model, Xtr, ytr, cv=cv)
    print(f"\n5-fold CV  R2 = {r2_score(ytr, yhat):.3f} | "
          f"RMSE = {mean_squared_error(ytr, yhat)**0.5:.3f} pIC50 units")

    # fit on all labelled data, predict the screened ferulic library
    model.fit(Xtr, ytr)
    lib = pd.read_csv(os.path.join(ROOT, "data", "library.csv"))
    Xlib, keepl = featurize(lib["smiles"])
    pred = model.predict(Xlib)
    out = lib.iloc[keepl].copy()
    out["pred_pIC50"] = pred.round(2)
    out = out.sort_values("pred_pIC50", ascending=False)
    dest = os.path.join(ROOT, "data", "qsar_predictions.csv")
    out[["id","name","klass","pred_pIC50"]].to_csv(dest, index=False)
    print("\nTop predicted (mode-dependent):")
    print(out[["id","name","pred_pIC50"]].head(5).to_string(index=False))
    print(f"\nSaved: {os.path.relpath(dest, ROOT)}")
    if a.demo or not a.data:
        print("\nReminder: replace with real ChEMBL COX-2 data before using any of this.")

if __name__ == "__main__":
    main()
