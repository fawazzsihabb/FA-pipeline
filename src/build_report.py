#!/usr/bin/env python3
"""Assemble the FA19 virtual-screening report PDF from computed results."""
import json, datetime
import pandas as pd
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak, HRFlowable, KeepTogether)
import rdkit, sklearn

import os as _os
OUT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
df = pd.read_csv(f"{OUT}/data/screening_results.csv")
summ = json.load(open(f"{OUT}/data/summary.json"))
cs = pd.read_csv(f"{OUT}/data/chemical_space.csv")
lib = df[df["is_reference"]!=True].copy()

NAVY=colors.HexColor("#1f2d4d"); TEAL=colors.HexColor("#2a9d8f")
AMBER=colors.HexColor("#e9a522"); LGREY=colors.HexColor("#eef1f5")
RED=colors.HexColor("#c44536"); DGREY=colors.HexColor("#444")

ss=getSampleStyleSheet()
def S(name,**kw):
    base=kw.pop("parent",ss["Normal"]); return ParagraphStyle(name,parent=base,**kw)
title   = S("t",fontName="Helvetica-Bold",fontSize=21,textColor=NAVY,leading=25,spaceAfter=2)
subtitle= S("s",fontName="Helvetica",fontSize=11.5,textColor=DGREY,leading=15,spaceAfter=2)
meta    = S("m",fontName="Helvetica",fontSize=8.5,textColor=colors.grey,leading=12)
h1      = S("h1",fontName="Helvetica-Bold",fontSize=13.5,textColor=NAVY,spaceBefore=14,spaceAfter=5,leading=16)
h2      = S("h2",fontName="Helvetica-Bold",fontSize=11,textColor=TEAL,spaceBefore=9,spaceAfter=3)
body    = S("b",fontSize=9.7,leading=14,alignment=TA_JUSTIFY,spaceAfter=5)
bullet  = S("bu",fontSize=9.7,leading=13.5,leftIndent=12,bulletIndent=2,spaceAfter=2)
cap     = S("c",fontSize=8,textColor=colors.grey,alignment=TA_CENTER,spaceBefore=2,spaceAfter=10,leading=10)
note    = S("n",fontSize=8.7,leading=12,textColor=DGREY,leftIndent=8,
            borderColor=AMBER,borderWidth=0,backColor=colors.HexColor("#fdf6e7"),
            borderPadding=6,spaceAfter=8)
small   = S("sm",fontSize=8,leading=10.5,textColor=colors.grey)

st=[]
def rule(c=NAVY,w=1.2): st.append(HRFlowable(width="100%",thickness=w,color=c,spaceBefore=3,spaceAfter=8))
def img(path,w):
    from reportlab.lib.utils import ImageReader
    iw,ih=ImageReader(path).getSize(); return Image(path,width=w,height=w*ih/iw)

# ---------- header ----------
st.append(Paragraph("Virtual Screening of Ferulic Acid Derivatives Targeting COX-2",title))
st.append(Paragraph("Computational property profiling, drug-likeness filtering, and chemical-space "
                    "analysis of a 24-compound library, with selection of lead candidate FA19",subtitle))
rule()
today=datetime.date.today().strftime("%d %B %Y")
st.append(Paragraph(f"DMTA phase: Design / Test &nbsp;|&nbsp; Target: cyclooxygenase-2 (COX-2) &nbsp;|&nbsp; "
                    f"Indication: breast cancer &nbsp;|&nbsp; Report generated: {today}",meta))
st.append(Spacer(1,8))

# ---------- executive summary ----------
st.append(Paragraph("Executive summary",h1))
lead=lib[lib["id"]=="FA19"].iloc[0]
st.append(Paragraph(
    f"A representative library of <b>{summ['n_library']} ferulic acid derivatives</b> spanning alkyl/aralkyl "
    f"esters, amides, phenolic acylations, bulky aroyl esters, and ring-fused analogues was profiled in silico. "
    f"All compounds were scored for physicochemical properties, drug-likeness (QED), synthetic accessibility "
    f"(SA score), and structural-alert liabilities (PAINS, Brenk), then ranked through a transparent screening "
    f"funnel. <b>{summ['n_pass_prefilter']} of {summ['n_library']} compounds cleared the drug-likeness gate</b> "
    f"(&le;1 Lipinski violation, Veber-compliant, no PAINS alerts), indicating the ferulic scaffold is "
    f"intrinsically well-behaved; the two exceptions (FA2, a catechol, and FA24) were flagged by PAINS, not by physicochemical rules.",body))
st.append(Paragraph(
    f"The lead candidate <b>FA19</b> (4-(2-carboxyvinyl)-2-methoxyphenyl 4-methoxybenzoate, C<sub>18</sub>H<sub>16</sub>O<sub>6</sub>) "
    f"ranks <b>#{summ['lead_rank_by_druglikeness']} of {summ['n_library']}</b> on drug-likeness alone "
    f"(QED {summ['lead_QED']:.2f}, SA score {summ['lead_SAscore']:.2f}, {summ['lead_Lipinski_violations']} Lipinski violations). "
    f"This is the intended result: FA19 is <b>not</b> selected for being the most drug-like molecule, but for its "
    f"predicted binding to the COX-2 hydrophobic accessory pocket (docking, reported separately). The screening "
    f"establishes that FA19 sits in clean, synthesizable, drug-like chemical space, qualifying it for the docking "
    f"and molecular-dynamics stages that justify synthesis.",body))

# KPI strip
kpi=[["Library","Pass drug-likeness gate","FA19 QED","FA19 SA score","FA19 Lipinski viol."],
     [str(summ['n_library']),f"{summ['n_pass_prefilter']} / {summ['n_library']}",
      f"{summ['lead_QED']:.2f}",f"{summ['lead_SAscore']:.2f}",str(summ['lead_Lipinski_violations'])]]
t=Table(kpi,colWidths=[3.0*cm]*5)
t.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(-1,0),NAVY),("TEXTCOLOR",(0,0),(-1,0),colors.white),
    ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,0),7.5),
    ("FONTNAME",(0,1),(-1,1),"Helvetica-Bold"),("FONTSIZE",(0,1),(-1,1),13),
    ("TEXTCOLOR",(0,1),(-1,1),TEAL),("ALIGN",(0,0),(-1,-1),"CENTER"),
    ("VALIGN",(0,0),(-1,-1),"MIDDLE"),("TOPPADDING",(0,0),(-1,-1),6),
    ("BOTTOMPADDING",(0,0),(-1,-1),6),("GRID",(0,0),(-1,-1),0.5,colors.white),
    ("BACKGROUND",(0,1),(-1,1),LGREY)]))
st.append(t); st.append(Spacer(1,4))

# ---------- methods ----------
st.append(Paragraph("1 · Objective and workflow",h1))
st.append(Paragraph(
    "This work executes the Design and Test stages of an industry-standard DMTA (Design–Make–Test–Analyse) "
    "drug-discovery cycle. The objective is to triage a focused library of ferulic acid derivatives and identify "
    "a single lead worth synthesizing, using reproducible open-source cheminformatics rather than relying solely "
    "on a commercial docking score. The funnel is intentionally ordered: <b>(i)</b> compute properties, "
    "<b>(ii)</b> apply rule-based drug-likeness and structural-alert filters, <b>(iii)</b> assess synthesizability, "
    "<b>(iv)</b> map chemical-space diversity, then <b>(v)</b> dock survivors and <b>(vi)</b> validate the winner "
    "by molecular dynamics.",body))

st.append(Paragraph("2 · Compound library",h1))
cls_counts=lib["klass"].value_counts()
ctab=[["Structural class","n"]]+[[k,str(v)] for k,v in cls_counts.items()]+[["Total",str(len(lib))]]
t=Table(ctab,colWidths=[7*cm,2*cm])
t.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(-1,0),TEAL),("TEXTCOLOR",(0,0),(-1,0),colors.white),
    ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,-1),9),
    ("ALIGN",(1,0),(1,-1),"CENTER"),("ROWBACKGROUNDS",(0,1),(-1,-2),[colors.white,LGREY]),
    ("FONTNAME",(0,-1),(-1,-1),"Helvetica-Bold"),("BACKGROUND",(0,-1),(-1,-1),colors.HexColor("#dfe5ec")),
    ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#ccd3dc")),("TOPPADDING",(0,0),(-1,-1),3.5),
    ("BOTTOMPADDING",(0,0),(-1,-1),3.5),("LEFTPADDING",(0,0),(-1,-1),6)]))
st.append(Paragraph(
    f"The library comprises {len(lib)} derivatives built on the ferulic acid core "
    "((E)-3-(4-hydroxy-3-methoxyphenyl)prop-2-enoic acid), functionalised at the carboxylic acid (esters, amides) "
    "and the phenolic hydroxyl (acylations, bulky aroyl esters), plus a ring-fused analogue. The parent ferulic "
    "acid is retained as a control and celecoxib as a reference COX-2 inhibitor.",body))
st.append(t); st.append(Spacer(1,4))
st.append(Paragraph(
    "<b>Note —</b> this is a representative demonstration library. To screen the exact professor-supplied FA1–FA24 set, "
    "replace <font face='Courier'>data/library.csv</font> with those SMILES; the pipeline is otherwise unchanged.",note))

st.append(Paragraph("3 · Computational methods",h1))
for txt in [
    f"<b>Descriptors.</b> RDKit {rdkit.__version__} computed molecular weight, cLogP (Crippen), H-bond donors/acceptors, "
    "topological polar surface area, rotatable bonds, aromatic-ring count, fraction sp3 carbon, molar refractivity, and QED.",
    "<b>Drug-likeness filters.</b> Lipinski rule-of-five (violations counted), Veber (rotatable bonds &le; 10 and "
    "TPSA &le; 140), and Egan absorption rules. The pre-docking gate requires &le; 1 Lipinski violation, Veber pass, "
    "and zero PAINS alerts.",
    "<b>Structural alerts.</b> PAINS (480 filters) and Brenk catalogues via the RDKit FilterCatalog flag "
    "pan-assay-interference and unstable/reactive substructures.",
    "<b>Synthetic accessibility.</b> Ertl–Schuffenhauer SA score (1 = easy, 10 = hard) estimates synthetic feasibility, "
    "directly informing the Make stage.",
    f"<b>Chemical space / unsupervised ML.</b> 2048-bit Morgan fingerprints (radius 2) feed principal-component analysis "
    f"(scikit-learn {sklearn.__version__}) and Butina clustering (Tanimoto distance, cutoff 0.45). Tanimoto similarity to "
    "celecoxib benchmarks novelty against a known inhibitor.",
    "<b>Composite drug-likeness score.</b> A transparent, normalised index = 0.55 x QED + 0.30 x synthetic-ease "
    "+ 0.15 x gate-pass, used only to order the pre-docking pool (binding affinity is the decisive downstream filter).",
]:
    st.append(Paragraph("&bull;&nbsp; "+txt,bullet))

st.append(PageBreak())

# ---------- results ----------
st.append(Paragraph("4 · Results",h1))

st.append(Paragraph("4.1 Property distributions and the drug-likeness gate",h2))
st.append(Paragraph(
    f"Across the library, molecular weight, cLogP, and TPSA fall well within Lipinski/Veber limits; the lead FA19 "
    f"(dashed red line) sits inside every distribution. {summ['n_pass_prefilter']} of {summ['n_library']} compounds pass the drug-likeness gate, "
    "so only the two PAINS-flagged compounds are set aside here — diversity and binding, not gross property liabilities, drive the "
    "downstream selection.",body))
st.append(img(f"{OUT}/figures/fig1_property_distributions.png",15.5*cm))
st.append(Paragraph("Figure 1. Distributions of key physicochemical properties (n=24). Dotted lines mark rule limits; "
                    "dashed red marks FA19.",cap))

st.append(KeepTogether([
    Paragraph("4.2 Drug-likeness versus synthesizability",h2),
    Paragraph(
    "Plotting QED against SA score separates candidates that are both drug-like and easy to make (upper-left) from "
    "harder or less drug-like ones. The entire library is highly synthesizable (SA score 1.7–2.1), which de-risks the "
    "Make stage. FA19 is moderate on QED but among the easier syntheses, consistent with the green-chemistry, "
    "microwave-assisted route used to prepare it.",body),
    img(f"{OUT}/figures/fig2_qed_vs_sa.png",13.5*cm),
    Paragraph("Figure 2. QED vs. SA score; marker size scales with molecular weight, colour with composite score. "
              "FA19 starred.",cap)]))

st.append(KeepTogether([
    Paragraph("4.3 Chemical-space coverage (PCA + clustering)",h2),
    Paragraph(
    f"PCA on Morgan fingerprints (PC1+PC2 capture ~43% of variance) shows the derivatives spreading by structural "
    f"class, with Butina clustering resolving 5 chemotype families at a Tanimoto cutoff of 0.45. Celecoxib lies far "
    f"from every derivative, confirming the ferulic series occupies novel chemical space rather than mimicking the "
    f"reference inhibitor. FA19 clusters with the bulky aroyl esters that were designed to reach the COX-2 accessory "
    f"pocket.",body),
    img(f"{OUT}/figures/fig3_chemical_space_pca.png",13.8*cm),
    Paragraph("Figure 3. Chemical space by PCA on Morgan fingerprints, coloured by structural class. "
              "FA19 (star) and celecoxib (diamond) highlighted.",cap)]))

st.append(KeepTogether([
    Paragraph("4.4 Candidate ranking",h2),
    img(f"{OUT}/figures/fig4_ranking.png",13.5*cm),
    Paragraph("Figure 4. Top-15 candidates by composite drug-likeness score (pre-docking). FA19 highlighted in red.",cap)]))

# ---------- results table ----------
st.append(Paragraph("4.5 Screening table (top candidates + lead)",h2))
show=lib.sort_values("rank").head(12).copy()
if "FA19" not in show["id"].values:
    show=pd.concat([show,lib[lib["id"]=="FA19"]],ignore_index=True)
hdr=["#","ID","Name","Formula","MW","cLogP","QED","SA","Ro5","Score"]
data=[hdr]
for _,r in show.iterrows():
    data.append([int(r["rank"]),r["id"],r["name"][:24],r["formula"],
                 f'{r["MW"]:.0f}',f'{r["LogP"]:.1f}',f'{r["QED"]:.2f}',
                 f'{r["SAscore"]:.2f}',int(r["Lipinski_violations"]),f'{r["DrugLikeScore"]:.2f}'])
t=Table(data,colWidths=[0.7*cm,1.4*cm,4.3*cm,2.2*cm,1.2*cm,1.3*cm,1.1*cm,1.0*cm,0.9*cm,1.2*cm],repeatRows=1)
tstyle=[("BACKGROUND",(0,0),(-1,0),NAVY),("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,-1),8),
        ("ALIGN",(0,0),(0,-1),"CENTER"),("ALIGN",(4,0),(-1,-1),"CENTER"),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,LGREY]),
        ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#ccd3dc")),
        ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),
        ("LEFTPADDING",(0,0),(-1,-1),4)]
# highlight FA19 row
for i,r in enumerate(show.sort_values("rank").itertuples(),start=1):
    if r.id=="FA19":
        tstyle += [("BACKGROUND",(0,i),(-1,i),colors.HexColor("#fbe3df")),
                   ("FONTNAME",(0,i),(-1,i),"Helvetica-Bold"),("TEXTCOLOR",(1,i),(2,i),RED)]
t.setStyle(TableStyle(tstyle))
st.append(t)
st.append(Paragraph("MW in Da; Ro5 = Lipinski violations; Score = composite pre-docking drug-likeness index. "
                    "Full table with all descriptors in data/screening_results.csv.",small))

# ---------- rationale ----------
st.append(Paragraph("5 · Why FA19 is the lead",h1))
st.append(Paragraph(
    "FA19 introduces a bulky 4-methoxybenzoyl group on the ferulic phenol, designed to occupy the hydrophobic "
    "accessory pocket that distinguishes COX-2 from COX-1 and underlies coxib selectivity. In MOE, this drops the "
    "predicted binding energy from approximately -5.5 kcal/mol (parent control) to -8.3 kcal/mol. The role of this "
    "screening report is to show that the binding gain is not bought at the expense of developability: FA19 remains "
    "drug-like (0 Lipinski violations, Veber- and Egan-compliant, no PAINS/Brenk alerts) and highly synthesizable "
    "(SA score 1.85). The combination — strong predicted binding from docking plus clean properties from screening — "
    "is what qualifies FA19 for molecular-dynamics validation and wet-lab synthesis.",body))

st.append(Paragraph("6 · Limitations and honest caveats",h1))
for txt in [
    "<b>QED definition.</b> The QED reported here (RDKit default, ~0.50 for FA19) is the reproducible open-source value. "
    "A previously cited 0.74 likely came from a different tool or weighting; any portfolio claim should state the method.",
    "<b>Representative library.</b> Results illustrate the pipeline on a representative derivative set; the exact "
    "FA1–FA24 conclusions follow once the real SMILES are loaded.",
    "<b>Docking is pending in this report.</b> The decisive affinity ranking is produced by the separate AutoDock Vina "
    "re-docking pipeline, which must include a celecoxib positive-control redock (RMSD &lt; 2 angstrom) to be trusted.",
    "<b>No experimental bioactivity.</b> No measured IC50/inhibition data are included; predicted properties are not a "
    "substitute for assays. A supervised QSAR extension is provided but requires labelled COX-2 data (e.g. ChEMBL).",
]:
    st.append(Paragraph("&bull;&nbsp; "+txt,bullet))

st.append(Paragraph("7 · Next steps",h1))
for txt in [
    "Dock all gate-passing candidates with the AutoDock Vina notebook; confirm protocol with celecoxib redock RMSD.",
    "Carry the best-binding, property-clean candidate (FA19) into YASARA molecular dynamics; extract RMSD/RMSF stability.",
    "Document the green-chemistry synthesis, purification, and NMR/FTIR characterisation as the Make-stage feasibility evidence.",
    "Compile design + dock + MD + synthesis into the industry pitch (developability + manufacturability narrative).",
]:
    st.append(Paragraph("&bull;&nbsp; "+txt,bullet))

rule(TEAL,0.8)
st.append(Paragraph(
    f"Reproducibility — generated by screen.py + plots.py (RDKit {rdkit.__version__}, scikit-learn {sklearn.__version__}). "
    "Morgan radius 2 / 2048 bits; PCA random_state=0; Butina cutoff 0.45; composite score weights stated in §3. "
    "Re-run: <font face='Courier'>python screen.py &amp;&amp; python plots.py &amp;&amp; python build_report.py</font>.",small))

doc=SimpleDocTemplate(f"{OUT}/FA19_virtual_screening_report.pdf",pagesize=A4,
                      topMargin=1.4*cm,bottomMargin=1.4*cm,leftMargin=1.6*cm,rightMargin=1.6*cm,
                      title="FA19 Virtual Screening Report",author="")
doc.build(st)
print("PDF written:",f"{OUT}/FA19_virtual_screening_report.pdf")
