#!/usr/bin/env python3
"""Figures + unsupervised cheminformatics (PCA chemical space, Butina clustering)."""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina
from sklearn.decomposition import PCA

import os as _os
OUT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
df = pd.read_csv(f"{OUT}/data/screening_results.csv")
lib = df[df["is_reference"]!=True].copy()
ref = df[df["is_reference"]==True].iloc[0]

# style
plt.rcParams.update({
    "figure.dpi":150,"savefig.dpi":150,"font.size":11,"axes.edgecolor":"#444",
    "axes.linewidth":0.8,"axes.grid":True,"grid.alpha":0.25,"grid.linewidth":0.6,
    "axes.titleweight":"bold","axes.titlesize":12,"figure.facecolor":"white",
    "axes.spines.top":False,"axes.spines.right":False,
})
NAVY,TEAL,AMBER,RED,GREY = "#1f2d4d","#2a9d8f","#e9a522","#c44536","#9aa0a6"

# ---- FIG 1: property distributions with thresholds -------------------------
fig,axs = plt.subplots(2,2,figsize=(10,7.2))
specs = [("MW","Molecular weight (Da)",500,axs[0,0]),
         ("LogP","cLogP",5,axs[0,1]),
         ("TPSA","TPSA (Å²)",140,axs[1,0]),
         ("QED","QED (drug-likeness)",None,axs[1,1])]
for col,label,thr,ax in specs:
    ax.hist(lib[col],bins=12,color=TEAL,edgecolor="white",alpha=0.85)
    ax.axvline(lib[lib["id"]=="FA19"][col].iloc[0],color=RED,lw=2,ls="--",label="FA19 lead")
    if thr is not None:
        ax.axvline(thr,color=NAVY,lw=1.4,ls=":",label="rule limit")
    ax.set_xlabel(label); ax.set_ylabel("count"); ax.legend(fontsize=8,frameon=False)
fig.suptitle("Molecular property distributions — ferulic acid derivative library (n=24)",
             fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(f"{OUT}/figures/fig1_property_distributions.png"); plt.close(fig)

# ---- FIG 2: QED vs SA (drug-likeness vs synthesizability) ------------------
fig,ax = plt.subplots(figsize=(8.6,6.2))
sizes = (lib["MW"]/lib["MW"].max()*260)+30
sc = ax.scatter(lib["SAscore"],lib["QED"],s=sizes,c=lib["DrugLikeScore"],
                cmap="viridis",edgecolor="#333",linewidth=0.6,alpha=0.9,zorder=3)
plt.colorbar(sc,label="composite drug-likeness score")
# quadrant guides
ax.axvline(lib["SAscore"].median(),color=GREY,ls=":",lw=1)
ax.axhline(0.5,color=GREY,ls=":",lw=1)
# annotate top-5 and FA19
top5 = lib.nsmallest(0,"rank") if False else lib.sort_values("rank").head(5)
for _,r in top5.iterrows():
    ax.annotate(r["id"],(r["SAscore"],r["QED"]),fontsize=8,fontweight="bold",
                xytext=(4,4),textcoords="offset points")
fa = lib[lib["id"]=="FA19"].iloc[0]
ax.scatter([fa["SAscore"]],[fa["QED"]],s=320,marker="*",color=RED,
           edgecolor="black",linewidth=1.0,zorder=5,label="FA19 lead")
ax.annotate("FA19",(fa["SAscore"],fa["QED"]),fontsize=10,fontweight="bold",
            color=RED,xytext=(6,-12),textcoords="offset points")
ax.set_xlabel("Synthetic accessibility (SA score — lower = easier)")
ax.set_ylabel("QED (drug-likeness)")
ax.set_title("Drug-likeness vs. synthesizability\n(top-left = ideal: easy to make + drug-like)")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/figures/fig2_qed_vs_sa.png"); plt.close(fig)

# ---- Fingerprints for ML ---------------------------------------------------
def fp(smi):
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi),2,nBits=2048)
all_df = pd.concat([lib,pd.DataFrame([ref])],ignore_index=True)
fps = [fp(s) for s in all_df["smiles"]]
arr = np.zeros((len(fps),2048))
for i,f in enumerate(fps): DataStructs.ConvertToNumpyArray(f,arr[i])

# ---- FIG 3: PCA chemical space ---------------------------------------------
pca = PCA(n_components=2, random_state=0)
xy = pca.fit_transform(arr)
all_df["PC1"],all_df["PC2"] = xy[:,0],xy[:,1]
classes = sorted(lib["klass"].unique())
palette = plt.cm.tab10(np.linspace(0,1,len(classes)))
cmap = dict(zip(classes,palette))
fig,ax = plt.subplots(figsize=(9,6.6))
for k in classes:
    sub = all_df[(all_df["klass"]==k)&(all_df["is_reference"]!=True)]
    ax.scatter(sub["PC1"],sub["PC2"],s=90,color=cmap[k],edgecolor="#333",
               linewidth=0.5,alpha=0.9,label=k,zorder=3)
# reference + lead overlays
r = all_df[all_df["is_reference"]==True].iloc[0]
ax.scatter([r["PC1"]],[r["PC2"]],s=240,marker="D",color="black",zorder=6)
ax.annotate("Celecoxib (ref)",(r["PC1"],r["PC2"]),fontsize=9,xytext=(8,4),textcoords="offset points")
fa = all_df[all_df["id"]=="FA19"].iloc[0]
ax.scatter([fa["PC1"]],[fa["PC2"]],s=340,marker="*",color=RED,edgecolor="black",linewidth=1,zorder=7)
ax.annotate("FA19",(fa["PC1"],fa["PC2"]),fontsize=10,fontweight="bold",color=RED,
            xytext=(8,-12),textcoords="offset points")
var = pca.explained_variance_ratio_*100
ax.set_xlabel(f"PC1 ({var[0]:.0f}% variance)"); ax.set_ylabel(f"PC2 ({var[1]:.0f}% variance)")
ax.set_title("Chemical space (PCA on Morgan fingerprints)")
ax.legend(fontsize=8,frameon=False,loc="best",ncol=2)
fig.tight_layout(); fig.savefig(f"{OUT}/figures/fig3_chemical_space_pca.png"); plt.close(fig)

# ---- Butina clustering -----------------------------------------------------
n = len(fps)
dists=[]
for i in range(1,n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
    dists.extend([1-s for s in sims])
clusters = Butina.ClusterData(dists,n,0.45,isDistData=True)
cluster_of = {}
for ci,members in enumerate(clusters):
    for m in members: cluster_of[m]=ci+1
all_df["cluster"] = [cluster_of[i] for i in range(n)]
all_df[["id","name","klass","cluster","PC1","PC2","DrugLikeScore"]].to_csv(
    f"{OUT}/data/chemical_space.csv",index=False)

# ---- FIG 4: ranked funnel bar ----------------------------------------------
topN = lib.sort_values("DrugLikeScore",ascending=True).tail(15)
fig,ax = plt.subplots(figsize=(8.6,6.4))
colors = [RED if i=="FA19" else TEAL for i in topN["id"]]
ax.barh(topN["id"]+" — "+topN["name"].str.slice(0,22),topN["DrugLikeScore"],
        color=colors,edgecolor="white")
ax.set_xlabel("composite drug-likeness score")
ax.set_title("Candidate ranking (top 15) — FA19 highlighted")
for y,(v) in enumerate(topN["DrugLikeScore"]):
    ax.text(v+0.005,y,f"{v:.2f}",va="center",fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/figures/fig4_ranking.png"); plt.close(fig)

print("Figures written:")
for f in sorted(os.listdir(f"{OUT}/figures")): print("  figures/"+f)
print(f"\nPCA variance explained: PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%  (cumulative {var[:2].sum():.1f}%)")
print(f"Butina clusters (cutoff 0.45): {len(clusters)}")
print("FA19 cluster:",all_df[all_df['id']=='FA19']['cluster'].iloc[0])
