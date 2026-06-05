#!/usr/bin/env python3
"""Analyse the 50 ns YASARA MD trajectory of the FA19-COX2 complex.
Parses docking3_analysis.tab -> RMSD/Rg time series, H-bond occupancy, stats."""
import re, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os
OUT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
os.makedirs(f"{OUT}/figures", exist_ok=True)
TAB = _os.path.join(OUT,"data","docking3_analysis.tab")

rows = open(TAB).read().splitlines()
hdr = rows[0].split()
col = {n:i for i,n in enumerate(hdr)}
atm_cols = [i for i,n in enumerate(hdr) if re.match(r'HB\d+Atm[12]$', n)]
ncol = len(hdr)

def getf(f, name): return float(f[col[name]])

t=[]; rmsd_ca=[]; rmsd_bb=[]; rmsd_all=[]; lig_move=[]; lig_conf=[]; rg=[]
res_per_frame=[]
for line in rows[1:]:
    f=line.split()
    if len(f)<ncol: continue
    try: float(f[0])
    except: continue
    t.append(getf(f,"Time[ns]"))
    rmsd_ca.append(getf(f,"RMSDCa")); rmsd_bb.append(getf(f,"RMSDBb"))
    rmsd_all.append(getf(f,"RMSDAll"))
    lig_move.append(getf(f,"RMSDLigMove")); lig_conf.append(getf(f,"RMSDLigConf"))
    rg.append(getf(f,"RadGyration"))
    seen=set()
    for c in atm_cols:
        m=re.match(r'[A-Za-z0-9]+\.([A-Z])(\d+)\.A$', f[c])
        if m: seen.add(f"{m.group(1)}{m.group(2)}")
    res_per_frame.append(seen)

t=np.array(t); n=len(t)
arrs={k:np.array(v) for k,v in dict(rmsd_ca=rmsd_ca,rmsd_bb=rmsd_bb,rmsd_all=rmsd_all,
                                    lig_move=lig_move,lig_conf=lig_conf,rg=rg).items()}
# equilibrated window: t >= 10 ns
eq = t>=10.0
def stat(a): return float(np.mean(a[eq])), float(np.std(a[eq]))

# H-bond occupancy
from collections import Counter
occ=Counter()
for s in res_per_frame:
    for r in s: occ[r]+=1
occ_pct={r:100*c/n for r,c in occ.items()}

summary={
  "sim_length_ns": float(t[-1]), "frames": n, "frame_spacing_ps": float((t[1]-t[0])*1000),
  "equil_window_ns": [10.0, float(t[-1])],
  "RMSD_Ca_mean_sd":   stat(arrs["rmsd_ca"]),
  "RMSD_backbone_mean_sd": stat(arrs["rmsd_bb"]),
  "RMSD_allatom_mean_sd":  stat(arrs["rmsd_all"]),
  "ligand_drift_mean_sd":  stat(arrs["lig_move"]),
  "ligand_conf_mean_sd":   stat(arrs["lig_conf"]),
  "RadGyration_mean_sd":   stat(arrs["rg"]),
  "hbond_occupancy_pct": {r: round(occ_pct[r],1) for r in sorted(occ_pct,key=occ_pct.get,reverse=True)[:6]},
}
json.dump(summary, open(f"{OUT}/data/md_summary.json","w"), indent=2)

# ---- styling ----
plt.rcParams.update({"figure.dpi":150,"savefig.dpi":150,"font.size":11,
    "axes.grid":True,"grid.alpha":0.25,"axes.titleweight":"bold","figure.facecolor":"white",
    "axes.spines.top":False,"axes.spines.right":False})
NAVY,TEAL,AMBER,RED,PURP="#1f2d4d","#2a9d8f","#e9a522","#c44536","#6a4c93"

# ---- combined 2x2 panel ----
fig,ax=plt.subplots(2,2,figsize=(11,8))
# (a) protein RMSD
ax[0,0].plot(t,arrs["rmsd_ca"],color=NAVY,lw=1,label="C-alpha")
ax[0,0].plot(t,arrs["rmsd_bb"],color=TEAL,lw=0.9,alpha=0.8,label="backbone")
ax[0,0].plot(t,arrs["rmsd_all"],color=AMBER,lw=0.8,alpha=0.7,label="all atoms")
m,s=summary["RMSD_Ca_mean_sd"]
ax[0,0].axhline(m,color=RED,ls="--",lw=1,alpha=0.8)
ax[0,0].set_title("(a) Protein backbone RMSD"); ax[0,0].set_xlabel("time (ns)")
ax[0,0].set_ylabel("RMSD (Angstrom)"); ax[0,0].legend(fontsize=8,frameon=False)
ax[0,0].text(0.97,0.06,f"Ca mean {m:.2f}±{s:.2f} A (10-50 ns)",transform=ax[0,0].transAxes,
             ha="right",fontsize=8,color=RED)
# (b) ligand RMSD
ax[0,1].plot(t,arrs["lig_move"],color=PURP,lw=1,label="ligand drift (vs receptor)")
ax[0,1].plot(t,arrs["lig_conf"],color=TEAL,lw=1,label="ligand internal conformation")
ax[0,1].set_title("(b) Ligand RMSD"); ax[0,1].set_xlabel("time (ns)")
ax[0,1].set_ylabel("RMSD (Angstrom)"); ax[0,1].legend(fontsize=8,frameon=False)
# (c) Rg
ax[1,0].plot(t,arrs["rg"],color=NAVY,lw=1)
mr,sr=summary["RadGyration_mean_sd"]
ax[1,0].set_title("(c) Radius of gyration (compactness)"); ax[1,0].set_xlabel("time (ns)")
ax[1,0].set_ylabel("Rg (Angstrom)")
ax[1,0].text(0.97,0.9,f"{mr:.2f}±{sr:.2f} A",transform=ax[1,0].transAxes,ha="right",fontsize=8)
# (d) H-bond occupancy
keep=[(r,occ_pct[r]) for r in sorted(occ_pct,key=occ_pct.get,reverse=True) if occ_pct[r]>=1][:6]
labs=[r for r,_ in keep]; vals=[v for _,v in keep]
bars=ax[1,1].bar(labs,vals,color=[TEAL if v>=50 else AMBER for v in vals],edgecolor="white")
ax[1,1].set_title("(d) Protein-ligand H-bond occupancy"); ax[1,1].set_ylabel("% of frames")
ax[1,1].set_ylim(0,100)
for b,v in zip(bars,vals): ax[1,1].text(b.get_x()+b.get_width()/2,v+1.5,f"{v:.0f}%",ha="center",fontsize=8.5,fontweight="bold")
fig.suptitle("Molecular-dynamics validation of the FA19-COX2 complex (50 ns, YASARA)",fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(f"{OUT}/figures/md1_stability_panel.png"); plt.close(fig)

# ---- standalone H-bond figure (for portfolio highlight) ----
fig,ax=plt.subplots(figsize=(7.5,5))
canon={"R120":"Arg120","Y355":"Tyr355","S530":"Ser530","R513":"Arg513","Y385":"Tyr385","H90":"His90"}
keep=[(canon.get(r,r),occ_pct[r]) for r in sorted(occ_pct,key=occ_pct.get,reverse=True)][:6]
labs=[r for r,_ in keep]; vals=[v for _,v in keep]
bars=ax.bar(labs,vals,color=[TEAL if v>=50 else AMBER for v in vals],edgecolor="white")
for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+1.5,f"{v:.0f}%",ha="center",fontweight="bold")
ax.set_ylabel("H-bond occupancy (% of 50 ns)"); ax.set_ylim(0,100)
ax.set_title("FA19 engages the canonical COX-2 active-site triad")
ax.text(0.5,-0.16,"Arg120 / Tyr355 = carboxylate constriction site   |   Ser530 = aspirin-acetylation residue",
        transform=ax.transAxes,ha="center",fontsize=8,color="#555")
fig.tight_layout(); fig.savefig(f"{OUT}/figures/md2_hbond_occupancy.png"); plt.close(fig)

print("MD analysis complete.")
print(json.dumps(summary,indent=2))
