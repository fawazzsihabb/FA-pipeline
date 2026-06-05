#!/usr/bin/env python3
"""Synthesis optimisation figure for FA19 (microwave-assisted Doebner condensation).
Data: Airlangga microwave campaign (Entries 2-6) + Hoshi conventional reference."""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import os as _os
OUT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
NAVY,TEAL,AMBER,RED="#1f2d4d","#2a9d8f","#e9a522","#c44536"
plt.rcParams.update({"figure.dpi":150,"savefig.dpi":150,"font.size":11,
    "axes.grid":True,"grid.alpha":0.25,"axes.axisbelow":True,"figure.facecolor":"white",
    "axes.titleweight":"bold","axes.spines.top":False,"axes.spines.right":False})

# Real data (RECAP.pdf; Entry 2 = des-methoxy benzoyloxy analogue / morpholine only)
entries = ["E2","E3","E4","E5","E6"]
crude   = [50.8, 40.0, 44.0, 64.0, 93.0]          # % crude yield
note    = ["morpholine only\n(4-formyl remains)",
           "+ pyridine\n40% crude /17% pure",
           "+ pyridine\n4-formyl gone",
           "+ pyridine\nsuspended 12 h",
           "+ pyridine\n20 min, 1-pot"]
HOSHI = 81.99   # conventional pyrrolidine route, crude
colors = [AMBER] + [TEAL]*4   # E2 differs in base + substrate

fig,ax=plt.subplots(figsize=(9,5.8))
bars=ax.bar(entries, crude, color=colors, edgecolor="white", width=0.62, zorder=3)
ax.axhline(HOSHI, color=NAVY, ls="--", lw=1.3, zorder=2)
ax.text(4.35, HOSHI+1.2, f"Hoshi conventional route ({HOSHI:.0f}% crude)",
        ha="right", fontsize=8.5, color=NAVY)
for b,v,nt in zip(bars,crude,note):
    ax.text(b.get_x()+b.get_width()/2, v+1.5, f"{v:.0f}%", ha="center",
            fontweight="bold", fontsize=11)
    ax.text(b.get_x()+b.get_width()/2, 4, nt, ha="center", va="bottom",
            fontsize=7.2, color="#333")
# pyridine breakthrough arrow
ax.annotate("pyridine co-solvent\ndrives reaction to completion",
            xy=(1,46), xytext=(1.9,22), fontsize=8.5, color=RED, ha="center",
            arrowprops=dict(arrowstyle="->", color=RED, lw=1.3))
ax.set_ylim(0,100); ax.set_ylabel("Crude yield (%)")
ax.set_xlabel("Microwave-assisted synthesis entry (Universitas Airlangga)")
ax.set_title("Optimisation of the FA19 microwave-assisted synthesis")
fig.tight_layout(); fig.savefig(f"{OUT}/figures/make1_yield_optimization.png"); plt.close(fig)
print("Saved figures/make1_yield_optimization.png")
print("Crude yield trajectory E2->E6:", crude, "| Hoshi ref:", HOSHI)
