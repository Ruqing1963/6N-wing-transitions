import csv, numpy as np, matplotlib, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
from math import isqrt
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})

def rd(n): return list(csv.DictReader(open(os.path.join(_DATA,n))))
lad=rd("dodge_ladder.csv")

# small primes
P=400; s=np.ones(P+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(P)+1):
    if s[p]: s[p*p::p]=False
sp=[int(q) for q in np.nonzero(s)[0] if q>2]

# ============ FIGURE 1 ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) dodge map: which primes can divide n^(2^k)+1
degs=[2,4,8]; mods={2:4,4:8,8:16}; ys={2:3,4:2,8:1}
col={2:"#1f77b4",4:"#2ca02c",8:"#d62728"}
for d in degs:
    hit=[q for q in sp if q%mods[d]==1 and q<=200]
    dodge=[q for q in sp if q%mods[d]!=1 and q<=200]
    ax[0].scatter(hit,[ys[d]]*len(hit),s=30,color=col[d],marker="s",label=r"$n^{%d}+1$ can divide ($q\equiv1\,\mathrm{mod}\,%d$)"%(d,mods[d]))
    ax[0].scatter(dodge,[ys[d]]*len(dodge),s=8,color="0.8",marker=".")
ax[0].set_yticks([3,2,1]); ax[0].set_yticklabels([r"$n^2{+}1$",r"$n^4{+}1$",r"$n^8{+}1$"])
ax[0].set_ylim(0.5,3.5); ax[0].set_xlim(0,200)
ax[0].set_xlabel(r"prime $q$"); ax[0].set_title("(A) the dodge tightens with degree",fontsize=9)
ax[0].legend(fontsize=6,loc="upper center",ncol=1)
# (B) richness C vs degree (non-monotonic)
dd=np.array([int(r["degree"]) for r in lad]); Cp=np.array([float(r["C_product"]) for r in lad]); Ce=np.array([float(r["C_empirical"]) for r in lad])
ax[1].plot(dd,Cp,"o-",ms=7,color="#1f77b4",label=r"$C=\prod_q g(q)$ (product)")
ax[1].plot(dd,Ce,"D",ms=6,color="#d62728",label=r"measured $Q/\!\int dt/\log(t^{d}{+}1)$")
for d,c in zip(dd,Cp): ax[1].annotate(r"$n^{%d}{+}1$"%d,(d,c),fontsize=7,xytext=(6,4),textcoords="offset points")
ax[1].set_xscale("log",base=2); ax[1].set_xticks(dd); ax[1].set_xticklabels(["2","4","8"])
ax[1].set_xlabel(r"degree $d=2^k$"); ax[1].set_ylabel("richness $C$")
ax[1].set_title("(B) richness peaks at $n^4{+}1$ (non-monotonic)",fontsize=9)
ax[1].legend(fontsize=7,loc="lower center"); ax[1].set_ylim(1.2,2.9)
fig.suptitle(r"The Fermat-type dodge ladder $n^{2^k}+1$ on the $6N$ skeleton",fontsize=10)
fig.savefig(os.path.join(_FIG,"p33_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) character spectrum g(q) for n^4+1 split by q mod 8
deg=4; mod=8
qs=np.array([q for q in sp if q<=160])
g=np.where(qs%mod==1,(qs-deg)/(qs-1),qs/(qs-1))
hit=qs%mod==1
ax2[0].axhline(1.0,color="0.6",lw=0.8,ls="--")
ax2[0].plot(qs[~hit],g[~hit],"o",ms=4,color="#2ca02c",label=r"$q\not\equiv1\,(8)$: dodged, $\frac{q}{q-1}>1$")
ax2[0].plot(qs[hit],g[hit],"s",ms=5,color="#d62728",label=r"$q\equiv1\,(8)$: $\omega=4$, $\frac{q-4}{q-1}<1$")
ax2[0].set_xlabel(r"prime $q$"); ax2[0].set_ylabel(r"local factor $g(q)$ for $n^4+1$")
ax2[0].set_title(r"(A) the mod-8 quartic-character spectrum",fontsize=9)
ax2[0].legend(fontsize=7,loc="lower right")
# (B) wings: r/l for the three (degree-invariant mask)
xb=np.arange(3); rr=np.array([int(r["right_6Np1"]) for r in lad]); ll=np.array([int(r["left_6Nm1"]) for r in lad])
frac_r=rr/(rr+ll); frac_l=ll/(rr+ll)
ax2[1].bar(xb,frac_r,color="#1f77b4",label="right $6N{+}1$")
ax2[1].bar(xb,frac_l,bottom=frac_r,color="#2ca02c",label="left $6N{-}1$")
ax2[1].axhline(1/3,color="0.3",lw=0.8,ls=":"); ax2[1].axhline(1.0,color="0.3",lw=0.4)
ax2[1].set_xticks(xb); ax2[1].set_xticklabels([r"$n^2{+}1$",r"$n^4{+}1$",r"$n^8{+}1$"])
ax2[1].set_ylabel("fraction of primes"); ax2[1].set_ylim(0,1.0)
ax2[1].set_title("(B) the $1{:}2$ wing split is degree-invariant",fontsize=9)
ax2[1].legend(fontsize=7,loc="center right")
ax2[1].text(1,0.16,"right $=1/3$",fontsize=7,ha="center",color="white")
fig2.savefig(os.path.join(_FIG,"p33_fig2.pdf")); print("fig2 done")
