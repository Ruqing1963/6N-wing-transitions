import csv, numpy as np, matplotlib, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})

def rd(n): return list(csv.DictReader(open(os.path.join(_DATA,n))))

# ============ FIGURE 1 : memory structure ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) step autocorrelation vs lag
ac=rd("markov_autocorr.csv")
lag=np.array([int(r["lag"]) for r in ac]); r=np.array([float(r["step_autocorr"]) for r in ac])
cols=["#d62728" if abs(v)>0.01 else ("#ff7f0e" if abs(v)>0.002 else "0.6") for v in r]
ax[0].axhline(0,color="0.5",lw=0.8)
ax[0].bar(lag,r,color=cols,width=0.6)
ax[0].annotate("lag 1\n$-0.033$ (168$\\sigma$)",(1,-0.033),fontsize=7,ha="center",va="top",color="#d62728")
ax[0].annotate("lag 2\n$+0.005$ (24$\\sigma$)",(2,0.005),fontsize=7,ha="center",va="bottom",color="#ff7f0e")
ax[0].annotate("lag $\\geq3$: $\\approx 0$",(5.5,-0.006),fontsize=7.5,ha="center",color="0.4")
ax[0].set_xlabel("lag (steps)"); ax[0].set_ylabel("switch-sequence autocorrelation")
ax[0].set_title("(A) memory decay: dominant at lag 1",fontsize=9)
ax[0].set_xticks(lag)
# (B) incremental memory by order
md=rd("markov_memory_depth.csv")
order=np.array([int(r["order"]) for r in md]); im=np.array([float(r["incremental_memory_mean"]) for r in md])
ax[1].plot(order,im,"o-",ms=7,color="#1f77b4")
for o,v in zip(order,im): ax[1].annotate("%.4f"%v,(o,v),fontsize=7,xytext=(5,4),textcoords="offset points")
ax[1].axhline(0.0009,color="0.6",lw=0.8,ls=":")
ax[1].annotate("finite-sample floor",(4.5,0.0012),fontsize=6.5,ha="center",color="0.5")
ax[1].set_xlabel("Markov order $m$"); ax[1].set_ylabel("incremental memory (oldest symbol)")
ax[1].set_title("(B) memory halves per order (not abrupt)",fontsize=9)
ax[1].set_xticks(order); ax[1].set_ylim(0,0.036)
fig.suptitle(r"Wing-transition memory on the $6N$ skeleton: not first-order Markov",fontsize=10)
fig.savefig(os.path.join(_FIG,"p34_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 : unification + scale decay ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
# (A) gap conditional curve E[g_{n+1}|g_n]  (anti-correlation / memory)
gc=rd("markov_gapcond.csv")
xm=np.array([(float(r["g_lo"])+float(r["g_hi"]))/2 for r in gc]); en=np.array([float(r["E_next_g"]) for r in gc])
ax2[0].axhline(1.0,color="0.5",lw=0.8,ls="--",label="independent (Cramer)")
ax2[0].plot(xm,en,"o-",ms=5,color="#2ca02c")
ax2[0].set_xlabel(r"current normalized gap $g_n=d_n/\log p_n$")
ax2[0].set_ylabel(r"$E[\,g_{n+1}\mid g_n\,]$")
ax2[0].set_title("(A) gap memory: big gap $\\to$ smaller next",fontsize=9)
uni=rd("markov_unification.csv"); ud={r["quantity"]:float(r["value"]) for r in uni}
ax2[0].text(0.5,0.06,"step lag-1 $=%.4f$\ngap lag-1 $=%.4f$\n(same effect)"%(ud["step_lag1_autocorr"],ud["gap_normalized_lag1_corr"]),
            transform=ax2[0].transAxes,fontsize=7,va="bottom",ha="left",
            bbox=dict(boxstyle="round",fc="white",ec="0.7",alpha=0.9))
ax2[0].legend(fontsize=7,loc="upper right")
# (B) scale decay of 1st & 2nd order signals
dec=rd("markov_decay.csv")
mid=np.array([ (float(r["win_lo"])*float(r["win_hi"]))**0.5 for r in dec])
psame=np.array([float(r["P_same"]) for r in dec]); s2=np.array([abs(float(r["second_order_signal"])) for r in dec])
ax2[1].semilogx(mid,0.5-psame,"o-",ms=5,color="#1f77b4",label=r"1st order $\frac{1}{2}-P(\mathrm{same})$")
ax2[1].semilogx(mid,s2,"s-",ms=5,color="#d62728",label=r"2nd order $|P(W|W)-P(W|S)|$")
ax2[1].set_xlabel(r"prime size (window centre)"); ax2[1].set_ylabel("signal strength")
ax2[1].set_title("(B) both biases fade slowly with scale",fontsize=9)
ax2[1].legend(fontsize=7,loc="upper right"); ax2[1].set_ylim(0,0.105)
fig2.savefig(os.path.join(_FIG,"p34_fig2.pdf")); print("fig2 done")
