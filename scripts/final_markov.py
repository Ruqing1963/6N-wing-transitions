import numpy as np, csv, os, time, math
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
from math import isqrt
t0=time.time()
N=500_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
print("sieve done, primes=%d (%.1fs)"%(len(primes),time.time()-t0))
pr=primes[primes>3]
w=(pr%6==1).astype(np.int8)               # 1=R(6N+1), 0=L(6N-1)
n=len(w)
print("primes>3=%d  R=%d L=%d  (%.1fs)"%(n,int(w.sum()),n-int(w.sum()),time.time()-t0))

# ---------- 1st/2nd/3rd order transition table ----------
def patP(order):
    # P(next=R | pattern of length `order`) for all 2^order patterns
    rows=[]
    base=w[order:]                         # the "next" symbol aligned
    cols=[w[i:n-order+i] for i in range(order)]   # pattern symbols
    for pat in range(2**order):
        bits=[(pat>>(order-1-j))&1 for j in range(order)]  # MSB = oldest
        mask=np.ones(len(base),dtype=bool)
        for j in range(order): mask &= (cols[j]==bits[j])
        cnt=int(mask.sum())
        pr_=float(base[mask].mean()) if cnt else float('nan')
        lab="".join('R' if b else 'L' for b in bits)
        rows.append((lab,cnt,pr_))
    return rows

with open(os.path.join(_DATA,"markov_transitions.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["order","pattern_oldest_to_newest","count","P_next_is_R"])
    for order in (1,2,3,4):
        for lab,cnt,p in patP(order):
            wcsv.writerow([order,lab,cnt,"%.6f"%p])
print("transition table done (%.1fs)"%(time.time()-t0))

# ---------- incremental memory at each order: effect of the OLDEST symbol given the rest ----------
# order k: avg over 2^(k-1) "rest" patterns of |P(R| R,rest) - P(R| L,rest)|
def incmem(order):
    tab={lab:(cnt,p) for lab,cnt,p in patP(order)}
    diffs=[]; wsum=0.0; wnum=0.0
    for rest in range(2**(order-1)):
        restbits="".join('R' if (rest>>(order-2-j))&1 else 'L' for j in range(order-1)) if order>1 else ""
        pr_R=tab["R"+restbits][1]; pr_L=tab["L"+restbits][1]
        nR=tab["R"+restbits][0]; nL=tab["L"+restbits][0]
        diffs.append(abs(pr_R-pr_L))
        # weight by samples for a pooled estimate
        wsum+=abs(pr_R-pr_L)*(nR+nL); wnum+=(nR+nL)
    return float(np.mean(diffs)), float(wsum/wnum)
with open(os.path.join(_DATA,"markov_memory_depth.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["order","incremental_memory_mean","incremental_memory_weighted","ratio_to_prev"])
    print("\nincremental memory by order (effect of oldest symbol given the rest):")
    prev=None
    for order in (2,3,4,5):
        mean_d,wd=incmem(order)
        ratio=(mean_d/prev) if prev else float('nan')
        wcsv.writerow([order,"%.6f"%mean_d,"%.6f"%wd,"%.4f"%ratio])
        print("  order %d: incremental memory = %.5f   ratio-to-prev = %s"%(order,mean_d,"%.3f"%ratio if prev else "-"))
        prev=mean_d

# ---------- step (switch) sequence autocorrelation to lag 8 ----------
step=(w[1:]!=w[:-1]).astype(np.float64)
ps=step.mean(); sm=step-ps; denom=float(np.sum(sm*sm)); m=len(step)
with open(os.path.join(_DATA,"markov_autocorr.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["lag","step_autocorr","sigma_significance"])
    print("\nstep autocorrelation (P(switch)=%.5f):"%ps)
    for lag in range(1,9):
        ac=float(np.sum(sm[:-lag]*sm[lag:])/denom)
        sig=abs(ac)*math.sqrt(m-lag)        # approx sigma for autocorr
        wcsv.writerow([lag,"%.6f"%ac,"%.1f"%sig])
        print("  lag %d: r=%+.5f  (%.0f sigma)"%(lag,ac,sig))

# ---------- decay of 1st & 2nd order signals by window ----------
s0=step[:-1]; s1=step[1:]; posn=pr[1:len(step)]    # position proxy = p_{n} at step i
a=w[:-1]; b=w[1:]
with open(os.path.join(_DATA,"markov_decay.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["win_lo","win_hi","n","P_same","first_order_bias","P_W_after_W","P_W_after_S","second_order_signal"])
    edges=[10**4,10**5,10**6,10**7,10**8,3*10**8,5*10**8]
    print("\nwindow            P(same)  2nd-order P(W|W)-P(W|S)")
    for i in range(len(edges)-1):
        lo,hi=edges[i],edges[i+1]
        mm=(pr[:-1]>=lo)&(pr[:-1]<hi)
        psame=float((a[mm]==b[mm]).mean()) if mm.sum() else float('nan')
        idx=np.nonzero((posn>=lo)&(posn<hi))[0]
        idx=idx[idx+1<len(step)]
        ww=float(s1[idx][s0[idx]==1].mean()); wss=float(s1[idx][s0[idx]==0].mean())
        wcsv.writerow([lo,hi,int(mm.sum()),"%.5f"%psame,"%+.5f"%(psame-0.5),"%.5f"%ww,"%.5f"%wss,"%+.5f"%(ww-wss)])
        print("  [%.0e,%.0e) %.5f   %+.5f"%(lo,hi,psame,ww-wss))

# ---------- unification: step vs gap anticorrelation ----------
d=np.diff(pr).astype(np.float64)
g=d/np.log(pr[:-1].astype(np.float64))
gap_corr=float(np.corrcoef(g[:-1],g[1:])[0,1])
step_lag1=float(np.sum(sm[:-1]*sm[1:])/denom)
# step is exactly the indicator [gap mod 6 != 0]; verify
stp_from_gap=((d.astype(np.int64))%6!=0).astype(np.int8)
agree=float((stp_from_gap==step.astype(np.int8)).mean())
with open(os.path.join(_DATA,"markov_unification.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["quantity","value","note"])
    wcsv.writerow(["step_lag1_autocorr","%.5f"%step_lag1,"wing switch-sequence lag-1"])
    wcsv.writerow(["gap_normalized_lag1_corr","%.5f"%gap_corr,"corr(g_n,g_{n+1}), g=d/log p"])
    wcsv.writerow(["raw_gap_lag1_corr","%.5f"%float(np.corrcoef(d[:-1],d[1:])[0,1]),"corr(d_n,d_{n+1})"])
    wcsv.writerow(["step_equals_gap_mod6_indicator","%.5f"%agree,"step_i == [gap_i not divisible by 6]"])
print("\nstep lag-1 = %+.5f   gap(norm) lag-1 = %+.5f   (step==[gap%%6!=0]: %.4f)"%(step_lag1,gap_corr,agree))

# gap conditional curve E[g_{n+1} | g_n bucket]  (for unification figure)
xg=g[:-1]; yg=g[1:]
with open(os.path.join(_DATA,"markov_gapcond.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["g_lo","g_hi","E_next_g","n"])
    for lo,hi in [(0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.0,2.5),(2.5,3.0),(3.0,4.0),(4.0,10.0)]:
        mm=(xg>=lo)&(xg<hi)
        if mm.sum()>1000:
            wcsv.writerow(["%.2f"%lo,"%.2f"%hi,"%.5f"%float(yg[mm].mean()),int(mm.sum())])

with open(os.path.join(_DATA,"markov_summary.csv"),"w",newline="") as fh:
    wcsv=csv.writer(fh); wcsv.writerow(["parameter","value"])
    t1=patP(1); t2=patP(2)
    PR_R=[p for l,c,p in t1 if l=='R'][0]; PR_L=[p for l,c,p in t1 if l=='L'][0]
    d2={l:p for l,c,p in t2}
    for k,v in [("Nmax",N),("n_primes_gt3",n),("P_next_R_given_R","%.5f"%PR_R),
        ("P_next_R_given_RR","%.5f"%d2['RR']),("P_next_R_given_LR","%.5f"%d2['LR']),
        ("second_order_split_RR_vs_LR","%.5f"%(d2['RR']-d2['LR'])),
        ("step_lag1_autocorr","%.5f"%step_lag1),("gap_norm_lag1_corr","%.5f"%gap_corr),
        ("verdict","wing chain is NOT first-order Markov; effective second-order (lag>=3 step autocorr vanishes)")]:
        wcsv.writerow([k,v])
print("\nTOTAL %.1fs"%(time.time()-t0))
