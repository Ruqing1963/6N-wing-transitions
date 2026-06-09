import numpy as np, time
from math import isqrt, log
t0=time.time()
N=100_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
pr=primes[primes>3]                      # primes >3
print("primes>3 up to %d: %d  (%.1fs)"%(N,len(pr),time.time()-t0))

# wing symbol: R=1 (6N+1), L=0 (6N-1)
wing=(pr%6==1).astype(np.int8)           # 1=R, 0=L
nR=int(wing.sum()); nL=len(wing)-nR
print("R(6N+1)=%d  L(6N-1)=%d  R-L=%d  R/total=%.5f (Chebyshev: slight L excess)"%(nR,nL,nR-nL,nR/len(wing)))

# transition matrix on consecutive primes
a=wing[:-1]; b=wing[1:]
RR=int(((a==1)&(b==1)).sum()); RL=int(((a==1)&(b==0)).sum())
LR=int(((a==0)&(b==1)).sum()); LL=int(((a==0)&(b==0)).sum())
tot=RR+RL+LR+LL
print("\nTransition counts (consecutive primes):")
print("  R->R %d   R->L %d   L->R %d   L->L %d"%(RR,RL,LR,LL))
# conditional probabilities
print("\nConditional P(next | current)  [fair coin = 0.5]:")
print("  P(R->R)=%.5f  P(R->L)=%.5f"%(RR/(RR+RL),RL/(RR+RL)))
print("  P(L->L)=%.5f  P(L->R)=%.5f"%(LL/(LL+LR),LR/(LL+LR)))
psame=(RR+LL)/tot; pdiff=(RL+LR)/tot
print("  P(same wing)=%.5f   P(switch)=%.5f   [fair coin: 0.5/0.5]"%(psame,pdiff))
print("  => repetition bias: same - 0.5 = %+.5f"%(psame-0.5))

# run-length distribution vs geometric 2^-k
runs={}
cur=wing[0]; ln=1
for x in wing[1:]:
    if x==cur: ln+=1
    else: runs[ln]=runs.get(ln,0)+1; cur=x; ln=1
runs[ln]=runs.get(ln,0)+1
totruns=sum(runs.values())
print("\nRun-length distribution (consecutive same wing):")
print("  k   observed_frac   fair-coin(geometric)=2^-k*?  ratio")
# for fair coin, P(run=k)= (1/2)^k normalized over runs: actually P(run length=k)=2^{-k}
for k in range(1,9):
    obs=runs.get(k,0)/totruns
    geo=2.0**(-k)
    print("  %d   %.5f        %.5f          %.3f"%(k,obs,geo,obs/geo if geo>0 else 0))
print("  mean run length: obs=%.4f  fair-coin=2.0"%(len(wing)/totruns))

# bias decay: measure P(same) in windows to see trend with size
print("\nP(same wing) by decade window:")
edges=[10**k for k in range(4,9)]
for i in range(len(edges)-1):
    lo,hi=edges[i],edges[i+1]
    m=(pr[:-1]>=lo)&(pr[:-1]<hi)
    aa=a[m]; bb=b[m]
    if len(aa)>1000:
        ps=int(((aa==bb)).sum())/len(aa)
        print("  [%.0e,%.0e): P(same)=%.5f  (n=%d)"%(lo,hi,ps,len(aa)))
print("\nTOTAL %.1fs"%(time.time()-t0))
