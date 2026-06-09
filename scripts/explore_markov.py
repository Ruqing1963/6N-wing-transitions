import numpy as np, time
from math import isqrt
t0=time.time()
N=200_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
pr=primes[primes>3]
w=(pr%6==1).astype(np.int8)              # 1=R, 0=L
n=len(w)
print("primes>3 up to %d: %d  (%.1fs)"%(N,n,time.time()-t0))

# --- 1st order (symmetry check) ---
a=w[:-1]; b=w[1:]
PRR=((a==1)&(b==1)).mean()/(a==1).mean()
PLL=((a==0)&(b==0)).mean()/(a==0).mean()
q=((a!=b)).mean()                         # overall switch prob
print("\n1st order: P(R|R)=%.5f  P(L|L)=%.5f  (symmetry ok if equal)"%(PRR,PLL))
print("overall P(switch q)=%.5f  P(same)=%.5f"%(q,1-q))

# --- step (switch) sequence: step_i = 1 if wing changes ---
step=(w[1:]!=w[:-1]).astype(np.int8)      # length n-1
m=len(step)
ps=step.mean()
# autocorrelation: P(step_{i+1}=1 | step_i=value)
s0=step[:-1]; s1=step[1:]
PW_afterW=s1[s0==1].mean()                # P(switch | prev was switch)
PW_afterS=s1[s0==0].mean()                # P(switch | prev was same)
print("\n--- Markov-order test via step (switch) sequence ---")
print("P(switch)=%.5f"%ps)
print("P(switch | prev step = switch)=%.5f"%PW_afterW)
print("P(switch | prev step = same  )=%.5f"%PW_afterS)
print("=> 2nd-order signal: P(W|W)-P(W|S) = %+.5f   (Markov-1 symmetric => 0)"%(PW_afterW-PW_afterS))
# significance
import math
se=math.sqrt(ps*(1-ps)*(1/(s0==1).sum()+1/(s0==0).sum()))
print("   approx SE of that difference ~ %.5f  => %.0f sigma"%(se,abs(PW_afterW-PW_afterS)/se))
# lag-2, lag-3 autocorr of step sequence
sm=step.astype(np.float64)-ps
for lag in (1,2,3,4):
    ac=np.sum(sm[:-lag]*sm[lag:])/np.sum(sm*sm)
    print("   step autocorr lag %d = %+.5f"%(lag,ac))

# --- triples: P(z | x,y) vs Markov-1 prediction P(z|y) ---
print("\n--- triples P(next=R | prev,cur) ; Markov-1: depends only on cur ---")
x=w[:-2]; y=w[1:-1]; z=w[2:]
def cond(px,py):
    mm=(x==px)&(y==py)
    return z[mm].mean(), int(mm.sum())
for px in (1,0):
    for py in (1,0):
        pr_,cnt=cond(px,py)
        lab={1:'R',0:'L'}
        print("   P(R | %s%s)=%.5f  (n=%d)"%(lab[px],lab[py],pr_,cnt))
print("   [Markov-1: P(R|RR)=P(R|LR) and P(R|RL)=P(R|LL)]")

# --- decay of 2nd-order signal by decade ---
print("\n2nd-order signal P(W|W)-P(W|S) by window:")
edges=[10**k for k in range(4,9)]+[2*10**8]
base=pr[1:-1]   # align with step pairs (step index i corresponds to gap between w[i],w[i+1])
pp=pr[1:len(step)]   # position proxy
for i in range(len(edges)-1):
    lo,hi=edges[i],edges[i+1]
    mm=(pr[1:len(step)]>=lo)&(pr[1:len(step)]<hi)
    s0w=step[:-1][mm[:-1]] if False else None
    # simpler: recompute on masked step pairs
    idx=np.nonzero((pr[1:-1]>=lo)&(pr[1:-1]<hi))[0]
    if len(idx)>5000:
        s0m=step[idx]; s1m=step[idx+1]
        ww=s1m[s0m==1].mean(); ws=s1m[s0m==0].mean()
        print("   [%.0e,%.0e): P(W|W)-P(W|S)=%+.5f  (n=%d)"%(lo,hi,ww-ws,len(idx)))
print("\nTOTAL %.1fs"%(time.time()-t0))
