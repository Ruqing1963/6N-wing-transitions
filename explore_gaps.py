import numpy as np, time
from math import isqrt, log
t0=time.time()
N=100_000_000
s=np.ones(N+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(N)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
pr=primes[primes>3].astype(np.float64)
d=np.diff(pr)                              # gaps d_n = p_{n+1}-p_n
print("gaps: %d  mean=%.3f  (%.1fs)"%(len(d),d.mean(),time.time()-t0))

# raw correlation of consecutive gaps
x=d[:-1]; y=d[1:]
c_raw=np.corrcoef(x,y)[0,1]
print("\nraw corr(d_n, d_{n+1}) = %+.5f"%c_raw)

# normalized gaps g = d/log(p) (Cramer: mean ~1, exponential, independent)
g=d/np.log(pr[:-1])
xg=g[:-1]; yg=g[1:]
c_norm=np.corrcoef(xg,yg)[0,1]
print("normalized corr(g_n, g_{n+1}) = %+.5f   [Cramer independent: 0]"%c_norm)
print("normalized gap mean=%.4f  var=%.4f  [exp: mean1 var1]"%(g.mean(),g.var()))

# E[g_{n+1} | g_n bucket]  -- does a big gap predict the next?
print("\nE[g_{n+1} | g_n in bucket]   (fair/independent: flat ~%.3f):"%g.mean())
buckets=[(0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.0,3.0),(3.0,10)]
for lo,hi in buckets:
    m=(xg>=lo)&(xg<hi)
    if m.sum()>500:
        print("  g_n in [%.1f,%.1f): E[g_{n+1}]=%.4f  (n=%d)"%(lo,hi,yg[m].mean(),int(m.sum())))

# small-gap conditioning (raw): after gap d_n=2 (twin), what is next gap?
print("\nE[d_{n+1} | d_n = g] for small gaps  (overall mean d=%.3f):"%d.mean())
for gg in [2,4,6,8,10,12,14]:
    m=(x==gg)
    if m.sum()>500:
        print("  d_n=%2d : E[d_{n+1}]=%.3f  P(d_{n+1}=d_n)=%.4f  (n=%d)"%(
            gg,y[m].mean(),float((y[m]==gg).sum())/m.sum(),int(m.sum())))

# gap mod 6 transition (connects to wing chain): gaps are even, mod6 in {0,2,4}
print("\ngap mod 6 distribution and 'after gap≡0' vs others:")
dm=(d.astype(np.int64))%6
for r in (0,2,4):
    print("  d≡%d mod6: frac=%.4f"%(r,(dm==r).mean()))
# does d_n≡0 mod6 (same wing) cluster? P(d_{n+1}≡0 | d_n≡0) vs P(d_{n+1}≡0)
dm0=dm[:-1]; dm1=dm[1:]
base=(dm1==0).mean()
cond=(dm1[dm0==0]==0).mean()
print("  P(d_{n+1}≡0 mod6)=%.4f ;  P(d_{n+1}≡0 | d_n≡0)=%.4f  (excess %+.4f)"%(base,cond,cond-base))
print("\nTOTAL %.1fs"%(time.time()-t0))
