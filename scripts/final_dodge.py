import numpy as np, csv, os, time
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
from math import log, isqrt
import sympy
t0=time.time()

_B=[2,3,5,7,11,13,17,19,23,29,31,37]
def is_prime_det(n):                 # deterministic for n < 3.317e24
    if n<2: return False
    for p in _B:
        if n%p==0: return n==p
    d=n-1; r=0
    while d%2==0: d//=2; r+=1
    for a in _B:
        x=pow(a,d,n)
        if x==1 or x==n-1: continue
        for _ in range(r-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

def count(deg,Nmax,prime_fn):
    Q=right=left=other=0
    for n in range(2,Nmax+1,2):       # n even (odd -> value even)
        v=n**deg+1
        if prime_fn(v):
            Q+=1; m=v%6
            if m==1: right+=1
            elif m==5: left+=1
            else: other+=1
    # Bateman-Horn integral  int dt/log(t^deg+1) over even n only? No: over all n, p=2 factor handles parity.
    tt=np.arange(2,Nmax+1,dtype=np.float64); val=tt**deg+1.0
    BH=float(np.sum(1.0/np.log(val)))
    return Q,right,left,other,BH

# primes to 1e7 for the constant product
P=10_000_000
s=np.ones(P+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(P)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
oddp=primes[primes>2]
print("primes to %d ready (%.1fs)"%(P,time.time()-t0))

def Cproduct(deg):
    mod=2*deg                         # q=1 mod 2*deg can divide (omega=deg)
    g=np.where(oddp%mod==1,(oddp-deg)/(oddp-1),oddp/(oddp-1))
    return float(np.prod(g))

# ladder rows
configs=[(2,1_000_000,is_prime_det),(4,1_000_000,is_prime_det),(8,100_000,sympy.isprime)]
rows=[]
print("\n deg  modulus  Nmax     Q      right   left   r/l    C_prod  C_emp")
for deg,Nmax,fn in configs:
    Q,right,left,other,BH=count(deg,Nmax,fn)
    Cp=Cproduct(deg); Ce=Q/BH
    rl=right/left if left else 0
    rows.append((deg,2*deg,Nmax,Q,right,left,other,Cp,Ce))
    print(" %2d    %3d   %8d %7d %6d %6d %.4f  %.4f  %.4f"%(deg,2*deg,Nmax,Q,right,left,rl,Cp,Ce))
    print("    (%.1fs)"%(time.time()-t0))

with open(os.path.join(_DATA,"dodge_ladder.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["degree","divisor_modulus","Nmax","Q","right_6Np1","left_6Nm1","other","C_product","C_empirical"])
    for r in rows: w.writerow([r[0],r[1],r[2],r[3],r[4],r[5],r[6],"%.5f"%r[7],"%.5f"%r[8]])

# omega(q) ladder over small primes
def omega(q,deg): return sum(1 for n in range(q) if (pow(n,deg,q)+1)%q==0)
with open(os.path.join(_DATA,"dodge_omega.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["q","q_mod16","omega_n2","omega_n4","omega_n8"])
    for q in [2,3,5,7,11,13,17,29,37,41,73,89,97,113,193,257]:
        w.writerow([q,q%16,omega(q,2),omega(q,4),omega(q,8)])

with open(os.path.join(_DATA,"dodge_summary.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["parameter","value"])
    for r in rows:
        w.writerow(["n^%d+1_C_empirical"%r[0],"%.4f"%r[8]])
        w.writerow(["n^%d+1_right_over_left"%r[0],"%.4f"%(r[4]/r[5] if r[5] else 0)])
    w.writerow(["dodge_rule","n^(2^k)+1 divisible only by q=2 or q=1 mod 2^(k+1), omega=2^k"])
print("\nTOTAL %.1fs"%(time.time()-t0))
