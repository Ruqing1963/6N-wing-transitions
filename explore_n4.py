import time
from math import log
t0=time.time()

# deterministic Miller-Rabin (first 12 primes -> correct for n < 3.317e24)
_BASES=[2,3,5,7,11,13,17,19,23,29,31,37]
def is_prime(n):
    if n<2: return False
    for p in _BASES:
        if n%p==0: return n==p
    d=n-1; r=0
    while d%2==0: d//=2; r+=1
    for a in _BASES:
        x=pow(a,d,n)
        if x==1 or x==n-1: continue
        for _ in range(r-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

# omega(q): number of solutions of n^(2^k)+1 == 0 mod q  (brute, small q)
def omega(q,deg):
    return sum(1 for n in range(q) if (pow(n,deg,q)+1)%q==0)

for deg,modulus in [(2,4),(4,8),(8,16)]:
    print("\n=== n^%d+1 : expect divisible only by q=2 or q=1 mod %d, omega=%d ==="%(deg,modulus,deg))
    bad=[]
    for q in [2,3,5,7,11,13,17,19,23,29,31,37,41,73,89,97,113]:
        om=omega(q,deg)
        if om>0: bad.append((q,q%modulus,om))
    print("  primes q with omega>0:", bad)

# wings & count for n^4+1
def analyze(deg,Nmax):
    Q=0; right=0; left=0; other=0
    for n in range(2,Nmax+1):
        if n%2: continue                 # n odd -> n^deg+1 even
        v=n**deg+1
        if is_prime(v):
            Q+=1; m=v%6
            if m==1: right+=1
            elif m==5: left+=1
            else: other+=1
    return Q,right,left,other

for deg in (2,4):
    Nmax=200000
    Q,right,left,other=analyze(deg,Nmax)
    rl=right/left if left else 0
    print("\nn^%d+1, n<=%d: Q=%d right=%d left=%d r/l=%.4f other=%d (%.1fs)"%(
        deg,Nmax,Q,right,left,rl,other,time.time()-t0))

print("\nTOTAL %.1fs"%(time.time()-t0))
