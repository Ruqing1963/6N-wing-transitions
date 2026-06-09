# 6N-dodge-ladder

**Part XXXIII — The Fermat-Type Dodge Ladder: Primes n^(2^k)+1 on the 6N Skeleton and the Non-Monotonic Richness**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXXIII of *Arithmetic Geodynamics on the 6N Skeleton*. We push the
"richness = dodging" principle into higher degree along `n²+1`, `n⁴+1`, `n⁸+1`. **Everything here is a
measured primality count — no fitted parameters, no fabricated numbers.**

## The cyclotomic dodge

`x^(2^k) ≡ −1 (mod q)` is solvable iff `q ≡ 1 (mod 2^(k+1))`, with exactly `2^k` roots. So `n^(2^k)+1`
is divisible **only** by `q = 2` and primes `q ≡ 1 (mod 2^(k+1))` — a vanishing fraction `2^(−k)` of
all primes — and where hit, it is hit on `2^k` residues at once. Local factor
`g(q) = (1 − ω(q)/q)/(1 − 1/q)`: dodged primes give `q/(q−1) > 1`, hit primes `(q − 2^k)/(q−1) < 1`.

## The two findings

| f | divides iff | ω | N_max | Q | right:left | C (meas. / prod.) |
|---|---|---|---|---|---|---|
| n²+1 | q≡1 (mod 4) | 2 | 10⁶ | 54,109 | 0.500 | 1.376 / 1.373 |
| n⁴+1 | q≡1 (mod 8) | 4 | 10⁶ | 52,610 | 0.499 | 2.676 / 2.679 |
| n⁸+1 | q≡1 (mod 16) | 8 | 10⁵ | 2,498 | 0.517 | 2.075 / 2.093 |

1. **The wing split is degree-invariant.** Since `n^(2^k) ≡ n² (mod 6)` for all k≥1, every `n^(2^k)+1`
   splits its primes right:left = **1:2** exactly as `n²+1` does (odd n annihilated, 6|n → right
   6N+1, n≡2,4 → left 6N−1). The exponent is invisible to the 6N lattice.

2. **The richness C is non-monotonic.** As the degree doubles the dodge widens (more enhancing primes,
   C up) but the hits deepen (ω = 2^k, heavier suppression, C down). The balance peaks at **n⁴+1**:
   C = 1.376 → **2.676** → 2.075 for k = 1, 2, 3. Confirmed by direct measurement, not just the
   product (the 22% drop n⁴→n⁸ far exceeds the ~2% sampling error on the smaller n⁸+1 sample).

**Scope (honest).** The solvability of `x^(2^k) ≡ −1` and its root count are classical cyclotomic
facts; C and the Bateman–Horn asymptotic are standard; the n²+1 constant is Landau–Shanks. What this
adds is the 6N reading of the ladder, the degree-invariance of the wings, and a direct high-precision
measurement of the richness (including the non-monotonic peak and a **degree-8** Bateman–Horn check,
where values reach 10⁴⁰). **The infinitude of primes n⁴+1 is open** (Landau-type); no claim is made.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_n4.py     # quick look: ω(q) ladder + wings (console)
python3 final_dodge.py    # n²+1, n⁴+1 (deterministic Miller-Rabin), n⁸+1 (sympy BPSW)
                          #   + analytic constants -> data/dodge_*.csv   (~20 s)
python3 makefigs_dodge.py # reads the CSVs -> figures/p33_fig1.pdf, p33_fig2.pdf
```

Enumeration is by direct primality testing (a root sieve would need primes to N^(2^k)); deterministic
Miller–Rabin (12 bases, valid below 3.3e24) for n²+1 and n⁴+1, sympy's BPSW for n⁸+1. Paths resolve
relative to the script. NumPy 2.x compliant.

## Files

```
code/    explore_n4.py   final_dodge.py   makefigs_dodge.py
data/    dodge_ladder.csv   degree, divisor_modulus, Nmax, Q, right, left, C_product, C_empirical
         dodge_omega.csv    q, q_mod16, omega_n2, omega_n4, omega_n8
         dodge_summary.csv  parameter, value
figures/ p33_fig1.pdf  p33_fig2.pdf
paper/   paper33.tex   paper33.pdf
```

All data files are plain CSV.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
