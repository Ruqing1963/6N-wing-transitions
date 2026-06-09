# 6N-wing-transitions

**Part XXXIV — Topological Strain Memory on the 6N Skeleton: the Effective Second-Order Markov Dynamics of Prime Wing Transitions**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXXIV of *Arithmetic Geodynamics on the 6N Skeleton*. We ask how deep
the *memory* of the prime wing string runs. **Everything here is a measured count over the primes — no
fitted parameters, no fabricated numbers.**

## The question

Each prime >3 sits on a wing: R = 6N+1 (≡1 mod 3) or L = 6N−1 (≡2 mod 3). The primes form a string
…RLRRLR…. Lemke Oliver & Soundararajan (2016) showed consecutive primes avoid repeating their class.
We ask the sharper question: **is the wing string first-order Markov, or does it remember further back —
and how far?** Measured over all 26,355,865 primes 3 < p ≤ 5×10⁸.

## Findings (all measured)

1. **Not first-order Markov.** Holding the current wing at R, the next wing still depends on the
   previous one: `P(R|RR) = 0.4252` vs `P(R|LR) = 0.4580` — a gap of 0.033 at **168σ**. A first-order
   chain requires these equal.

2. **Memory depth.** The switch (step) sequence `σ_i = 1[w_{i+1} ≠ w_i]` would be independent noise for
   a symmetric first-order chain. Instead its autocorrelation is `ρ₁ = −0.0328 (168σ)`, `ρ₂ = +0.0046
   (24σ)`, `ρ₃ ≈ 0`. A switch suppresses the next switch ("strain memory"). The incremental memory per
   added prime **halves per order** (0.033 → 0.016 → 0.009), descending to the finite-sample floor by
   order ~5 — so the chain is *effectively* (not exactly) second-order: dominant nearest-neighbour
   memory with a geometrically decaying tail.

3. **One effect, two projections.** The switch is *identically* the event "gap not divisible by 6"
   (agreement 1.0000), and the step lag-1 autocorrelation (−0.0328) equals the normalised-gap lag-1
   correlation (−0.0332). The wing memory and the known consecutive-gap anti-correlation are the same
   anti-clustering tendency at two resolutions.

4. **Slow fade with scale.** P(same) rises 0.404 → 0.446 from <10⁵ to ~5×10⁸; the second-order signal
   shrinks in step — the finite-size behaviour the Hardy–Littlewood framework requires.

**Scope (honest).** The anti-repetition bias is Lemke Oliver–Soundararajan; the heuristic that predicts
every consecutive-pattern frequency (hence all structure here) is Hardy–Littlewood. **We prove no
theorem and propose no mechanism.** What this adds is the explicit memory-depth characterisation (order,
sign, decay) in 6N wing language, and the identification of the wing memory with the gap
anti-correlation. "Strain"/"anti-clustering" name measured correlations, not a physical force.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_wings.py   # wing transitions + run lengths (console)
python3 explore_gaps.py    # consecutive-gap correlations (console)
python3 explore_markov.py  # Markov-order probe (console)
python3 final_markov.py    # full measurement to N=5e8 -> data/markov_*.csv   (~17 s, ~0.7 GB RAM)
python3 makefigs_markov.py # reads the CSVs -> figures/p34_fig1.pdf, p34_fig2.pdf
```

Paths resolve relative to the script. NumPy 2.x compliant; single-threaded.

## Files

```
code/    explore_wings.py  explore_gaps.py  explore_markov.py  final_markov.py  makefigs_markov.py
data/    markov_transitions.csv   order, pattern, count, P(next=R)   [orders 1-4]
         markov_autocorr.csv      lag, step_autocorr, sigma          [lags 1-8]
         markov_memory_depth.csv  order, incremental_memory, ratio
         markov_decay.csv         window, P_same, second_order_signal
         markov_unification.csv   step vs gap lag-1; step==[gap%6!=0]
         markov_gapcond.csv       E[g_{n+1} | g_n bucket]
         markov_summary.csv       parameter, value
figures/ p34_fig1.pdf  p34_fig2.pdf
paper/   paper34.tex   paper34.pdf
```

All data files are plain CSV.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
