# Redundancy thresholds κ_s(t) for non-systematic function-correcting codes

Exact values and the growth rate of the redundancy-s thresholds κ_s(t) defined in

> C. Rajput, K. L. Jagatti, B. Sundar Rajan, *On Optimal Encodings and Systematicity in Function-Correcting Codes*, arXiv:2610.09011v1 (2026), §III-A, Definition 1.

There, κ_1(1) = 10 is determined, and for t ≥ 2 only 4t − 1 ≤ κ_1(t) ≤ ⌈(8/3)(2t+1)²⌉ − 1 is known
(7 ≤ κ_1(2) ≤ 66); the growth of κ_1(t) (linear or quadratic) and lower bounds on κ_s(t) for s ≥ 2 are stated as open.

## Results

**Exact values** (computed from Theorem 3 of the paper; every k up to the bound of Theorem 5 is checked):

| t \ s | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 10 | | | | | | | | |
| 2 | **36** | **12** | 6 | | | | | | |
| 3 | 80 | 27 | 14 | 9 | 5 | | | | |
| 4 | 142 | 49 | 26 | 17 | 12 | 7 | 4 | | |
| 5 | 221 | 75 | 42 | 27 | 19 | 14 | 10 | 7 | 4 |

Also κ_1(6) = 318, κ_1(7) = 432, κ_1(8) = 564, κ_1(9) = 713, κ_1(10) = 880.

For each value, `results/witnesses-t1-5.txt` gives, at k = κ − 1, a class size M ≤ 2^{k−1} with
|N_{2t}(I_M)| − M > 2^{k+s} − 2^k (so redundancy s fails there). For example, for t = 2, s = 1, k = 35:
M = 15771159390, |N_4(I_M)| = 50460379332.

**Growth.** With z_s = Φ^{−1}(1 − 2^{−s−1}), for each fixed s ≥ 1,

  κ_s(t) / t² → 4 / z_s²  (t → ∞, s fixed),

e.g. κ_1(t) ~ 8.792 t², κ_2(t) ~ 3.023 t². So κ_1(t) grows quadratically, and for each fixed s ≥ 2,
κ_s(t) ≥ (4/z_s² − o(1)) t² as t → ∞ (an asymptotic lower bound, not uniform in s). A Berry–Esseen refinement gives
κ_s(t) = 4t²/z_s² + O_s(t). The proof (two elementary lemmas on Hamming balls plus the de Moivre–Laplace theorem) is in
[`asymptotics.md`](asymptotics.md). The lemmas show, with D_s(n) = n − 2 j_s(n) and j_s(n) the least j with
V(n, j) > 2^{n−s−1}: at k = n − s, redundancy s suffices if D_s(n) ≥ 2t + 2 and fails if D_s(n) ≤ 2t − 2.
Numerically κ_1(t) − 4t²/z_1² stays between 0.7 and 1.5 for t ≤ 10.

## Method

By Theorem 3 of the paper, redundancy s suffices for all nonconstant f on F_2^k iff
max_{1 ≤ M ≤ 2^{k−1}} (|N_{2t}(I_M)| − M) ≤ 2^{k+s} − 2^k. Writing M = V(n, j−1) + a, the set I_M is the ball B_{j−1}
plus the first a weight-j vectors in the simplicial order of the paper (lexicographic order of supports,
coordinate 1 first), and |N_ρ(I_M)| = V(n, j−1+ρ) + U(n, j, a), where U is the
size of the ρ-th upper shadow of that lex initial segment. U satisfies

  U(n, m, a) = U(n−1, m−1, a) if a ≤ C(n−1, m−1), else C(n−1, m−1+ρ) + U(n−1, m, a − C(n−1, m−1)),

and max_a (U − a) satisfies the matching recursion in (n, m), so the maximum over all M is computed exactly
(big integers) without enumerating M.

## Verification

- `python kappa.py selftest`: the formula for |N_ρ(I_M)| and the maximum agree with direct Hamming-distance
  computation for every M, n ≤ 8, ρ ≤ 4.
- `python kappa.py table 1 1 11` reproduces Table III of the paper (t = 1, k ≤ 9) and κ_1(1) = 10.
- `brute.py` (independent of `kappa.py`) builds I_M vector by vector and checks Theorem 3's condition directly;
  it confirms failure at k = κ − 1 and success at k = κ for 7 small cases (`results/brute-boundaries.txt`).
- An independent implementation written separately (GPT, without access to this code) obtained the same 25 values,
  κ_1(6), κ_1(7) and κ_1(8), and the same witness M for t = 2, s = 1.

Reproduce: `python kappa.py kappa T S` (Python 3, standard library only); `brute.py` needs numpy.
Run times: t ≤ 5 in seconds; κ_1(8) about 7 minutes, κ_1(9) about 23 minutes. For t ≥ 10 the recursion is deeper
than Python's default stack allows; `python3 run_big.py 10 1` runs it in a thread with a large stack (Linux;
κ_1(10) took about 50 minutes).

## Note on how this was produced

The computations and the proof were produced with AI assistance (Anthropic's Claude and OpenAI's GPT),
checked against each other and against brute force as described above. The proof in `asymptotics.md` was also
reviewed independently by GPT and by a separate Claude session. No human mathematician has checked it yet.

Hirotaka Shimizu
