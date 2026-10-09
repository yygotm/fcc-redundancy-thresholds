# Growth of the redundancy thresholds κ_s(t)

Setting and notation of arXiv:2610.09011 (Rajput–Jagatti–Sundar Rajan). Fix integers s ≥ 1, t ≥ 1, put ρ = 2t, and
let n ≥ s + 1, k = n − s ≥ 1. Let W ~ Bin(n, 1/2), F_n(r) = P(W ≤ r) = V(n, r)/2^n, p = 2^{−s−1}. Define

  j_s(n) = min{ j ≥ 0 : V(n, j) > 2^{n−s−1} },   D_s(n) = n − 2 j_s(n).

Since p 2^n = 2^{k−1} ≥ 1, we have j_s(n) ≥ 1, and F_n(j_s(n) − 1) ≤ p < F_n(j_s(n)).

Theorem 3 of the paper: redundancy s suffices for every nonconstant f on F_2^k iff
for every 1 ≤ M ≤ 2^{k−1} = 2^{n−s−1}:  |N_{2t}(I_M)| − M ≤ 2^n − 2^k = 2^n − 2^{n−s}.   (★)

## Lemma A (failure from a Hamming ball)
Assume 2t ≥ D_s(n) + 2. Then (★) fails for M = V(n, j−1), j := j_s(n).

Proof. 1 ≤ M ≤ 2^{n−s−1} by the choice of j. I_M = B_{j−1}, so N_{2t}(I_M) = B_{j−1+2t}, and (★) reads
2^n F̄ + V(n, j−1) ≥ 2^{n−s}, where F̄ = P(W ≥ j + 2t). Since j + 2t ≥ n − j + 2,
P(W ≥ j+2t) ≤ P(W ≥ n−j+2) = P(W ≤ j−2) < P(W ≤ j−1) ≤ 2^{−s−1} (strict because C(n, j−1) > 0).
Hence 2^n F̄ + V(n, j−1) < 2^{n−s−1} + 2^{n−s−1} = 2^{n−s}, so (★) fails. ∎
(With 2t = D_s(n) + 1 the left side is exactly 2V(n, j−1), so (★) holds for this ball only if V(n, j−1) = 2^{n−s−1}.)

## Lemma B (success)
Assume 2t ≤ D_s(n) − 2. Then (★) holds for all 1 ≤ M ≤ 2^{n−s−1}.

Proof. Write j* = j_s(n); then 2j* + 2t + 2 ≤ n, so j* + 1 ≤ n/2.
Given M, let j ≥ 0 be least with V(n, j) ≥ M (so V(n, j−1) < M, with V(n, −1) = 0). As M ≤ 2^{n−s−1} < V(n, j*), j ≤ j*.
Since I_M ⊆ B_j, |N_{2t}(I_M)| ≤ V(n, j+2t), and |N_{2t}(I_M)| − M < V(n, j+2t) − V(n, j−1) = 2^n w(j),
w(j) := P(j ≤ W ≤ j + 2t).
w(j+1) − w(j) = 2^{−n}[C(n, j+2t+1) − C(n, j)] ≥ 0 when j + 2t + 1 ≤ n − j (by symmetry and unimodality of the
binomial coefficients), which holds for all j < j* because 2j* + 2t + 1 ≤ n. So w(j) ≤ w(j*), and it suffices that
w(j*) ≤ 1 − 2^{−s}, i.e.  P(W ≤ j*−1) + P(W ≥ j*+2t+1) ≥ 2^{−s}.
Put m = n − j* − 2t − 1 ≥ j* + 1. By symmetry P(W ≥ j*+2t+1) = P(W ≤ m) ≥ F_n(j*+1)
= F_n(j*−1) + 2^{−n}[C(n, j*) + C(n, j*+1)].
Also F_n(j*−1) + 2^{−n}C(n, j*) = F_n(j*) > 2^{−s−1}, and C(n, j*+1) ≥ C(n, j*) as j*+1 ≤ n/2. Therefore
F_n(j*−1) + P(W ≤ m) ≥ 2F_n(j*−1) + 2^{−n}[C(n,j*) + C(n,j*+1)] ≥ 2F_n(j*) > 2^{−s}. ∎

So, for k = n − s: redundancy s suffices if D_s(n) ≥ 2t + 2 and fails if D_s(n) ≤ 2t − 2.
Only three values D_s(n) ∈ {2t−1, 2t, 2t+1} are left to the exact computation.

## Theorem (growth)
Let z_s = Φ^{−1}(1 − 2^{−s−1}) (z_1 ≈ 0.6745, z_2 ≈ 1.1503, z_3 ≈ 1.5341). For each fixed integer s ≥ 1,

  κ_s(t) / t² → 4 / z_s²   as t → ∞   (4/z_1² ≈ 8.792, 4/z_2² ≈ 3.023, 4/z_3² ≈ 1.700).

In particular κ_1(t) grows quadratically, and for each fixed s ≥ 2 we get κ_s(t) ≥ (4/z_s² − o(1)) t² as t → ∞.
This lower bound is asymptotic, with no explicit range of t, and no uniformity in s is claimed.

Proof. Quantile limit. By the de Moivre–Laplace theorem, F_n(n/2 + x√n/2) → Φ(x) uniformly in x. Put q = −z_s =
Φ^{−1}(p). For fixed η > 0 and large n, the binomial CDF at the standardized points q − η and q + η is below and above
p respectively; since rounding to integers moves the standardized coordinate by O(n^{−1/2}), j_s(n) lies between them.
Hence (2 j_s(n) − n)/√n → −z_s, i.e. D_s(n) = (z_s + o(1))√n.

Lower bound. For 0 < ε < 1 and n = ⌊(1−ε) 4t²/z_s²⌋ (so n ≥ s + 1 for large t),
D_s(n) + 2 ≤ (1+o(1)) z_s √n + 2 ≤ (1−ε)^{1/2}(1+o(1)) 2t + 2 ≤ 2t for t large; by Lemma A redundancy s fails at
k = n − s, so κ_s(t) > n − s (Definition 1 requires success at every k ≥ κ_s(t)).

Upper bound. Fix ε > 0 and 0 < δ < 1 − (1+ε)^{−1/2}. Choose N with D_s(m) ≥ (1 − δ) z_s √m for all m ≥ N, and put
n_0 = ⌈(1+ε) 4t²/z_s²⌉. For t large, n_0 ≥ N and every m ≥ n_0 has D_s(m) ≥ (1−δ)(1+ε)^{1/2} 2t ≥ 2t + 2, so by Lemma B
redundancy s suffices at every k = m − s ≥ n_0 − s. Hence κ_s(t) ≤ n_0 − s ≤ (1+ε) 4t²/z_s² + 1. ∎

(D_s(n) is not monotone in n, e.g. D_1(4) = 2, D_1(5) = 1; the argument does not use monotonicity.)

## Remark (error term)
By Berry–Esseen, sup_r |F_n(r) − Φ((2r − n)/√n)| ≤ C_0/√n, and since the normal density is at most 1/√(2π), the same
holds with (2r + 1 − n)/√n for a larger constant C. Put e_n = C/√n, j = j_s(n), x_n = (2j − n)/√n. Applying this to
F_n(j−1) ≤ p < F_n(j) gives, for large n,

  Φ^{−1}(p − e_n) − 1/√n < x_n ≤ Φ^{−1}(p + e_n) + 1/√n.

Φ^{−1} is Lipschitz near p (the normal density at −z_s is positive), so x_n = −z_s + O_s(n^{−1/2}), i.e.
D_s(n) = z_s √n + O_s(1). If |D_s(n) − z_s √n| ≤ A_s for all large n, put
n_− = ⌊((2t − A_s − 2)/z_s)²⌋ and n_+ = ⌈((2t + A_s + 2)/z_s)²⌉. Then D_s(n_−) ≤ 2t − 2 (failure by Lemma A) and
D_s(n) ≥ 2t + 2 for all n ≥ n_+ (success by Lemma B), so n_− − s < κ_s(t) ≤ n_+ − s and

  κ_s(t) = 4t²/z_s² + O_s(t).

Numerically κ_1(t) − 4t²/z_1² = 1.2, 0.8, 0.9, 1.3, 1.2, 1.5, 1.2, 1.3, 0.8 for t = 1..9,
much smaller than the O(t) bound above.

## Exact values (computed; see kappa.py)
| t \ s | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 10 (paper) | | | | | | | | |
| 2 | 36 | 12 | 6 | | | | | | |
| 3 | 80 | 27 | 14 | 9 | 5 | | | | |
| 4 | 142 | 49 | 26 | 17 | 12 | 7 | 4 | | |
| 5 | 221 | 75 | 42 | 27 | 19 | 14 | 10 | 7 | 4 |
κ_1(6) = 318, κ_1(7) = 432, κ_1(8) = 564, κ_1(9) = 713.
