"""Independent brute force: for F_2^n with n = k+s, build I_M in simplicial order one vector at a time,
maintain the distance of every y to I_M, and test Theorem 3's condition |N_2t(I_M)| <= 2^n - 2^k + M
for every 1 <= M <= 2^(k-1). Prints whether redundancy s suffices for all f on F_2^k.
Does not use kappa.py."""
import sys
import numpy as np

def run(t, s, k):
    n = k + s; rho = 2 * t
    N = 1 << n
    ys = np.arange(N, dtype=np.int64)
    pc = np.array([bin(i).count("1") for i in range(N)], dtype=np.int16)
    # simplicial order: by weight, then lexicographically with coordinate 1 = most significant bit, 1 before 0
    order = sorted(range(N), key=lambda x: (pc[x], -x))
    dist = np.full(N, 99, dtype=np.int16)
    for M in range(1, (1 << (k - 1)) + 1):
        x = order[M - 1]
        np.minimum(dist, pc[ys ^ x], out=dist)
        if (dist <= rho).sum() > N - (1 << k) + M:
            return False, M
    return True, None

for spec in sys.argv[1:]:
    t, s, k = map(int, spec.split(","))
    okk, M = run(t, s, k)
    print(f"t={t} s={s} k={k}: {'all M ok' if okk else f'fails at M={M}'}")
