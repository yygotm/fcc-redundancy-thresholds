"""For given t, s, k: find M (smaller class size, 1 <= M <= 2^(k-1)) maximising |N_2t(I_M)| - M in F_2^(k+s),
and report whether Theorem 3's condition fails there."""
import sys
from functools import lru_cache
from kappa import make, V, C

t, s, k = map(int, sys.argv[1:4])
rho = 2 * t; n = k + s
U, G, H = make(rho)

@lru_cache(maxsize=None)
def Garg(n, m):  # (max_{0<=a<=C(n,m)} U-a, argmax a)
    if m == 0:
        return (C(n, rho) - 1, 1) if C(n, rho) - 1 > 0 else (0, 0)
    if m > n:
        return (0, 0)
    c = C(n - 1, m - 1)
    v1, a1 = Garg(n - 1, m - 1)
    v2, a2 = Garg(n - 1, m)
    v2 += C(n - 1, m - 1 + rho) - c
    return (v1, a1) if v1 >= v2 else (v2, c + a2)

def Harg(n, m, A):
    best, besta, base, off = 0, 0, 0, 0
    while True:
        if A == 0:
            return (base, off) if base > best else (best, besta)
        if m == 0:
            v = base + C(n, rho) - 1
            return (v, off + 1) if v > best else (best, besta)
        c = C(n - 1, m - 1)
        if A <= c:
            n, m = n - 1, m - 1
            continue
        v, a = Garg(n - 1, m - 1)
        if base + v > best: best, besta = base + v, off + a
        base += C(n - 1, m - 1 + rho) - c
        off += c; n, A = n - 1, A - c

Mmax = 1 << (k - 1); bound = (1 << n) - (1 << k)
worst = (V(n, rho) - 1, 1)
for j in range(1, n + 1):
    lo = V(n, j - 1)
    if lo >= Mmax: break
    v, a = Harg(n, j, min(C(n, j), Mmax - lo))
    if a == 0: continue
    val = V(n, j - 1 + rho) - lo + v
    if val > worst[0]: worst = (val, lo + a)
val, M = worst
# recompute directly from the formula for the found M
j = 0
while V(n, j) < M: j += 1
a = M - V(n, j - 1)
N = V(n, rho) if j == 0 else V(n, j - 1 + rho) + U(n, j, a)
assert N - M == val
print(f"t={t} s={s} k={k} n={n}: worst M={M} (= V(n,{j-1}) + {a}), |N_{rho}(I_M)|={N}, "
      f"|N|-M={val} {'>' if val > bound else '<='} 2^n-2^k={bound}  -> {'NO FCC of length k+s' if val > bound else 'all M ok'}")
assert val == max(V(n, rho) - 1, val)
