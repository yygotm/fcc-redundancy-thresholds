"""Redundancy thresholds kappa_s(t) of arXiv:2610.09011 (non-systematic FCCs, two-valued f).

Theorem 3 there: an (f,t)-FCC of length n exists iff |N_{2t}(I_M)| <= 2^n - 2^k + M,
I_M = first M vectors of F_2^n in the simplicial order, M = smaller preimage size (1 <= M <= 2^(k-1)).

Write M = V(n,j-1) + a with 1 <= a <= C(n,j). Then I_M = B_{j-1} u L_a, with L_a the first a
weight-j vectors (lex order of supports, element 1 first), and
|N_rho(I_M)| = V(n, j-1+rho) + U(n, j, a), U = size of the rho-th upper shadow of L_a.
U and max_a (U - a) satisfy the recursions below, so no loop over M is needed.

Usage: python kappa.py selftest | table T S KMAX | kappa T S
"""
import sys
from functools import lru_cache
from math import comb

sys.setrecursionlimit(100000)


def C(n, r):
    return comb(n, r) if 0 <= r <= n else 0


def V(n, r):
    if r < 0:
        return 0
    return sum(C(n, i) for i in range(min(r, n) + 1))


def make(rho):
    @lru_cache(maxsize=None)
    def U(n, m, a):
        if a == 0:
            return 0
        if m == 0:
            return C(n, rho)
        c = C(n - 1, m - 1)
        if a <= c:
            return U(n - 1, m - 1, a)
        return C(n - 1, m - 1 + rho) + U(n - 1, m, a - c)

    @lru_cache(maxsize=None)
    def G(n, m):  # max over 0 <= a <= C(n,m) of U(n,m,a) - a
        if m == 0:
            return max(0, C(n, rho) - 1)
        if m > n:
            return 0
        c = C(n - 1, m - 1)
        return max(G(n - 1, m - 1), C(n - 1, m - 1 + rho) - c + G(n - 1, m))

    def H(n, m, A):  # max over 0 <= a <= A of U(n,m,a) - a  (A <= C(n,m))
        best, base = 0, 0
        while True:
            if A == 0:
                return max(best, base)
            if m == 0:
                return max(best, base + C(n, rho) - 1)
            c = C(n - 1, m - 1)
            if A <= c:
                n, m = n - 1, m - 1
                continue
            best = max(best, base + G(n - 1, m - 1))
            base += C(n - 1, m - 1 + rho) - c
            n, A = n - 1, A - c

    return U, G, H


def max_excess(n, k, rho, fns):
    """max over 1 <= M <= 2^(k-1) of |N_rho(I_M)| - M, in F_2^n."""
    U, G, H = fns
    Mmax = 1 << (k - 1)
    best = V(n, rho) - 1  # M = 1 (j = 0)
    j = 1
    while j <= n and V(n, j - 1) < Mmax:
        lo = V(n, j - 1)
        A = min(C(n, j), Mmax - lo)
        best = max(best, V(n, j - 1 + rho) - lo + H(n, j, A))
        j += 1
    return best


def ok(k, s, t, fns):
    n = k + s
    return max_excess(n, k, 2 * t, fns) <= (1 << n) - (1 << k)


def kappa(t, s):
    """least kappa with ok(k) for all k >= kappa; uses Theorem 5's bound as the cutoff."""
    from math import ceil
    from fractions import Fraction
    bound = ceil(Fraction(2 * (2 * t + 1) ** 2) / (3 * (1 - Fraction(1, 2 ** s)) ** 2)) - s
    fns = make(2 * t)
    # scan downward from the bound: the first failing k is the largest one
    for k in range(max(bound, 1), 0, -1):
        if not ok(k, s, t, fns):
            return k + 1, bound
    return 1, bound


# ---------- brute force for checking ----------
def simplicial_order(n):
    vecs = list(range(1 << n))
    # bit i (from the left) = coordinate i+1; compare by weight, then x < y iff first differing coord is in supp(x)
    def key(x):
        return (bin(x).count("1"), [-((x >> (n - 1 - i)) & 1) for i in range(n)])
    return sorted(vecs, key=key)


def brute_excess_all(n, k, rho):
    order = simplicial_order(n)
    dist = lambda x, y: bin(x ^ y).count("1")
    res = {}
    S = []
    for M in range(1, (1 << (k - 1)) + 1):
        S.append(order[M - 1])
        N = sum(1 for y in range(1 << n) if any(dist(x, y) <= rho for x in S))
        res[M] = N - M
    return res


def selftest():
    for n in range(2, 9):
        for rho in range(1, 5):
            fns = make(rho)
            U = fns[0]
            order = simplicial_order(n)
            for M in range(1, 1 << n):
                S = order[:M]
                N = sum(1 for y in range(1 << n) if any(bin(x ^ y).count("1") <= rho for x in S))
                j = 0
                while V(n, j) < M:
                    j += 1
                if j == 0:
                    pred = V(n, rho)
                else:
                    pred = V(n, j - 1 + rho) + U(n, j, M - V(n, j - 1))
                assert pred == N, (n, rho, M, pred, N)
            for k in range(1, n + 1):
                bf = max(brute_excess_all(n, k, rho).values())
                assert bf == max_excess(n, k, rho, fns), (n, k, rho)
    print("selftest OK (n<=8, rho<=4: |N_rho(I_M)| for all M, and max excess)")


def table(t, s, kmax):
    rho = 2 * t
    fns = make(rho)
    U = fns[0]
    for k in range(1, kmax + 1):
        n = k + s
        good = []
        for M in range(1, (1 << (k - 1)) + 1):
            j = 0
            while V(n, j) < M:
                j += 1
            N = V(n, rho) if j == 0 else V(n, j - 1 + rho) + U(n, j, M - V(n, j - 1))
            if N <= (1 << n) - (1 << k) + M:
                good.append(M)
        print(k, compress(good), "ALL" if len(good) == 1 << (k - 1) else "")


def compress(xs):
    out, i = [], 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[j] + 1:
            j += 1
        out.append(f"{xs[i]}" if i == j else f"{xs[i]}-{xs[j]}")
        i = j + 1
    return ",".join(out) if out else "none"


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "selftest":
        selftest()
    elif cmd == "table":
        table(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "kappa":
        t, s = int(sys.argv[2]), int(sys.argv[3])
        kap, bound = kappa(t, s)
        print(f"kappa_{s}({t}) = {kap}   (checked k = 1..{bound}, Theorem 5 bound {bound})")
