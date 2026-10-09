"""Run `kappa.py kappa T S` in a thread with a large stack (deep recursion for t >= 10).

Usage: python run_big.py T S
"""
import sys
import threading
import time

import kappa


def main():
    t, s = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    kap, bound = kappa.kappa(t, s)
    print(f"kappa_{s}({t}) = {kap}   (checked k = 1..{bound}, Theorem 5 bound {bound})  {time.time() - t0:.0f}s",
          flush=True)


sys.setrecursionlimit(1000000)
threading.stack_size(256 * 1024 * 1024)
th = threading.Thread(target=main)
th.start()
th.join()
