#!/usr/bin/env python3
"""
Reproduce the computations recorded in docs/07_data.md.

  python3 scripts/experiments.py fields   --Dmax 200        # period, iota, kappa_sq, N(eps), Reg
  python3 scripts/experiments.py winner   19 22 31          # max rank-2 Gram on balanced class reps
  python3 scripts/experiments.py analyze  19:600 31:1400    # represented classes / lin-irred / n*
  python3 scripts/experiments.py frame57                    # imported Delta=57 certificate
  python3 scripts/experiments.py fan      66 146            # s=2 fan: winner + true Ind
  python3 scripts/experiments.py convergents 19 31 43 46    # joint rank-2 compatibility of odd convergents

All outputs are LOWER bounds for S(K,2) unless stated (see find_winner docstring).
"""
import os, sys, json, argparse, time
from math import isqrt, log
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rqf import QF, cf_quadratic, fundamental_unit, tp_unit
from indec import indecomposables, square_class_reps
from lattice import (balanced_classes, find_winner, realized_module, analyze_lattice)


def squarefree(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


def period(D):
    P, Q = (1, 2) if D % 4 == 1 else (0, 1)
    pre, per = cf_quadratic(P, Q, D)
    if not pre:                       # purely periodic (only D=5): [a0, a1..] = [a0; a1..,a0]
        return per[0], per[1:] + per[:1]
    return pre[0], per


def cmd_fields(a):
    print(f"{'D':>5} {'s':>3} {'iota':>5} {'kappa':>5} {'N(eps)':>6} {'Reg':>8}")
    for D in range(2, a.Dmax + 1):
        if not squarefree(D):
            continue
        u0, per = period(D)
        eps, n = fundamental_unit(D)
        x = eps.emb()[0]
        print(f"{D:>5} {len(per):>3} {len(indecomposables(D)):>5} {len(square_class_reps(D)):>5} "
              f"{n:>6} {log(x):>8.3f}")


def reps_for(D, kind):
    """'balanced' = one small-trace rep per square class; 'orbit' = one rep per O^{x,+}-orbit
    (normalized indecomposables(D)) -- the latter reproduces the 'orbit plateau' data."""
    return balanced_classes(D) if kind == "balanced" else indecomposables(D)


def cmd_winner(a):
    for D in map(int, a.D):
        t = time.time()
        M = find_winner(D, reps_for(D, a.reps))
        print(f"D={D}: winner size {len(M)}  diag={[str(M[i][i]) for i in range(len(M))]}  "
              f"({time.time()-t:.1f}s)")


def _report(D, r):
    cl = sorted(int(c.trace()) for c, _ in r["classes"])
    li = sorted(int(c.trace()) for c, _ in r["lin_irred_classes"])
    print(f"D={D} cap={r['cap']}: #classes={len(cl)} traces={cl} | #orbits={len(r['orbits'])} "
          f"| lin-irred classes={len(li)} {li} | n*={r['nstar']} | #short={r['n_short']}")


def cmd_analyze(a):
    for tok in a.spec:
        D, cap = map(int, tok.split(":"))
        M = find_winner(D, reps_for(D, a.reps))
        Gp, vs, basis = realized_module(D, M)
        _report(D, analyze_lattice(D, Gp, basis, cap))


def cmd_frame57(a):
    D = 57
    w, one, z = QF(D, F(1, 2), F(1, 2)), QF(D, 1), QF(D, 0)
    a1, a2, b = one * 10 + w * 3, one * 23 + w * 7, -(one * 13 + w * 4)
    delta = a1 * a2 - b * b
    print("delta =", delta, " = tp unit:", delta == tp_unit(D))
    basis = [(one, z), (w, z), (z, one), (z, w)]
    _report(D, analyze_lattice(D, [[a1, b], [b, a2]], basis, 1500))


def cmd_fan(a):
    for D in map(int, a.D):
        u0, per = period(D)
        assert len(per) == 2, "fan experiment is for period-2 fields"
        fan = [QF(D, 1 + r * u0, r) for r in range(per[0])]   # 1 + r*alpha0, alpha0 = u0+sqrt D
        M = find_winner(D, fan)
        Gp, vs, basis = realized_module(D, M)
        _report(D, analyze_lattice(D, Gp, basis, 30 * max(int(M[i][i].trace()) for i in range(len(M)))))


def cmd_convergents(a):
    for D in map(int, a.D):
        assert D % 4 != 1
        u0, per = period(D)
        s = len(per)
        us = [u0] + per
        p_, p, q_, q = 1, u0, 0, 1
        conv = [QF(D, p, q)]
        for i in range(1, s):
            p_, p = p, us[i] * p + p_
            q_, q = q, us[i] * q + q_
            conv.append(QF(D, p, q))
        odd = [conv[i] for i in range(s) if i % 2 == 1]
        M = find_winner(D, odd)
        print(f"D={D} s={s}: odd convergents traces {[int(x.trace()) for x in odd]} -> "
              f"max jointly rank-2 compatible {len(M)}/{len(odd)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("fields"); p.add_argument("--Dmax", type=int, default=100); p.set_defaults(f=cmd_fields)
    for name, fn in [("winner", cmd_winner), ("fan", cmd_fan), ("convergents", cmd_convergents)]:
        p = sp.add_parser(name); p.add_argument("D", nargs="+"); p.set_defaults(f=fn)
        p.add_argument("--reps", choices=["balanced", "orbit"], default="balanced")
    p = sp.add_parser("analyze"); p.add_argument("spec", nargs="+"); p.set_defaults(f=cmd_analyze)
    p.add_argument("--reps", choices=["balanced", "orbit"], default="balanced")
    p = sp.add_parser("frame57"); p.set_defaults(f=cmd_frame57)
    a = ap.parse_args()
    a.f(a)
