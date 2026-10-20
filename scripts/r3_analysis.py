#!/usr/bin/env python3
"""Tables of docs/11_R3.md sections 4 and 6 from data/r3_scan_1e6.jsonl.gz and data/r3_top_1e7.jsonl.
  avg := 2 * sum_f nu(f) / h(-8p)   (average of n(c) over A = support group, Theorem 11.1(b))
  x   := iota / sqrt p
  L8  := pi h(-8p) / sqrt(8p) = L(1, chi_{-8p})
Usage: python3 scripts/r3_analysis.py"""
import os, sys, json, gzip, math, statistics as st
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def load(name):
    p = os.path.join(DATA, name)
    op = gzip.open if p.endswith(".gz") else open
    return [json.loads(l) for l in op(p, "rt")]


def enrich(r):
    p = r["D"]
    r["x"] = r["iota"] / math.sqrt(p)
    r["avg"] = 2 * r["nu_sum"] / r["h8"]
    r["L8"] = math.pi * r["h8"] / math.sqrt(8 * p)
    return r


def corr(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return num / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))


def table(rs, label):
    print(f"\n{label}: range | #p | mean nu* | max nu* | mean avg | mean x | mean avg/x | mean nu*/avg"
          f" | [h(K)=1: #, mean nu*, mean x]")
    lo = 1000
    while lo < 2 * max(r["D"] for r in rs):
        b = [r for r in rs if lo <= r["D"] < 10 * lo]
        if b:
            h1 = [r for r in b if r["hK"] == 1]
            ok = [r for r in b if r["x"] > 0.5]
            print(f"[{lo:.0e},{10 * lo:.0e}) | {len(b)} | {st.mean(r['nu_max'] for r in b):.2f} | "
                  f"{max(r['nu_max'] for r in b)} | {st.mean(r['avg'] for r in b):.2f} | "
                  f"{st.mean(r['x'] for r in b):.2f} | {st.mean(r['avg'] / r['x'] for r in ok):.2f} | "
                  f"{st.mean(r['nu_max'] / r['avg'] for r in b if r['avg'] > 0):.2f} | "
                  f"{len(h1)}, {st.mean(r['nu_max'] for r in h1) if h1 else 0:.2f}, "
                  f"{st.mean(r['x'] for r in h1) if h1 else 0:.2f}")
        lo *= 10


def main():
    rs = [enrich(r) for r in load("r3_scan_1e6.jsonl.gz")]
    rs.sort(key=lambda r: r["D"])
    print(f"{len(rs)} primes p = 3 mod 4 in [{rs[0]['D']}, {rs[-1]['D']}]")
    # records
    best, rec = 0, []
    for r in rs:
        if r["nu_max"] > best:
            best = r["nu_max"]
            rec.append((r["D"], r["nu_max"], r["s"], r["iota"], r["h8"], r["hK"], round(r["x"], 2)))
    print("records (p, nu*, s, iota, h(-8p), h(K), iota/sqrt p):")
    for t in rec:
        print("  ", t)
    table(rs, "p <= 10^6")
    # identities (Zagier, Hirzebruch) on h(K) = 1
    h1 = [r for r in rs if r["hK"] == 1 and r["D"] > 3]
    print(f"\nh(K)=1 primes: {len(h1)} of {len(rs)};  Zagier sum a_i = T/2: "
          f"{sum(r['sigma_odd'] + r['sigma_even'] == r['T'] // 2 for r in h1)};  Hirzebruch "
          f"even - odd = 3h(-p): {sum(r['sigma_even'] - r['sigma_odd'] == 3 * r['hm'] for r in h1)}")
    print(f"h(K) odd for all: {all(r['hK'] % 2 for r in rs)}")
    # the cancellation: avg / x against L(1, chi_{-8p}), and avg against x
    big = [r for r in rs if r["D"] > 10 ** 5 and r["x"] > 0.5]
    y = [r["avg"] / r["x"] for r in big]
    print(f"\np in (1e5, 1e6], x > 0.5 ({len(big)}): corr(avg, x) = {corr([r['avg'] for r in big], [r['x'] for r in big]):.3f};"
          f" corr(avg/x, L(1,chi_-8p)) = {corr(y, [r['L8'] for r in big]):.3f};"
          f" corr(avg, L(1,chi_-8p)) = {corr([r['avg'] for r in big], [r['L8'] for r in big]):.3f};"
          f" avg/x: mean {st.mean(y):.3f}, sd {st.pstdev(y):.3f}")
    print(f"   corr(nu*, x) = {corr([r['nu_max'] for r in big], [r['x'] for r in big]):.3f};"
          f" corr(2 N_E/h8, x) = {corr([2 * r['nE'] / r['h8'] for r in big], [r['x'] for r in big]):.3f}")
    try:
        top = [enrich(r) for r in load("r3_top_1e7.jsonl")]
    except FileNotFoundError:
        return
    sel = [r for r in top if r.get("sample") == "top"]
    rnd = [r for r in top if r.get("sample") == "random"]
    print(f"\n(1e6, 1e7]: top-by-x {len(sel)}, random {len(rnd)}")
    for lab, b in (("top", sel), ("random", rnd)):
        if b:
            print(f"  {lab}: mean nu* {st.mean(r['nu_max'] for r in b):.2f}, max nu* {max(r['nu_max'] for r in b)},"
                  f" mean avg {st.mean(r['avg'] for r in b):.2f}, mean x {st.mean(r['x'] for r in b):.2f},"
                  f" mean avg/x {st.mean(r['avg'] / r['x'] for r in b if r['x'] > 0.5):.2f},"
                  f" #h(K)=1 {sum(r['hK'] == 1 for r in b)}")
    print("  best:", [(r["D"], r["nu_max"], r["s"], round(r["x"], 2), r["hK"]) for r in
                      sorted(top, key=lambda r: -r["nu_max"])[:12]])


if __name__ == "__main__":
    main()
