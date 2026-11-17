"""crosstalk-qa: flag baselines with a coherent delay line beyond the foreground horizon.

    python -m crosstalk.cli OBS.uvfits|OBS.ms [...] [--metafits OBS.metafits] [--out prefix]

Exit status: 0 ok, 2 if any input was refused (bandwidth too narrow).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys

import numpy as np

from . import core
from .api import FIELDS, analyse
from .io import load, read_metafits_rx


def fmt_row(r):
    rx = f" Rx{r['rx1']}/{r['slot1']}-Rx{r['rx2']}/{r['slot2']}" if r["rx1"] is not None else ""
    return (
        f"{r['pol']:2} {r['ant1']:>8}-{r['ant2']:<8} {r['length_m']:6.1f}m hor {r['horizon_ns']:5.0f}ns "
        f"peak {r['tau_peak_ns']:+7.0f}ns T {r['t']:9.1f} z {r['z']:5.1f} coh {r['coherence']:.2f} "
        f"Tinc {r['t_inc']:7.1f} sep {r['sep']:.2f} Tother {r['t_otherpol']:7.1f} "
        f"1side {r['onesided_db']:3.0f}dB raw {r['tau_raw_ns']:+5.0f}ns {r['raw_db']:4.0f}dB "
        f"raw1side {r['raw_onesided_db']:3.0f}dB{rx}"
        f"{' FLAG' if r['flag'] else ''}"
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vis", nargs="+", help="uvfits or measurement set(s) (Birli output)")
    ap.add_argument("--metafits", help="for Rx/slot attribution (TILEDATA Rx, Slot); one for all inputs")
    ap.add_argument("--out", help="write <out>.csv (all searchable rows) and <out>.json")
    ap.add_argument("--pols", default="XX,YY")
    ap.add_argument("--no-unproject", action="store_true", help="use the stored phase centre as-is")
    ap.add_argument("--flag-edge-chans", type=int, default=0, help="extra fine chans to flag at each coarse edge")
    d = core.Params()
    for name in ("buffer_ns", "noise_min_ns", "t_min", "z_min", "coherence_min", "sep_min",
                 "max_resolution_ns"):
        ap.add_argument("--" + name.replace("_", "-"), type=float, default=getattr(d, name))
    ap.add_argument("--fg-extra-modes", type=int, default=d.fg_extra_modes)
    ap.add_argument("--pad", type=int, default=d.pad)
    ap.add_argument("--k-neighbours", type=int, default=d.k_neighbours)
    ap.add_argument("--coarse-bw-hz", type=float, default=d.coarse_bw_hz, help="0 disables harmonic masks")
    ap.add_argument("--onesided-min-db", type=float, default=d.onesided_min_db, help="-1000 disables")
    ap.add_argument("--min-chan-time-frac", type=float, default=d.min_chan_time_frac)
    ap.add_argument("--asym-min", type=float, default=None, help="require T/T_otherpol >= this")
    ap.add_argument("--top", type=int, default=15, help="also print the top-N by T")
    ap.add_argument("--force", action="store_true", help="run even if the bandwidth is too narrow")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    params = core.Params(
        buffer_ns=a.buffer_ns, noise_min_ns=a.noise_min_ns, t_min=a.t_min, z_min=a.z_min,
        coherence_min=a.coherence_min, sep_min=a.sep_min, max_resolution_ns=a.max_resolution_ns,
        fg_extra_modes=a.fg_extra_modes, pad=a.pad, k_neighbours=a.k_neighbours,
        coarse_bw_hz=a.coarse_bw_hz or None, asym_min=a.asym_min, force=a.force,
        onesided_min_db=a.onesided_min_db, min_chan_time_frac=a.min_chan_time_frac,
    )
    rx = read_metafits_rx(a.metafits) if a.metafits else None
    results, rc = {}, 0
    for path in a.vis:
        cube = load(path, unproject=not a.no_unproject, pols=a.pols.split(","), flag_edge_chans=a.flag_edge_chans)
        bw = np.median(np.diff(cube.freqs_hz)) * len(cube.freqs_hz)
        print(f"# {path}: {len(cube.length_m)} bl x {len(cube.times_jd)} t x {len(cube.freqs_hz)} ch "
              f"({bw / 1e6:.2f} MHz, delay res {1e9 / bw:.1f} ns), frame {cube.frame}, pols {cube.pols}")
        try:
            res = analyse(cube, params, rx)
        except core.BandwidthError as e:
            print(f"REFUSED: {e}", file=sys.stderr)
            rc = 2
            continue
        for w in res.warnings:
            print(f"WARNING: {w}")
        print(f"## flagged {len(res.flagged)} (baseline, pol)")
        for r in res.flagged:
            print(fmt_row(r))
        print(f"## top {a.top} by T")
        for r in sorted(res.rows, key=lambda r: -np.nan_to_num(r["t"]))[: a.top]:
            print(fmt_row(r))
        print("## implicated tiles (HUB = likely leak source)")
        for pol, tiles in res.tiles.items():
            for t in tiles:
                rxs = f" Rx{t['rx']}/{t['slot']}" if t["rx"] is not None else ""
                print(f"{pol} {t['tile']:>8}{rxs} n={t['n_flagged']} same_rx={t['n_same_rx']}"
                      f"{' HUB' if t['hub'] else ''} partners={','.join(t['partners'])}")
        results[path] = res
    if a.out:
        with open(a.out + ".json", "w") as f:
            json.dump({p: {k: v for k, v in r.__dict__.items() if k != "rows"} | {"flagged": r.flagged}
                       for p, r in results.items()}, f, indent=1, default=str)
        with open(a.out + ".csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["file"] + FIELDS)
            w.writeheader()
            for path, res in results.items():
                for r in res.rows:
                    w.writerow({"file": path, **r})
    return rc


if __name__ == "__main__":
    sys.exit(main())
