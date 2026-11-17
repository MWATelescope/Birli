"""Scan every finished Birli MS of the night and summarise.

    python -m crosstalk.night [--ms-dir ../ms] [--raw-dir ../raw] [--out-dir ../xt-out] [--watch]

Per obs: <out-dir>/<obsid>.csv/.json (as the CLI) and a line in <out-dir>/night.tsv.
--watch keeps polling for new *.ms (ignores *.ms.part) until --until-count are done.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import time

from astropy.io import fits

from .api import FIELDS, scan

WATCH = [("HexE11", "HexE12"), ("HexE11", "HexE13"), ("HexE11", "HexE14"), ("HexE10", "HexE11"),
         ("HexE12", "HexE13")]


def one(ms, raw_dir, out_dir):
    obsid = os.path.basename(ms).split(".")[0]
    mf = os.path.join(raw_dir, obsid, f"{obsid}.metafits")
    mf = mf if os.path.exists(mf) else None
    hdr = fits.getheader(mf) if mf else {}
    t0 = time.time()
    res = scan(ms, metafits=mf)
    with open(os.path.join(out_dir, f"{obsid}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(res.rows)
    with open(os.path.join(out_dir, f"{obsid}.json"), "w") as f:
        json.dump({"warnings": res.warnings, "tiles": res.tiles, "flagged": res.flagged, "params": res.params,
                   "frame": res.frame, "gridnum": hdr.get("GRIDNUM"), "lst": hdr.get("LST")}, f, indent=1,
                  default=str)
    idx = {(r["pol"], r["ant1"], r["ant2"]): r for r in res.rows}
    watch = []
    for a, b in WATCH:
        r = idx.get(("XX", a, b)) or idx.get(("XX", b, a))
        watch.append(f"{r['t']:.0f}{'*' if r['flag'] else ''}@{r['tau_peak_ns']:+.0f}" if r else "-")
    hubs = {pol: ",".join(f"{t['tile']}({t['n_flagged']})" for t in tl if t["hub"])[:200]
            for pol, tl in res.tiles.items()}
    nflag = {pol: sum(1 for r in res.flagged if r["pol"] == pol) for pol in ("XX", "YY")}
    return [obsid, hdr.get("GRIDNUM"), round(hdr.get("LST", float("nan")), 2), nflag["XX"], nflag["YY"],
            *watch, hubs.get("XX", ""), hubs.get("YY", ""), round(time.time() - t0)]


def main(argv=None):
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--ms-dir", default=os.path.join(here, "ms"))
    ap.add_argument("--raw-dir", default=os.path.join(here, "raw"))
    ap.add_argument("--out-dir", default=os.path.join(here, "xt-out"))
    ap.add_argument("--watch", action="store_true")
    ap.add_argument("--until-count", type=int, default=71)
    ap.add_argument("--only", nargs="*", help="obsids to (re)do")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    tsv = os.path.join(a.out_dir, "night.tsv")
    if not os.path.exists(tsv):
        with open(tsv, "w") as f:
            f.write("\t".join(["obsid", "grid", "lst", "nflagXX", "nflagYY"]
                              + [f"{x}-{y}" for x, y in WATCH] + ["hubsXX", "hubsYY", "secs"]) + "\n")
    done = set()
    while True:
        todo = sorted(glob.glob(os.path.join(a.ms_dir, "*.ms")))
        if a.only:
            todo = [m for m in todo if os.path.basename(m).split(".")[0] in a.only]
        new = [m for m in todo if m not in done and (a.only or not os.path.exists(
            os.path.join(a.out_dir, os.path.basename(m).split(".")[0] + ".json")))]
        for m in new:
            try:
                row = one(m, a.raw_dir, a.out_dir)
            except Exception as e:  # keep the night going
                row = [os.path.basename(m), "ERROR", repr(e)[:200]]
            with open(tsv, "a") as f:
                f.write("\t".join(map(str, row)) + "\n")
            print("\t".join(map(str, row)), flush=True)
            done.add(m)
        done.update(todo)
        if not a.watch or len(glob.glob(os.path.join(a.out_dir, "*.json"))) >= a.until_count:
            break
        time.sleep(60)


if __name__ == "__main__":
    main()
