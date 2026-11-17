"""High-level entry point: scan a UVData / uvfits / MS for crosstalk delay lines.

    from crosstalk import scan
    res = scan("1472673248.ms", metafits="1472673248.metafits")
    res.flagged          # list of dict rows (pol, ant1, ant2, tau_peak_ns, t, z, ...)
    res.tiles["XX"]      # implicated-tile summary
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from . import core
from .io import Cube, load, read_metafits_rx

OTHER_POL = {"XX": "YY", "YY": "XX", "XY": "YX", "YX": "XY"}

FIELDS = [
    "pol", "ant1", "ant2", "length_m", "horizon_ns", "tau_peak_ns", "t", "z", "coherence",
    "t_inc", "sep", "t_otherpol", "onesided_db", "tau_raw_ns", "raw_db", "raw_onesided_db", "good_frac", "same_rx", "rx1", "slot1", "rx2", "slot2", "flag",
]


@dataclass
class Result:
    rows: list  # one dict per searchable (baseline, pol), keys = FIELDS
    warnings: list
    tiles: dict  # pol -> implicated-tile list
    frame: str = ""
    n_times: int = 0
    bandwidth_hz: float = 0.0
    params: dict = field(default_factory=dict)

    @property
    def flagged(self):
        return sorted((r for r in self.rows if r["flag"]), key=lambda r: -r["t"])


def analyse(cube: Cube, params: core.Params | None = None, rx: dict | None = None) -> Result:
    params = params or core.Params()
    warnings = core.check_bandwidth(cube.freqs_hz, params)
    nt = len(cube.times_jd)
    if nt < 4:
        warnings.append(f"only {nt} timesteps: coherence requirement skipped")
    horizon = cube.length_m / core.C_M_PER_NS
    out = {}
    for ip, pol in enumerate(cube.pols):
        out[pol] = core.run_pol(cube.vis[..., ip], cube.flags[..., ip], cube.freqs_hz, cube.length_m, params)
    rows = []
    none3 = (None, None, None)
    for pol, (s, spec, tau) in out.items():
        o = out.get(OTHER_POL.get(pol))
        t_other = core.value_at(o[1], o[2], s.tau_peak_ns) if o is not None else None
        flag = core.decide(s, nt, params, t_other)
        for i in np.flatnonzero(s.searchable):
            r1 = rx.get(cube.name1[i], none3) if rx else none3
            r2 = rx.get(cube.name2[i], none3) if rx else none3
            rows.append(dict(
                pol=pol, ant1=cube.name1[i], ant2=cube.name2[i],
                length_m=round(float(cube.length_m[i]), 2), horizon_ns=round(float(horizon[i]), 1),
                tau_peak_ns=round(float(s.tau_peak_ns[i]), 1), t=float(s.t[i]), z=float(s.z[i]),
                coherence=float(s.coherence[i]), t_inc=float(s.t_inc[i]), sep=float(s.sep[i]),
                t_otherpol=float(t_other[i]) if t_other is not None else float("nan"),
                tau_raw_ns=round(float(s.tau_raw_ns[i]), 1), raw_db=round(float(s.raw_db[i]), 1),
                onesided_db=round(float(s.onesided_db[i]), 1),
                raw_onesided_db=round(float(s.raw_onesided_db[i]), 1),
                good_frac=round(float(s.good_frac[i]), 3),
                same_rx=bool(r1[1] is not None and r1[1] == r2[1]),
                rx1=r1[1], slot1=r1[2], rx2=r2[1], slot2=r2[2], flag=bool(flag[i]),
            ))
    n_un = {p: int((~v[0].searchable).sum()) for p, v in out.items()}
    if any(n_un.values()):
        warnings.append(
            f"unsearchable baselines per pol {n_un} (of {len(cube.length_m)}): horizon+buffer plus "
            "coarse-harmonic masks cover all delays (|b| > ~100 m), or too few good channels")
    df = float(np.median(np.diff(cube.freqs_hz)))
    return Result(rows, warnings, tile_summary(rows), cube.frame, nt, df * len(cube.freqs_hz),
                  dict(params.__dict__))


def scan(source, params: core.Params | None = None, metafits: str | None = None, pols=("XX", "YY"),
         unproject: bool = True, flag_edge_chans: int = 0) -> Result:
    """source: uvfits / MS path or pyuvdata.UVData.  Raises core.BandwidthError when the
    band is too narrow (see core.check_bandwidth)."""
    cube = load(source, unproject=unproject, pols=pols, flag_edge_chans=flag_edge_chans)
    rx = read_metafits_rx(metafits) if metafits else None
    return analyse(cube, params, rx)


def tile_summary(rows):
    """Per pol: tiles ranked by number of flagged baselines.  A 'hub' (count >= 2 and
    strictly more than each partner) is the likely leak source: a leak of tile A's signal
    into inputs B, C, D appears as a fixed-phase line on A-B, A-C, A-D only."""
    out = {}
    for pol in sorted({r["pol"] for r in rows}):
        fl = [r for r in rows if r["pol"] == pol and r["flag"]]
        cnt, partners, same, rxs = defaultdict(int), defaultdict(list), defaultdict(int), {}
        for r in fl:
            for a, b, ra, sa in ((r["ant1"], r["ant2"], r["rx1"], r["slot1"]),
                                 (r["ant2"], r["ant1"], r["rx2"], r["slot2"])):
                cnt[a] += 1
                partners[a].append(b)
                same[a] += r["same_rx"]
                rxs[a] = (ra, sa)
        tiles = []
        for t in sorted(cnt, key=lambda t: -cnt[t]):
            hub = cnt[t] >= 2 and all(cnt[t] > cnt[p] for p in partners[t])
            tiles.append(dict(tile=t, rx=rxs[t][0], slot=rxs[t][1], n_flagged=cnt[t],
                              n_same_rx=same[t], hub=hub, partners=partners[t]))
        out[pol] = tiles
    return out
