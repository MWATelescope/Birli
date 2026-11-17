"""Load uvfits / MS (Birli output) into per-pol (Nbl, Nt, Nf) cubes via pyuvdata."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)

POL_NAMES = {-5: "XX", -6: "YY", -7: "XY", -8: "YX", 1: "I", 2: "Q", 3: "U", 4: "V"}


@dataclass
class Cube:
    vis: np.ndarray  # (Nbl, Nt, Nf, Npol) complex64
    flags: np.ndarray  # (Nbl, Nt, Nf, Npol) bool
    freqs_hz: np.ndarray  # (Nf,)
    pols: list[str]
    ant1: np.ndarray  # (Nbl,) antenna numbers
    ant2: np.ndarray
    name1: list[str]
    name2: list[str]
    length_m: np.ndarray  # (Nbl,)
    times_jd: np.ndarray
    frame: str


def coarse_edge_mask(freqs_hz: np.ndarray, n_edge: int, coarse_bw_hz: float = 1.28e6) -> np.ndarray:
    """True for the n_edge fine channels at each end of every MWA coarse channel
    (coarse centres at integer multiples of coarse_bw)."""
    m = np.zeros(len(freqs_hz), bool)
    if n_edge <= 0:
        return m
    cc = np.floor(freqs_hz / coarse_bw_hz + 0.5).astype(int)
    for c in np.unique(cc):
        idx = np.flatnonzero(cc == c)
        m[idx[:n_edge]] = True
        m[idx[-n_edge:]] = True
    return m


def load(source, unproject: bool = True, pols=("XX", "YY"), flag_edge_chans: int = 0) -> Cube:
    """source: path to uvfits / MS (anything pyuvdata reads) or a pyuvdata.UVData
    (copied, not modified)."""
    from pyuvdata import UVData

    uvd = source.copy() if isinstance(source, UVData) else UVData.from_file(source)
    frame = "as-is"
    if unproject:
        cats = [c.get("cat_type") for c in (uvd.phase_center_catalog or {}).values()]
        if any(c != "unprojected" for c in cats):
            uvd.unproject_phase()
            frame = "unprojected (zenith drift)"
        else:
            frame = "unprojected (already)"
    uvd.select(ant_str="cross")

    pol_idx = []
    for p in pols:
        hits = [i for i, n in enumerate(uvd.polarization_array) if POL_NAMES.get(int(n)) == p]
        if hits:
            pol_idx.append(hits[0])
        else:
            log.warning("pol %s not in file", p)
    ub, binv = np.unique(uvd.baseline_array, return_inverse=True)
    ut, tinv = np.unique(uvd.time_array, return_inverse=True)
    nbl, nt, nf = len(ub), len(ut), uvd.Nfreqs
    vis = np.zeros((nbl, nt, nf, len(pol_idx)), np.complex64)
    flags = np.ones((nbl, nt, nf, len(pol_idx)), bool)
    vis[binv, tinv] = uvd.data_array[:, :, pol_idx]
    flags[binv, tinv] = uvd.flag_array[:, :, pol_idx] | (uvd.nsample_array[:, :, pol_idx] <= 0)
    freqs = np.asarray(uvd.freq_array).ravel()
    order = np.argsort(freqs)
    if not np.all(order == np.arange(nf)):
        vis, flags, freqs = vis[:, :, order], flags[:, :, order], freqs[order]
    if flag_edge_chans:
        flags |= coarse_edge_mask(freqs, flag_edge_chans)[None, None, :, None]

    a1, a2 = uvd.baseline_to_antnums(ub)
    num2name = dict(zip(uvd.telescope.antenna_numbers, uvd.telescope.antenna_names))
    num2pos = dict(zip(uvd.telescope.antenna_numbers, uvd.telescope.antenna_positions))
    length = np.array([np.linalg.norm(num2pos[j] - num2pos[i]) for i, j in zip(a1, a2)])
    return Cube(
        vis=vis,
        flags=flags,
        freqs_hz=freqs,
        pols=[POL_NAMES[int(uvd.polarization_array[i])] for i in pol_idx],
        ant1=np.asarray(a1),
        ant2=np.asarray(a2),
        name1=[num2name[i] for i in a1],
        name2=[num2name[j] for j in a2],
        length_m=length,
        times_jd=ut,
        frame=frame,
    )


def read_metafits_rx(path: str) -> dict[str, tuple[int, int, int]]:
    """TileName -> (metafits Antenna index, Rx, Slot)."""
    from astropy.io import fits

    with fits.open(path) as h:
        t = h["TILEDATA"].data if "TILEDATA" in h else h[1].data
        out = {}
        for r in t:
            if str(r["Pol"]).strip() == "X":
                out[str(r["TileName"]).strip()] = (int(r["Antenna"]), int(r["Rx"]), int(r["Slot"]))
    return out
