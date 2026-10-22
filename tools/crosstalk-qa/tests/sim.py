"""Synthetic MWA-like visibilities: fringing point-source sky in the zenith drift frame,
per-antenna smooth bandpass x periodic coarse-channel (PFB-like) ripple, edge fine
channels flagged, thermal noise, and optional coherent delay lines on chosen baselines."""

from __future__ import annotations

import numpy as np

C = 299792458.0
LAT = np.deg2rad(-26.703)


def layout(n_hex_side=6, spacing=14.0, n_long=10, rng=None):
    rng = rng or np.random.default_rng(0)
    pos, names = [], []
    k = 0
    for r in range(n_hex_side):
        for c in range(n_hex_side):
            k += 1
            pos.append([c * spacing + (r % 2) * spacing / 2, r * spacing * np.sqrt(3) / 2, 0.0])
            names.append(f"Hex{k}")
    for i in range(n_long):
        rad = rng.uniform(150, 800)
        th = rng.uniform(0, 2 * np.pi)
        pos.append([rad * np.cos(th), rad * np.sin(th), rng.normal(0, 1)])
        names.append(f"Tile{i:03d}")
    return np.array(pos), names


def coarse_edge_flags(freqs, n_edge=1, cbw=1.28e6):
    cc = np.floor(freqs / cbw + 0.5).astype(int)
    m = np.zeros(len(freqs), bool)
    for c in np.unique(cc):
        idx = np.flatnonzero(cc == c)
        m[idx[:n_edge]] = m[idx[-n_edge:]] = True
    return m


def simulate(
    n_coarse=24,
    first_coarse=131,
    df=80e3,
    nt=15,
    dt=8.0,
    n_src=30,
    noise_jy=20.0,
    sky=True,
    lines=(),  # [(ant_i, ant_j, pol_index, tau_ns, amp_jy)]
    flag_coarse=(),  # coarse indices (0-based within band) fully flagged
    seed=1,
    n_hex_side=6,
    n_long=10,
    pfb_ripple=0.3,
    bright_low_jy=2000.0,
):
    rng = np.random.default_rng(seed)
    pos, names = layout(n_hex_side=n_hex_side, n_long=n_long, rng=rng)
    na = len(pos)
    nper = int(round(1.28e6 / df))
    nf = n_coarse * nper
    # fine channel centres: coarse centre +/- ... (MWAX: centred on coarse centre)
    freqs = (first_coarse * 1.28e6 - 0.64e6 + (np.arange(nf) + 0.5) * df)
    bl = [(i, j) for i in range(na) for j in range(i + 1, na)]
    a1 = np.array([b[0] for b in bl])
    a2 = np.array([b[1] for b in bl])
    bvec = pos[a2] - pos[a1]  # (Nbl, 3) ENU
    vis = np.zeros((len(bl), nt, nf, 2), np.complex128)
    if sky:
        dec = np.arcsin(rng.uniform(-1, 0.3, n_src))
        ha0 = rng.uniform(-1.2, 1.2, n_src)
        flux = 10 ** rng.uniform(0.5, 2.5, n_src)
        # one bright source low on the horizon
        dec[0], ha0[0], flux[0] = np.deg2rad(-26.7), np.deg2rad(75.0), bright_low_jy
        t = np.arange(nt) * dt
        ha = ha0[None, :] + 2 * np.pi * t[:, None] / 86164.0  # (nt, nsrc)
        e = -np.cos(dec) * np.sin(ha)
        n = np.cos(LAT) * np.sin(dec) - np.sin(LAT) * np.cos(dec) * np.cos(ha)
        u = np.sin(LAT) * np.sin(dec) + np.cos(LAT) * np.cos(dec) * np.cos(ha)
        up = (u > 0.02).all(0)
        s = np.stack([e, n, u], -1)[:, up]  # (nt, ns, 3)
        S = flux[up][None, :] * (freqs[:, None] / 150e6) ** -0.8  # (nf, ns)
        for p in range(2):
            beam = np.clip(s[..., 2], 0, 1) ** (2 + p)  # crude, pol-dependent
            for it in range(nt):
                tau = bvec @ s[it].T / C  # (Nbl, ns)
                ph = np.exp(-2j * np.pi * freqs[None, :, None] * tau[:, None, :])
                vis[:, it, :, p] = (ph * (S * beam[it])[None]).sum(-1)
        # uncalibrated instrument: smooth per-antenna bandpass x coarse-channel ripple
        x = (np.arange(nf) % nper + 0.5) / nper
        pfb = 1 - pfb_ripple * (2 * x - 1) ** 6
        for p in range(2):
            ph0 = rng.uniform(0, 2 * np.pi, na)
            amp = 1 + 0.1 * np.sin(2 * np.pi * freqs[None, :] / 25e6 + ph0[:, None])
            dly = rng.normal(0, 2e-9, na)  # residual cable delay errors (~2 ns)
            g = pfb[None] * amp * np.exp(2j * np.pi * freqs[None] * dly[:, None])
            vis[..., p] *= (g[a1] * np.conj(g[a2]))[:, None, :]
    vis += noise_jy * (rng.normal(size=vis.shape) + 1j * rng.normal(size=vis.shape)) / np.sqrt(2)
    idx = {(i, j): k for k, (i, j) in enumerate(bl)}
    for i, j, p, tau_ns, amp in lines:
        phase0 = rng.uniform(0, 2 * np.pi)
        line = amp * np.exp(1j * (2 * np.pi * freqs * tau_ns * 1e-9 + phase0))
        vis[idx[(i, j)], :, :, p] += line[None, :]
    flags = np.zeros(vis.shape, bool)
    flags[:, :, coarse_edge_flags(freqs)] = True
    for c in flag_coarse:
        flags[:, :, c * nper : (c + 1) * nper] = True
    vis[flags] = 0
    length = np.linalg.norm(bvec, axis=1)
    return dict(vis=vis.astype(np.complex64), flags=flags, freqs=freqs, a1=a1, a2=a2,
                names=names, pos=pos, length=length, idx=idx, nt=nt)
