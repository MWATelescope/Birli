"""Per-baseline crosstalk (fixed-phase delay line) detector -- array-level core.

Pure numpy (+ scipy's DPSS), no I/O, written to port 1:1 to Rust (NOTES.md,
"Crosstalk delay QA").

Prior art: first reported, without attribution to a tile or receiver, by Yoshiura et al.
(2026, PASJ 78, 1931, sec 5.3, Figs 12-14; doi:10.1093/pasj/psag110); a related but more
widespread, dual-polarisation non-fringing common mode in 2016 MWA Phase II data was
reported by Kolopanis et al. (2023, MNRAS 521, 5120; arXiv:2210.10885).  The
HexE11 / receiver-23 attribution is from our own investigation.

Model
-----
A receiver crosstalk leak of tile A's signal into input B adds eps*<x_A x_A*> to V_AB:
a fixed-phase, non-fringing term V_AB(t, f) += a * exp(2*pi*i*f*tau0) -- a narrow line
at tau0 in delay that survives coherent time averaging.  Sky foregrounds occupy
|tau| <= |b|/c in the zenith (drift) frame, broadened by spectral structure.  "Flag the
baseline if the FFT peak is off centre" is made rigorous as a test for one extra
delay line beyond horizon+buffer on top of a foreground band-limited to horizon+buffer.

Algorithm (per baseline, per pol)
---------------------------------
1. Coherent mean V(f) = weighted time mean over unflagged samples.  "Good" channels are
   unflagged in >= 50 % of the valid timesteps; others are excluded from the fit (not
   zero-filled into it) and zeroed in the transform.  For the per-timestep spectra,
   flagged samples of good channels are imputed with V(f).
2. Foreground subspace: DPSS (Slepian) vectors for Nf channels band-limited to
   |tau| <= tau_cut = |b|/c + buffer, K = ceil(2*B*tau_cut) + extra modes, restricted to the
   good channels and orthonormalised (QR) -> Q.  Residual r = (I - Q Q^H) V.  This removes
   the foreground *and* its window/flag sidelobes, unlike a delay cut on the raw spectrum
   (whose Blackman-Harris main lobe alone reaches ~4 bins = 130 ns beyond the horizon at
   30.72 MHz, at ~25 dB above the noise in the synthetic tests).
3. R(tau) = FFT(W r), W = Blackman-Harris taper (suppresses sidelobes of residual
   structure, e.g. PFB-ripple copies of the foreground), zero-padded x pad.
4. Separability sep(tau) = 1 - sum_k |FFT(W q_k)(tau)|^2 / sum W^2: the fraction of the
   noise in R(tau) that survives the projection (and approximately of a line's power).
   sep -> 0 inside the cut; this is the bandwidth limit made explicit (one 1.28 MHz
   coarse channel: sep ~ 0 out to ~800 ns).  T(tau) = |R|^2 / sep / floor is then
   exponentially distributed with unit mean under noise.
5. floor = median(|R|^2/sep over |tau| >= noise_min, outside masks) / ln 2.
6. Search region: |tau| > tau_cut and not within tau_cut of a coarse-channel harmonic
   k/coarse_bw (781.25 ns: periodic edge flags and the PFB passband copy the foreground
   there); sep >= sep_min (0.1) only for numerical safety.  Peak = highest *local
   maximum* of T in the region, with masked bins counting as +inf, so the bin next to
   the cut is never a peak and a residual skirt decaying from the cut is rejected.
7. Coherence at the peak from per-timestep residual spectra: |<R_t>|^2/<|R_t|^2>
   (noise ~ 1/Nt, fixed-phase line -> 1).
8. Robust z of log10(T) against the k_neighbours baselines nearest in length (same pol).
9. flag = T >= t_min and z >= z_min and coherence >= coherence_min
          and onesided >= 6 dB [and T / T_otherpol(same delay) >= asym_min].
   onesided = T(tau_peak) / max(T(-tau_peak), 1): an additive line adds power on one side
   of delay only, whereas a foreground skirt leaking past the cut is roughly two-sided.

Thresholds: under noise T ~ Exp(1), P(T > t) = e^-t.  A full MWA scan is ~8k baselines
x 2 pols x ~100 independent search bins ~ 2e6 trials, so t_min = 25 gives ~3e-5 expected
false flags from Gaussian noise alone.  z_min = 5 and coherence_min = 0.3 guard against
non-Gaussian systematics shared by similar baselines and incoherent (transient) RFI.

Delay sign: R(tau) = sum_f W r(f) exp(-2 pi i f tau), so a phase increasing with frequency
at 2*pi*tau0 peaks at +tau0.  Swapping the baseline order flips the sign.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C_M_PER_NS = 0.299792458
LN2 = np.log(2.0)
MWA_COARSE_BW_HZ = 1.28e6


class BandwidthError(ValueError):
    """Bandwidth too narrow to separate a near-horizon line from the foreground."""


@dataclass
class Params:
    buffer_ns: float = 50.0  # foreground support = |b|/c + buffer
    fg_extra_modes: int = 4  # DPSS modes beyond ceil(2*B*tau_cut)
    sep_min: float = 0.1  # numerical floor on separability for searched delays
    noise_min_ns: float = 3000.0  # floor region |tau| >= this
    coarse_bw_hz: float | None = MWA_COARSE_BW_HZ  # None -> no harmonic masking
    pad: int = 2  # FFT zero-padding factor
    t_min: float = 25.0  # T = power / (sep * floor) at the peak
    z_min: float = 5.0  # robust z vs similar-length baselines
    coherence_min: float = 0.3  # skipped if Nt < 4
    asym_min: float | None = None  # optional T / T_otherpol at the same delay
    k_neighbours: int = 50
    mad_floor_dex: float = 0.05
    min_noise_bins: int = 16
    max_resolution_ns: float = 100.0  # refuse B < 10 MHz
    min_chan_time_frac: float = 0.5  # channel used if unflagged in >= this fraction of timesteps
    onesided_min_db: float | None = 6.0  # T(tau_peak) over T(-tau_peak)
    force: bool = False


@dataclass
class Spectra:
    t: np.ndarray  # (Nbl, Ntau) normalised coherent statistic T(tau)
    sep: np.ndarray  # (Nbl, Ntau)
    coherence: np.ndarray  # (Nbl, Ntau)
    t_inc: np.ndarray  # (Nbl, Ntau) per-timestep mean power / its floor
    search: np.ndarray  # (Nbl, Ntau) bool
    noise: np.ndarray  # (Nbl, Ntau) bool
    lowsep: np.ndarray  # (Nbl, Ntau) bool: outside the masks but sep < sep_min
    raw: np.ndarray  # (Nbl, Ntau) windowed coherent power WITHOUT projection / its floor
    floor: np.ndarray  # (Nbl,)
    n_times: np.ndarray  # (Nbl,)
    good_frac: np.ndarray  # (Nbl,) fraction of channels used


@dataclass
class Scores:
    t: np.ndarray
    tau_peak_ns: np.ndarray
    sep: np.ndarray
    coherence: np.ndarray
    t_inc: np.ndarray
    floor: np.ndarray
    n_search: np.ndarray
    good_frac: np.ndarray
    searchable: np.ndarray
    tau_raw_ns: np.ndarray  # raw-spectrum maximum up to 4 bins inward of the peak (same side)
    raw_db: np.ndarray  # its power over the raw floor [dB]
    raw_onesided_db: np.ndarray  # raw power at tau_raw over power at -tau_raw [dB]
    onesided_db: np.ndarray  # T(tau_peak) / max(T(-tau_peak), 1) [dB], on the residual
    z: np.ndarray = field(default=None)


def check_bandwidth(freqs_hz: np.ndarray, params: Params) -> list[str]:
    """Delay resolution is 1/B.  Refuse when 1/B > max_resolution_ns (100 ns, B < 10 MHz):
    a single 1.28 MHz coarse channel has 781 ns resolution and cannot separate a +195 ns
    line from a foreground ending at 50-150 ns.  Warn when a line ~150 ns beyond the
    horizon (the HexE11-HexE12 case) would have < 3.5 bins clearance."""
    df = np.median(np.diff(freqs_hz))
    bw = df * len(freqs_hz)
    res_ns = 1e9 / bw
    warnings = []
    if res_ns > params.max_resolution_ns:
        msg = (
            f"bandwidth {bw / 1e6:.3f} MHz gives delay resolution {res_ns:.0f} ns > "
            f"{params.max_resolution_ns:.0f} ns: lines within ~{3.5 * res_ns:.0f} ns of the "
            f"horizon are inseparable from the foreground (need >= {1e3 / params.max_resolution_ns:.2f} "
            f"MHz = {int(np.ceil(1e9 / params.max_resolution_ns / MWA_COARSE_BW_HZ))} contiguous "
            "coarse channels; the full 24 is recommended)"
        )
        if not params.force:
            raise BandwidthError(msg)
        warnings.append("FORCED: " + msg)
    elif 3.5 * res_ns > 150.0:
        warnings.append(
            f"delay resolution {res_ns:.0f} ns: lines closer than ~{3.5 * res_ns:.0f} ns to the "
            "horizon lose sensitivity (HexE11-like +195 ns lines on 14 m baselines want >= 23 MHz)"
        )
    gaps = np.abs(np.diff(freqs_hz) - df) > 1e-3 * df
    if gaps.any():
        warnings.append("frequency axis is not contiguous; the delay transform assumes it is")
    return warnings


def blackman_harris(n: int) -> np.ndarray:
    """4-term Blackman-Harris (== scipy.signal.windows.blackmanharris(n, sym=True))."""
    a = (0.35875, 0.48829, 0.14128, 0.01168)
    if n == 1:
        return np.ones(1)
    x = 2 * np.pi * np.arange(n) / (n - 1)
    return a[0] - a[1] * np.cos(x) + a[2] * np.cos(2 * x) - a[3] * np.cos(3 * x)


def delay_axis_ns(freqs_hz: np.ndarray, pad: int = 1) -> np.ndarray:
    df = np.median(np.diff(freqs_hz))
    return np.fft.fftshift(np.fft.fftfreq(len(freqs_hz) * pad, d=df)) * 1e9


_DPSS_CACHE: dict = {}


def n_modes(nf: int, df_hz: float, tau_cut_ns: float, extra_modes: int) -> int:
    return int(np.ceil(2 * nf * df_hz * tau_cut_ns * 1e-9)) + extra_modes


def dpss_basis(nf: int, df_hz: float, tau_cut_ns: float, extra_modes: int) -> np.ndarray:
    """(nf, K) DPSS basis band-limited to |tau| <= tau_cut: NW = nf*df*tau_cut,
    K = ceil(2 NW) + extra.  Cached on tau_cut rounded up to 5 ns."""
    from scipy.signal.windows import dpss

    tq = 5.0 * np.ceil(tau_cut_ns / 5.0)
    key = (nf, round(df_hz, 3), tq, extra_modes)
    if key not in _DPSS_CACHE:
        nw = nf * df_hz * tq * 1e-9
        k = min(n_modes(nf, df_hz, tq, extra_modes), nf)
        _DPSS_CACHE[key] = np.atleast_2d(dpss(nf, max(nw, 1e-3), Kmax=k)).T
    return _DPSS_CACHE[key]


def excluded_mask(tau_ns, tau_cut_ns, params: Params):
    """(Nbl, Ntau) bool: |tau| <= cut, or within cut of a coarse-channel harmonic."""
    edge = tau_cut_ns[:, None]
    t = tau_ns[None, :]
    excl = np.abs(t) <= edge
    if params.coarse_bw_hz:
        period = 1e9 / params.coarse_bw_hz
        kmax = int(np.ceil(np.abs(tau_ns).max() / period)) + 1
        for k in range(1, kmax + 1):
            excl |= np.abs(t - k * period) <= edge
            excl |= np.abs(t + k * period) <= edge
    return excl


def delay_spectra(vis, flags, freqs_hz, horizon_ns, params: Params) -> Spectra:
    """vis, flags: (Nbl, Nt, Nf) for one pol.  horizon_ns: (Nbl,)."""
    nbl, nt, nf = vis.shape
    df = float(np.median(np.diff(freqs_hz)))
    nfft = nf * params.pad
    tau = delay_axis_ns(freqs_hz, params.pad)
    w = blackman_harris(nf)
    sw2 = np.sum(w**2)
    tau_cut = horizon_ns + params.buffer_ns
    excl = excluded_mask(tau, tau_cut, params)

    good = ~flags
    tvalid = good.any(axis=2)
    ntv = tvalid.sum(axis=1)
    nt_good = good.sum(axis=1)
    # channel usable if unflagged in >= min_chan_time_frac of the valid timesteps
    chan_ok = (nt_good >= params.min_chan_time_frac * ntv[:, None]) & (nt_good > 0)

    def fft(x):
        return np.fft.fftshift(np.fft.fft(x * w, n=nfft, axis=-1), axes=-1)

    T = np.zeros((nbl, nfft))
    sep = np.zeros((nbl, nfft))
    coh = np.zeros((nbl, nfft))
    tinc = np.zeros((nbl, nfft))
    floor = np.full(nbl, np.nan)
    noise = np.zeros((nbl, nfft), bool)
    search = np.zeros((nbl, nfft), bool)
    lowsep = np.zeros((nbl, nfft), bool)
    raw = np.zeros((nbl, nfft))
    for i in range(nbl):
        g = chan_ok[i]
        if g.sum() < nf // 4 or (~excl[i]).sum() == 0:
            continue
        if n_modes(nf, df, tau_cut[i], params.fg_extra_modes) > g.sum() // 2:
            continue  # foreground model would use > half the degrees of freedom
        basis = dpss_basis(nf, df, tau_cut[i], params.fg_extra_modes)[g]
        q, _ = np.linalg.qr(basis)  # (ng, K) orthonormal on the good channels
        vt = vis[i][tvalid[i]][:, g]  # (ntv, ng)
        gt = good[i][tvalid[i]][:, g]
        vbar = np.where(gt, vt, 0).sum(0) / gt.sum(0)  # weighted mean over unflagged samples
        vt = np.where(gt, vt, vbar[None, :])  # impute flagged samples (residual ~0 there)
        rb = np.zeros(nf, complex)
        rb[g] = vbar - q @ (q.T @ vbar)
        rt = np.zeros((len(vt), nf), complex)
        rt[:, g] = vt - (vt @ q) @ q.T
        qf = np.zeros((q.shape[1], nf))
        qf[:, g] = q.T
        s = np.clip(1.0 - np.sum(np.abs(fft(qf)) ** 2, axis=0) / sw2, 0.0, 1.0)
        ok = s >= params.sep_min
        nz = ~excl[i] & ok & (np.abs(tau) >= params.noise_min_ns)
        if nz.sum() < params.min_noise_bins:
            continue
        pc = np.abs(fft(rb)) ** 2
        vb0 = np.zeros(nf, complex)
        vb0[g] = vbar
        praw = np.abs(fft(vb0)) ** 2
        raw[i] = praw / (np.median(praw[nz]) / LN2)
        Rt = fft(rt)
        pinc = np.mean(np.abs(Rt) ** 2, axis=0)
        mc = np.abs(Rt.mean(axis=0)) ** 2
        sd = np.where(s > 0, s, 1.0)
        fl = np.median((pc / sd)[nz]) / LN2
        fli = np.median((pinc / sd)[nz]) / LN2
        T[i] = np.where(s > 0, pc / sd, 0.0) / fl
        tinc[i] = np.where(s > 0, pinc / sd, 0.0) / fli
        with np.errstate(divide="ignore", invalid="ignore"):
            coh[i] = np.where(pinc > 0, mc / pinc, 0.0)
        sep[i] = s
        floor[i] = fl
        noise[i] = nz
        search[i] = ~excl[i] & ok
        lowsep[i] = ~excl[i] & ~ok
    return Spectra(T, sep, coh, tinc, search, noise, lowsep, raw, floor, ntv, chan_ok.mean(axis=1))


def raw_feature(spec: Spectra, tau_ns, horizon_ns, tau_peak_ns, win_ns):
    """Diagnostics on the unprojected spectrum: the projection only sees the part of a
    feature beyond the cut, so tau_peak can sit on the skirt of a feature centred just
    inside it.  Report the raw maximum within 4 resolution bins inward of the peak (but not
    inside the horizon, same side), its level, and its one-sidedness vs the mirror delay."""
    n = len(tau_peak_ns)
    tr, db, one = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
    for i in np.flatnonzero(np.isfinite(tau_peak_ns)):
        sgn = np.sign(tau_peak_ns[i])
        lo = max(horizon_ns[i], abs(tau_peak_ns[i]) - win_ns)
        m = (np.sign(tau_ns) == sgn) & (np.abs(tau_ns) >= lo) & (np.abs(tau_ns) <= abs(tau_peak_ns[i]))
        if not m.any():
            continue
        j = np.flatnonzero(m)[np.argmax(spec.raw[i, m])]
        jm = np.argmin(np.abs(tau_ns + tau_ns[j]))
        tr[i] = tau_ns[j]
        with np.errstate(divide="ignore"):
            db[i] = 10 * np.log10(spec.raw[i, j])
            one[i] = 10 * np.log10(spec.raw[i, j] / max(spec.raw[i, jm], 1e-30))
    return tr, db, one


def score(spec: Spectra, tau_ns: np.ndarray, params: Params, horizon_ns=None) -> Scores:
    """Highest local maximum of T in the search region.  Excluded bins count as +inf so a
    bin adjacent to the cut (or to a harmonic mask) is never a local maximum."""
    T = spec.t
    nbl = T.shape[0]
    tin = np.where(spec.search | spec.lowsep, T, np.inf)
    left = np.concatenate([np.full((nbl, 1), np.inf), tin[:, :-1]], 1)
    right = np.concatenate([tin[:, 1:], np.full((nbl, 1), np.inf)], 1)
    locmax = (T > left) & (T >= right) & spec.search
    has = locmax.any(1)
    ipk = np.argmax(np.where(locmax, T, -np.inf), 1)
    rows = np.arange(nbl)
    searchable = np.isfinite(spec.floor) & spec.search.any(1)
    pk = searchable & has
    nan = np.nan
    tau_pk = np.where(pk, tau_ns[ipk], nan)
    if horizon_ns is None:
        horizon_ns = np.zeros(nbl)
    win_ns = 4 * (tau_ns[1] - tau_ns[0]) * params.pad  # 4 resolution bins
    tr, db, one = raw_feature(spec, tau_ns, horizon_ns, tau_pk, win_ns)
    imir = np.argmin(np.abs(tau_ns[None, :] + tau_ns[ipk][:, None]), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        os_db = np.where(pk, 10 * np.log10(T[rows, ipk] / np.maximum(T[rows, imir], 1.0)), nan)
    return Scores(
        t=np.where(searchable, np.where(has, T[rows, ipk], 0.0), nan),
        tau_peak_ns=tau_pk,
        sep=np.where(pk, spec.sep[rows, ipk], nan),
        coherence=np.where(pk, spec.coherence[rows, ipk], nan),
        t_inc=np.where(pk, spec.t_inc[rows, ipk], nan),
        floor=spec.floor,
        n_search=spec.search.sum(1),
        good_frac=spec.good_frac,
        searchable=searchable,
        tau_raw_ns=tr,
        raw_db=db,
        raw_onesided_db=one,
        onesided_db=os_db,
    )


def robust_z_by_length(length_m: np.ndarray, x: np.ndarray, k: int, mad_floor: float) -> np.ndarray:
    """Robust z of x against the k baselines nearest in length (excluding self).
    NaN in x are ignored as neighbours and give NaN z."""
    z = np.full(len(x), np.nan)
    idx = np.flatnonzero(np.isfinite(x))
    if len(idx) < 3:
        return z
    order = idx[np.argsort(length_m[idx], kind="stable")]
    xs = x[order]
    n = len(order)
    kk = min(k, n - 1)
    for j in range(n):
        lo = max(0, min(j - kk // 2, n - 1 - kk))
        nb = np.r_[xs[lo:j], xs[j + 1 : lo + kk + 1]]
        med = np.median(nb)
        mad = max(1.4826 * np.median(np.abs(nb - med)), mad_floor)
        z[order[j]] = (xs[j] - med) / mad
    return z


def decide(s: Scores, n_times: int, params: Params, t_other: np.ndarray | None = None) -> np.ndarray:
    """t_other: the other pol's T at this pol's peak delay (for the asymmetry test)."""
    with np.errstate(invalid="ignore"):
        f = s.searchable & (s.t >= params.t_min) & (s.z >= params.z_min)
        if n_times >= 4 and params.coherence_min:
            f &= s.coherence >= params.coherence_min
        if params.onesided_min_db is not None:
            f &= s.onesided_db >= params.onesided_min_db
        if params.asym_min and t_other is not None:
            f &= s.t / np.maximum(np.nan_to_num(t_other, nan=1.0), 1.0) >= params.asym_min
    return f


def run_pol(vis, flags, freqs_hz, length_m, params: Params, horizon_ns=None):
    """vis/flags (Nbl, Nt, Nf) -> (Scores, Spectra, tau_ns).  horizon defaults to |b|/c."""
    if horizon_ns is None:
        horizon_ns = length_m / C_M_PER_NS
    tau = delay_axis_ns(freqs_hz, params.pad)
    spec = delay_spectra(vis, flags, freqs_hz, horizon_ns, params)
    s = score(spec, tau, params, horizon_ns)
    s.z = robust_z_by_length(length_m, np.log10(np.maximum(s.t, 1e-3)), params.k_neighbours,
                             params.mad_floor_dex)
    return s, spec, tau


def value_at(spec: Spectra, tau_ns: np.ndarray, tau_peak_ns: np.ndarray) -> np.ndarray:
    """spec.t at the given per-baseline delays (NaN delay -> NaN)."""
    out = np.full(len(tau_peak_ns), np.nan)
    ok = np.flatnonzero(np.isfinite(tau_peak_ns))
    idx = np.clip(np.searchsorted(tau_ns, tau_peak_ns[ok] - 1e-6), 0, len(tau_ns) - 1)
    out[ok] = spec.t[ok, idx]
    return out
