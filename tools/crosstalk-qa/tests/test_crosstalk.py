import numpy as np
import pytest

from crosstalk import core
from crosstalk.api import Result, scan, tile_summary
from crosstalk.io import coarse_edge_mask
from sim import simulate

HUB = 7  # Hex8 in a 5x5 hex: row 1, column 2
# leak of tile HUB into its row neighbours: 14, 28, 42 m, XX only (pol index 0)
LINES = [(HUB, HUB + 1, 0, 195.0, 6.0), (HUB, HUB + 2, 0, 195.0, 6.0), (HUB, HUB + 3, 0, 228.0, 6.0)]
SIM = dict(n_hex_side=5, n_long=6, nt=10, n_src=15)


def run(d, params=None):
    params = params or core.Params()
    out = []
    for p in range(2):
        s, spec, tau = core.run_pol(d["vis"][..., p], d["flags"][..., p], d["freqs"], d["length"], params)
        out.append((s, core.decide(s, d["nt"], params), spec, tau))
    return out


def flagged_set(d, res, p):
    return {(int(d["a1"][i]), int(d["a2"][i])) for i in np.flatnonzero(res[p][1])}


@pytest.fixture(scope="module")
def with_lines():
    d = simulate(lines=LINES, **SIM)
    return d, run(d)


@pytest.fixture(scope="module")
def without_lines():
    d = simulate(lines=(), **SIM)
    return d, run(d)


def test_flags_exactly_the_injected_xx_baselines(with_lines):
    d, res = with_lines
    assert flagged_set(d, res, 0) == {(i, j) for i, j, *_ in LINES}
    assert flagged_set(d, res, 1) == set()
    s = res[0][0]
    for i, j, _, tau0, _ in LINES:
        k = d["idx"][(i, j)]
        assert abs(s.tau_peak_ns[k] - tau0) <= 1e9 / (d["freqs"].size * 80e3)  # within 1 bin
        assert s.coherence[k] > 0.7


def test_tile_summary_names_the_hub(with_lines):
    d, res = with_lines
    rows = []
    for p, pol in enumerate(["XX", "YY"]):
        for i in np.flatnonzero(res[p][1]):
            rows.append(dict(pol=pol, ant1=d["names"][d["a1"][i]], ant2=d["names"][d["a2"][i]], rx1=None,
                             slot1=None, rx2=None, slot2=None, same_rx=False, flag=True))
    tiles = tile_summary(rows)["XX"]
    assert tiles[0]["tile"] == d["names"][HUB] and tiles[0]["hub"] and tiles[0]["n_flagged"] == 3
    assert not any(t["hub"] for t in tiles[1:])


def test_foreground_and_bandpass_alone_not_flagged(without_lines):
    d, res = without_lines
    assert flagged_set(d, res, 0) == set() and flagged_set(d, res, 1) == set()


def test_raw_delay_cut_would_flag_foreground(without_lines):
    """Why the foreground is projected out: with a plain 50 ns cut on the windowed
    spectrum, the Blackman-Harris main lobe of the foreground is far above the noise just
    beyond the cut on the short baselines."""
    d, _ = without_lines
    k = d["idx"][(0, 1)]  # 14 m
    v = d["vis"][k, :, :, 0].mean(0) * core.blackman_harris(d["freqs"].size)
    p = np.abs(np.fft.fftshift(np.fft.fft(v))) ** 2
    tau = core.delay_axis_ns(d["freqs"])
    floor = np.median(p[np.abs(tau) > 3000]) / np.log(2)
    beyond = (np.abs(tau) > d["length"][k] / core.C_M_PER_NS + 50) & (np.abs(tau) < 200)
    assert p[beyond].max() / floor > 25


def test_long_baseline_foreground_not_flagged():
    """Long baselines with a bright near-horizon source and no PFB ripple (so they are
    searchable without harmonic masks)."""
    d = simulate(lines=(), pfb_ripple=0.0, bright_low_jy=5000.0, **SIM)
    params = core.Params(coarse_bw_hz=None)
    res = run(d, params)
    long_ = d["length"] > 150
    s = res[0][0]
    assert s.searchable[long_].sum() > 50
    assert not res[0][1].any() and not res[1][1].any()


def test_noise_only_statistic_is_unit_exponential():
    d = simulate(sky=False, lines=(), **SIM)
    res = run(d)
    assert not res[0][1].any() and not res[1][1].any()
    s, _, spec, _ = res[0]
    vals = spec.t[spec.noise]
    assert abs(np.mean(vals) - 1) < 0.1
    assert abs(np.mean(vals > 5) - np.exp(-5)) < 0.003


def test_sign_convention():
    rng = np.random.default_rng(3)
    freqs = 167.04e6 + (np.arange(384) + 0.5) * 80e3
    v = np.exp(2j * np.pi * freqs * 300e-9)[None, None, :] + 0.1 * (
        rng.normal(size=(1, 8, 384)) + 1j * rng.normal(size=(1, 8, 384)))
    f = np.zeros(v.shape, bool)
    L = np.array([14.0])
    p = core.Params(coarse_bw_hz=None)
    s, *_ = core.run_pol(v, f, freqs, L, p)
    assert abs(s.tau_peak_ns[0] - 300) < 20
    s, *_ = core.run_pol(v.conj(), f, freqs, L, p)
    assert abs(s.tau_peak_ns[0] + 300) < 20


def test_refuses_single_coarse_channel():
    freqs = 167.04e6 + (np.arange(16) + 0.5) * 80e3
    with pytest.raises(core.BandwidthError, match="781 ns"):
        core.check_bandwidth(freqs, core.Params())
    w = core.check_bandwidth(freqs, core.Params(force=True))
    assert w and w[0].startswith("FORCED")
    # 8 coarse channels (10.24 MHz) is accepted but warned about
    w = core.check_bandwidth(167.04e6 + (np.arange(128) + 0.5) * 80e3, core.Params())
    assert any("lose sensitivity" in x for x in w)
    assert core.check_bandwidth(167.04e6 + (np.arange(384) + 0.5) * 80e3, core.Params()) == []


def test_coarse_edge_mask():
    freqs = 167.04e6 + (np.arange(32) + 0.5) * 80e3
    m = coarse_edge_mask(freqs, 1)
    assert np.flatnonzero(m).tolist() == [0, 15, 16, 31]


def test_scan_uvfits_roundtrip(tmp_path):
    pyuvdata = pytest.importorskip("pyuvdata")
    from astropy.time import Time
    from pyuvdata import Telescope, UVData
    from pyuvdata.utils import ECEF_from_ENU

    d = simulate(lines=LINES, **SIM)
    known = Telescope.from_known_telescopes("mwa")
    loc = known.location
    ecef = ECEF_from_ENU(d["pos"], center_loc=loc)
    tel = Telescope.new(
        name="MWA", location=loc, instrument="MWA", antenna_positions=ecef - np.array(known._location.xyz()),
        antenna_names=list(d["names"]), antenna_numbers=np.arange(len(d["names"])), x_orientation="east",
        update_from_known=False,
    )
    times = Time(2460900.0, format="jd").jd + np.arange(d["nt"]) * 8.0 / 86400
    nbl = len(d["a1"])
    data = d["vis"].transpose(1, 0, 2, 3).reshape(-1, d["freqs"].size, 2)  # time-major blts
    flags = d["flags"].transpose(1, 0, 2, 3).reshape(data.shape)
    uvd = UVData.new(
        freq_array=d["freqs"], polarization_array=[-5, -6], telescope=tel, times=times,
        antpairs=list(zip(d["a1"].tolist(), d["a2"].tolist())), do_blt_outer=True,
        data_array=data, flag_array=flags, channel_width=80e3, integration_time=8.0,
        time_axis_faster_than_bls=False,
    )
    path = str(tmp_path / "sim.uvfits")
    uvd.write_uvfits(path, force_phase=True)  # uvfits must be phased; scan() unprojects
    res = scan(path)
    assert isinstance(res, Result)
    got = {(r["ant1"], r["ant2"]) for r in res.flagged}
    want = {(d["names"][i], d["names"][j]) for i, j, *_ in LINES}
    assert {r["pol"] for r in res.flagged} == {"XX"}
    assert got == want
    hub = [t for t in res.tiles["XX"] if t["hub"]]
    assert [t["tile"] for t in hub] == [d["names"][HUB]]
    # also accepts a UVData object directly
    assert {(r["ant1"], r["ant2"]) for r in scan(uvd).flagged} == want
