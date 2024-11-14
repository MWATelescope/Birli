from pyuvdata import UVData

default_read_args = dict(
    use_aoflagger_flags=False,
    remove_dig_gains=False,
    remove_coarse_band=False,
    correct_cable_len=False,
    correct_van_vleck=False,
    cheby_approx=False,
    flag_small_auto_ants=False,
    phase_to_pointing_center=False,
    propagate_coarse_flags=False,
    flag_init=False,
    remove_flagged_ants=False,
    run_check=False,
    fix_autos=False,
)


def dump_uvd(uvd: UVData, dump_csv, row_limit=250):
    row_limit = min(uvd.Nblts, row_limit)
    header = ["timestep", "baseline", "u", "v", "w", "pol", "type"]
    blts = range(min(uvd.Nblts, row_limit))
    times = uvd.time_array[:row_limit]
    # get the index of each antenna in the ant number array
    ant1s = [
        list(uvd.antenna_numbers).index(ant) + 1 for ant in uvd.ant_1_array[:row_limit]
    ]
    ant2s = [
        list(uvd.antenna_numbers).index(ant) + 1 for ant in uvd.ant_2_array[:row_limit]
    ]
    bls = [ant1 * 256 + ant2 for ant1, ant2 in zip(ant1s, ant2s)]
    uvws = -1 * uvd.uvw_array[:row_limit, :] / 299_792_458
    with open(dump_csv, "w") as dump:
        dump.write(",".join(header + [*map(str, range(uvd.Nfreqs))]) + "\n")
        for blt, time, bl, uvw in zip(blts, times, bls, uvws):
            row = [time, bl, *uvw]
            for pol_idx, pol in enumerate(["xx", "yy", "xy", "yx"]):
                vis = uvd.data_array[blt, :, pol_idx].tolist()
                dump.write(",".join(map(str, row[:] + [pol, "vis"] + vis)) + "\n")
    print(f"Wrote {dump.name}")
