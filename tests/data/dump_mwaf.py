#!/usr/bin/env python

from itertools import combinations_with_replacement, combinations
from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs

from astropy.io import fits
from pprint import pformat
import sys
from argparse import ArgumentParser
from tabulate import tabulate


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument("file", type=str)
    parser.add_argument("--timestep-limit", type=int, default=0)
    parser.add_argument("--baseline-limit", type=int, default=0)
    return parser.parse_args(argv)


def chunk(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)


def print_heading(heading):
    print("")
    print(f"{'-' * len(heading)}")
    print(heading)
    print(f"{'-' * len(heading)}")


def main(argv):
    args = parse_args(argv)
    print(f"-> args {vars(args)}")
    hdus = fits.open(args.file)
    print(f"-> hdus.info():")
    hdus.info()

    print_heading("HEADER")

    print(repr(hdus["PRIMARY"].header))

    num_scans = hdus["PRIMARY"].header["NSCANS"]
    print(f" -> num scans: hdus['PRIMARY'].header['NSCANS']={num_scans}")
    num_antenna = hdus["PRIMARY"].header["NANTENNA"]
    print(f" -> num scans: hdus['PRIMARY'].header['NANTENNA']={num_antenna}")
    num_chans = hdus["FLAGS"].header["NAXIS1"] * 8
    print(f" -> num chans: hdus['FLAGS'].header['NAXIS1']*8={num_chans}")
    num_rows = hdus["FLAGS"].header["NAXIS2"]
    assert num_rows % num_scans == 0
    num_baselines = num_rows // num_scans
    print(f" -> num baselines: {num_baselines}")

    ant_pairs = [*combinations_with_replacement(range(num_antenna), 2)]
    if len(ant_pairs) != num_baselines:
        ant_pairs_noautos = [*combinations(range(num_antenna), 2)]
        if len(ant_pairs_noautos) != num_baselines:
            raise ValueError(
                f"num_baselines={num_baselines}"
                f" != len(ant_pairs_noautos)={len(ant_pairs_noautos)}"
                f" or len(ant_pairs_autos)={len(ant_pairs)}"
            )
        ant_pairs = ant_pairs_noautos

    print_heading("FLAG DATA")

    print(repr(hdus["FLAGS"].header))

    flag_data = hdus["FLAGS"].data
    print(f"flags shape {flag_data.shape}")
    assert flag_data.shape[0] == num_rows

    timestep_limit = num_rows // num_baselines
    if num_rows != num_scans * num_baselines:
        print(f"num_scans ({num_scans}) * num_baselines({num_baselines}) != num_rows({num_rows})")
    print(f"actual num_scans: {timestep_limit}")
    if args.timestep_limit:
        timestep_limit = min(args.timestep_limit, timestep_limit)
    baseline_limit = num_baselines
    if args.baseline_limit:
        baseline_limit = min(args.baseline_limit, baseline_limit)
    print(f"-> limits: baseline={baseline_limit}, timestep={timestep_limit}")
    for baseline_idx, (ant1, ant2) in enumerate(ant_pairs):
        print(f"-> bl {baseline_idx:04d} ({ant1:03d},{ant2:03d}):")
        for timestep_idx in range(timestep_limit):
            flags = flag_data[timestep_idx * num_baselines + baseline_idx][0]
            flag_display = "".join([f"{'#' if flag else '.'}" for flag in flags])
            print(f" --> ts {timestep_idx:04d}: {flag_display}")

    if "CH_OCC" in hdus:
        print_heading("CH_OCC")
        print(repr(hdus["CH_OCC"].header))
        ch_occ_data = hdus["CH_OCC"].data
        print(f"ch_occ shape {ch_occ_data.shape}")
        assert ch_occ_data.shape[0] == num_chans
        print(tabulate(ch_occ_data, headers=["ch", "count", "occupancy"]))
    if "BL_OCC" in hdus:
        print_heading("BL_OCC")
        print(repr(hdus["BL_OCC"].header))
        bl_occ_data = hdus["BL_OCC"].data[0:baseline_limit]
        print(f"bl_occ shape {bl_occ_data.shape}")
        assert bl_occ_data.shape[0] == num_baselines
        print(tabulate(bl_occ_data, headers=["bl", "count", "occupancy"]))
    if "TILES" in hdus:
        print_heading("TILES")
        print(repr(hdus["TILES"].header))
        tiles_data = hdus["TILES"].data
        print(f"tiles shape {tiles_data.shape}")
        assert tiles_data.shape[0] == num_antenna
        print(tabulate(tiles_data, headers=["ant", "name"]))


if __name__ == "__main__":
    main(sys.argv[1:])
    # main(["tests/data/1196175296_mwa_ord/FlagfileCotter01.mwaf"])
