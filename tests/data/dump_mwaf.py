#!/usr/bin/env python

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
    parser.add_argument(
        "file",
        type=str
    )
    parser.add_argument(
        "--timestep-limit",
        type=int,
        default=0
    )
    parser.add_argument(
        "--baseline-limit",
        type=int,
        default=0
    )
    return parser.parse_args(argv)

def chunk(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)

def main(argv):
    args = parse_args(argv)
    print(f"-> args {vars(args)}")
    hdus = fits.open(args.file)
    print(f"-> hdus.info():")
    hdus.info()

    ####
    # Primary HDU
    #### 

    primary_hdu = hdus[0]
    print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")
    # hdus[0].info()
    num_scans = primary_hdu.header['NSCANS']
    num_antenna = primary_hdu.header['NANTENNA'] 
    num_baselines = (num_antenna * (num_antenna + 1) // 2)
    print(f" -> num baselines: {num_baselines}")

    # ####
    # # Flag Data
    # ####

    # flag_header = hdus[1].data.dtype.names
    flag_data = hdus[1].data
    print(f"flags shape {flag_data.shape}")
    rows = flag_data.shape[0]
    
    assert rows == num_scans * num_baselines
    timestep_limit = num_scans
    if args.timestep_limit:
        timestep_limit = min(args.timestep_limit, timestep_limit)
    baseline_limit = num_baselines
    if args.baseline_limit:
        baseline_limit = min(args.baseline_limit, baseline_limit)
    print(f"-> limits: baseline={baseline_limit}, timestep={timestep_limit}")
    for baseline_idx in range(baseline_limit):
        print(f"-> bl {baseline_idx:04d}:")
        for timestep_idx in range(timestep_limit):
            flags = flag_data[timestep_idx * num_baselines + baseline_idx][0]
            flag_display = "".join([f"{'#' if flag else '.'}" for flag in flags])
            print(f" --> ts {timestep_idx:04d}: {flag_display}")    

if __name__ == '__main__':
    main(sys.argv[1:])
    # main(["tests/data/1196175296_mwa_ord/FlagfileCotter01.mwaf"])

