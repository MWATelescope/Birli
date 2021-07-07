#!/usr/bin/env python

from copy import copy
from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs
# from tests.data.dump_mwaf import chunk

from astropy.io import fits
from pprint import pformat
import numpy as np
import pandas as pd
import math
import re
import itertools
import sys
from sys import stdout
from argparse import ArgumentParser
import pyuvdata
from pyuvdata import UVData
from tabulate import tabulate

from common import eprint


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "uvfile",
        type=str
    )
    parser.add_argument(
        "--timestep-limit",
        type=int,
        default=None
    )
    parser.add_argument(
        "--baseline-limit",
        type=int,
        default=None
    )
    parser.add_argument(
        "--chan-limit",
        type=int,
        default=None
    )
    parser.add_argument(
        "--antenna-limit",
        type=int,
        default=None
    )
    parser.add_argument(
        "--dump-csv",
        type=str,
        default=None
    )
    parser.add_argument(
        "--dump-vis",
        default=False,
        action="store_true"
    )
    parser.add_argument(
        "--dump-weights",
        default=False,
        action="store_true"
    )
    parser.add_argument(
        "--dump-mode",
        choices=["vis-weight", "vis-only", "weight-only"],
        default="vis-weight"
    )
    return parser.parse_args(argv)


def dump_uvfits_data_to_csv(hdus, ts_limit=None, bl_limit=None, chan_limit=None, dump_mode="vis-weight", file=None):
    if file is None:
        file = stdout
    else:
        file = open(file, "w")
    timestep_baselines_seen = {}

    # Analyse row 0

    row_0_shape = hdus[0].data[0].data.shape

    num_chans = row_0_shape[2]
    num_pols = row_0_shape[3]
    assert num_pols == 4
    pol_names = ["xx", "yy", "xy", "yx"]
    if dump_mode == 'weight-only':
        pol_names = ["xx"]

    _chan_limit = num_chans
    if chan_limit is not None:
        _chan_limit = min(chan_limit, num_chans)

    header = [
        "timestep", 
        "baseline", 
        "pol",
    ]

    if dump_mode == 'vis-weight':
        header.extend(["type"])

    header.extend(map(str, range(_chan_limit)))

    print(", ".join(header), file=file)

    for row in hdus[0].data:

        timestep = row['DATE']
        if timestep not in timestep_baselines_seen:
            if ts_limit is not None and len(timestep_baselines_seen) >= ts_limit:
                break
            timestep_baselines_seen[timestep] = set()

        baselines_seen = timestep_baselines_seen[timestep]
        baseline = row['BASELINE']
        if baseline not in baselines_seen:
            if bl_limit is not None and len(baselines_seen) >= bl_limit:
                continue
            baselines_seen.add(baseline)

        # eprint(f"timestep: {timestep}, baseline: {baseline}")

        # eprint(f"timestep_baselines_seen ({len(timestep_baselines_seen)}): {timestep_baselines_seen}")

        row_data = row.data

        row_out = [
            timestep, 
            baseline, 
        ]
            
        for pol_idx, pol_name in enumerate(pol_names):
            _row_out = copy(row_out)
            _row_out.extend(map(str, [pol_name]))


            if dump_mode in ['vis-weight', 'vis-only']:
                row_real_data = row_data[:, :, :, pol_idx, 0].reshape((num_chans, ))
                row_imag_data = row_data[:, :, :, pol_idx, 1].reshape((num_chans, ))
                row_vis_data = row_real_data - 1j * row_imag_data

                _row_out_vis = copy(_row_out)
                if dump_mode == 'vis-weight':
                    _row_out_vis.extend(["vis"])
                _row_out_vis.extend(row_vis_data[:_chan_limit])

                print(", ".join(map(str, _row_out_vis)), file=file)

            if dump_mode in ['vis-weight', 'weight-only']:
                row_flag_data = row_data[:, :, :, pol_idx, 2].reshape((num_chans, ))
                _row_out_weight = copy(_row_out)
                if dump_mode == 'vis-weight':
                    _row_out_weight.extend(["weight"])
                _row_out_weight.extend(row_flag_data[:_chan_limit])

                print(", ".join(map(str, _row_out_weight)), file=file)

        else:
            for chan_idx in range(_chan_limit):

                row_out = copy(row_out)


def main(argv):
    args = parse_args(argv)

    hdus = fits.open(args.uvfile, memmap=True)
    print(f"-> hdus.info():")
    hdus.info()

    print("")
    print("VISIBILITIES")
    print("")

    print(repr(hdus[0].header))

    print(f"-> data.shape {hdus[0].data.shape}")

    dump_uvfits_data_to_csv(
        hdus, 
        ts_limit=args.timestep_limit, 
        bl_limit=args.baseline_limit, 
        chan_limit=args.chan_limit,
        dump_mode=args.dump_mode,
        file=args.dump_csv
    )

    # Antenna info
    print("")
    print("ANTENNA INFO")
    print("")
    if len(hdus) < 2:
        print("no antenna hdu")
        exit()
    print(repr(hdus[1].header))

    print(f"-> data ({hdus[1].data.shape})")
    print(tabulate(hdus[1].data, headers=hdus[1].data.names))


if __name__ == '__main__':
    main(sys.argv[1:])
