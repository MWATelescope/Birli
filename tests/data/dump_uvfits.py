#!/usr/bin/env python

from copy import copy
from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs, chdir
from typing import OrderedDict
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
import itertools
import cmath

from common import eprint


def decode_uvfits_baseline(bl: int):
    if bl < 65_535:
        ant2 = bl % 256
        ant1 = (bl - ant2) // 256
    else:
        ant2 = (bl - 65_536) % 2048
        ant1 = (bl - ant2 - 65_536) // 2048
    return (ant1, ant2)


POL_NAMES = ["xx", "yy", "xy", "yx"]

CC_AXIS_CHOICES = [
    "ts_frac",
    "ts_mjd",
    "bl_frac",
    "bl_idx",
    "ch_frac",
    "ch_idx",
    "u",
    "v",
    "w",
]
for pol in POL_NAMES:
    CC_AXIS_CHOICES.extend([
        f"{pol}_real",
        f"{pol}_imag",
        f"{pol}_mag",
        f"{pol}_phase",
        f"{pol}_weight",
    ])


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "uvfile",
        type=str
    )
    parser.add_argument(
        "--timestep-offset",
        help="number of timesteps to skip",
        type=int,
        default=None
    )
    parser.add_argument(
        "--timestep-limit",
        help="number of timesteps to process",
        type=int,
        default=None
    )
    parser.add_argument(
        "--baseline-offset",
        help="number of baselines to skip",
        type=int,
        default=None
    )
    parser.add_argument(
        "--baseline-limit",
        help="number of baselines to process",
        type=int,
        default=None
    )
    parser.add_argument(
        "--chan-limit",
        help="number of frequency channels to process",
        type=int,
        default=None
    )
    parser.add_argument(
        "--dump-csv",
        type=str,
        default=None
    )
    parser.add_argument(
        "--dump-mode",
        help="flavour of CSV dump",
        choices=["vis-weight", "vis-only", "weight-only", "cloud-compare"],
        default="vis-weight"
    )

    # used in cloud compare dump modes only
    cloud_compare_group = parser.add_argument_group(
        "cloud compare",
    )
    cloud_compare_group.add_argument(
        "--vis-peak",
        help="approximate value for normalising visibility amplitudes",
        type=float,
        default=None
    )
    cloud_compare_group.add_argument(
        "--start-mjd",
        help="start time [MJD] for normalising timesteps",
        default=None
    )
    cloud_compare_group.add_argument(
        "--end-mjd",
        help="end time [MJD] for normalising timesteps",
        default=None
    )
    cloud_compare_group.add_argument(
        "--invert-y",
        help="invert the y axis",
        default=False,
        action="store_true"
    )
    # cloud_compare_group.add_argument(
    #     "--x-axis",
    #     help="name of parameter used as x-axis",
    #     choices=CC_AXIS_CHOICES,
    #     default="ts_frac"
    # )
    # cloud_compare_group.add_argument(
    #     "--y-axis",
    #     help="name of parameter used as y-axis",
    #     choices=CC_AXIS_CHOICES,
    #     default="ch_frac"
    # )
    # cloud_compare_group.add_argument(
    #     "--z-axis",
    #     help="name of parameter used as z-axis",
    #     choices=CC_AXIS_CHOICES,
    #     default="xx_mag"
    # )
    cloud_compare_group.add_argument(
        "--axes",
        help="names of parameters used as cloud compare axes",
        choices=CC_AXIS_CHOICES,
        nargs="+",
        default=["ts_frac", "ch_frac", "xx_mag", "xx_phase"]
    )

    return parser.parse_args(argv)


def dump_uvfits_data_to_csv(
    hdus,
    ts_offset=None,
    ts_limit=None,
    mjd_start=None,
    mjd_end=None,
    bl_offset=None,
    bl_limit=None,
    chan_limit=None,
    vis_scale=1,
    dump_mode="vis-weight",
    axes=None,
    file=None,
    invert_y=False,
):
    if file is None:
        file = stdout
    else:
        file = open(file, "w")

    if bl_limit is not None and bl_offset is not None:
        bl_limit = bl_offset + bl_limit

    if mjd_start is not None:
        mjd_start = float(mjd_start)
    if mjd_end is not None:
        mjd_end = float(mjd_end)

    timestep_baselines_seen = OrderedDict()
    # Assume baselines are seen in order
    all_baselines_seen = []

    # Analyse row 0

    row_0_shape = hdus[0].data[0].data.shape

    num_chans = row_0_shape[2]
    num_pols = row_0_shape[3]
    assert num_pols == 4
    pol_names = POL_NAMES
    if dump_mode == 'weight-only':
        pol_names = ["xx"]

    _chan_limit = num_chans
    if chan_limit is not None:
        _chan_limit = min(chan_limit, num_chans)

    if dump_mode in ['vis-weight', 'vis-only', 'weight-only']:
        header = [
            "timestep",
            "baseline",
            "u",
            "v",
            "w",
            "pol",
        ]

        if dump_mode == 'vis-weight':
            header.extend(["type"])

        header.extend(map(str, range(_chan_limit)))

    elif dump_mode in ["cloud-compare"]:
        header = axes

    print(", ".join(header), file=file)
    eprint(f"-> vis scale {vis_scale}")

    timestep_fraction = None

    vis_peak = 0

    row_dict = {}

    for row in hdus[0].data:

        row_dict['ts_mjd'] = float(row['DATE'])
        if row_dict['ts_mjd'] not in timestep_baselines_seen:
            timestep_baselines_seen[row_dict['ts_mjd']] = set()

        row_dict['ts_idx'] = list(
            timestep_baselines_seen.keys()).index(row_dict['ts_mjd'])
        # eprint(f"row_dict['ts_mjd'] index {row_dict['ts_idx']}")
        if ts_limit is not None and row_dict['ts_idx'] >= ts_limit:
            break
        if ts_offset is not None and row_dict['ts_idx'] < ts_offset:
            continue

        baselines_seen = timestep_baselines_seen[row_dict['ts_mjd']]
        baseline = int(row['BASELINE'])
        if baseline not in baselines_seen:
            baselines_seen.add(baseline)
        if baseline not in all_baselines_seen:
            all_baselines_seen.append(baseline)

        row_dict['bl_idx'] = all_baselines_seen.index(baseline)
        if bl_limit is not None and row_dict['bl_idx'] >= bl_limit:
            continue
        if bl_offset is not None and row_dict['bl_idx'] < bl_offset:
            continue

        ant1_idx, ant2_idx = decode_uvfits_baseline(row_dict['bl_idx'])

        if mjd_start is not None and mjd_end is not None:
            row_dict['ts_frac'] = (
                row_dict['ts_mjd'] - mjd_start) / (mjd_end - mjd_start)
            eprint(
                f"timestep fraction {row_dict['ts_frac']:09f}, bl {row_dict['bl_idx']} {(ant1_idx, ant2_idx)}", end="\r")

        if dump_mode in ["cloud-compare"]:
            row_dict['ts_mjd'] = row_dict['ts_frac']
        else:
            row_dict['bl_idx'] = baseline

        row_dict['u'] = row['UU']
        row_dict['v'] = row['VV']
        row_dict['w'] = row['WW']

        row_data = row.data

        if dump_mode in ['vis-weight', 'vis-only', 'weight-only']:
            row_out = [
                row_dict['ts_mjd'],
                row_dict['bl_idx'],
                row_dict['u'],
                row_dict['v'],
                row_dict['w'],
            ]
            for pol_idx, pol_name in enumerate(pol_names):
                row_out_pol = copy(row_out)
                row_out_pol.extend([pol_name])

                if dump_mode in ['vis-weight', 'vis-only']:
                    row_real_data = row_data[:, :, :,
                                             pol_idx, 0].reshape((num_chans, ))
                    row_imag_data = row_data[:, :, :,
                                             pol_idx, 1].reshape((num_chans, ))
                    row_vis_data = row_real_data - 1j * row_imag_data
                    row_out_vis = copy(row_out_pol)

                    if dump_mode == 'vis-weight':
                        row_out_vis.extend(["vis"])
                    row_out_vis.extend(row_vis_data[:_chan_limit])

                    print(", ".join(map(str, row_out_vis)), file=file)

                    vis_peak = max(vis_peak, *map(abs, row_vis_data))

                if dump_mode in ['vis-weight', 'weight-only']:
                    row_weight_data = row_data[:, :, :,
                                               pol_idx, 2].reshape((num_chans, ))
                    row_out_weight = copy(row_out_pol)
                    if dump_mode == 'vis-weight':
                        row_out_weight.extend(["weight"])
                    row_out_weight.extend(row_weight_data[:_chan_limit])

                    print(", ".join(map(str, row_out_weight)), file=file)

        elif dump_mode in ["cloud-compare"]:
            for chan_idx in range(_chan_limit):

                row_dict["ch_idx"] = chan_idx
                row_dict["ch_frac"] = chan_idx / max(1, num_chans - 1)

                for pol_idx, pol_name in enumerate(pol_names):
                    # header.extend([f"{pol}_re", f"{pol}_im"])
                    vis = complex(
                        row_data[0, 0, chan_idx, pol_idx, 0] * vis_scale,
                        row_data[0, 0, chan_idx, pol_idx, 1] * -vis_scale
                    )
                    weight = row_data[0, 0, chan_idx, pol_idx, 2]

                    phase = cmath.phase(vis)

                    # normalise phase from (-pi,pi) to  (0,1)
                    if phase < 0:
                        phase += math.tau
                    phase /= math.tau

                    row_dict[f"{pol_name}_real"] = vis.real
                    row_dict[f"{pol_name}_imag"] = vis.imag
                    row_dict[f"{pol_name}_mag"] = abs(vis)
                    row_dict[f"{pol_name}_phase"] = phase
                    row_dict[f"{pol_name}_weight"] = weight

                    vis_peak = max(vis_peak, abs(vis))

                row_out = list(map(lambda ax: row_dict[ax], axes))
                if invert_y:
                    row_out[1] = -1 * row_out[1]
                print(", ".join(map(str, row_out)), file=file)

    print(f"aggregate info")
    print(f"-> vis peak {vis_peak}")
    mjds = timestep_baselines_seen.keys()
    print(f"-> {len(mjds)} mjds from {min (mjds)} to {max(mjds)}")


def main(argv):
    args = parse_args(argv)

    eprint(f"-> args {vars(args)}")

    hdus = fits.open(args.uvfile, memmap=True)
    print(f"-> hdus.info():")
    hdus.info()

    print("")
    print("VISIBILITIES")
    print("")

    print(repr(hdus[0].header))

    print(f"-> data.shape {hdus[0].data.shape}")

    if args.dump_csv is not None:
        kwargs = {}
        for key, arg in [
            ("ts_offset", args.timestep_offset),
            ("ts_limit", args.timestep_limit),
            ("mjd_start", args.start_mjd),
            ("mjd_end", args.end_mjd),
            ("bl_offset", args.baseline_offset),
            ("bl_limit", args.baseline_limit),
            ("chan_limit", args.chan_limit),
            ("invert_y", args.invert_y),
            ("axes", args.axes),
        ]:
            if arg is not None:
                kwargs[key] = arg
        if args.vis_peak is not None:
            kwargs["vis_scale"] = 1/args.vis_peak
        dump_uvfits_data_to_csv(
            hdus,
            dump_mode=args.dump_mode,
            file=args.dump_csv,
            **kwargs
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
