#!/usr/bin/env python

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
from argparse import ArgumentParser
import pyuvdata
from pyuvdata import UVData
from tabulate import tabulate


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "uvfile",
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
    parser.add_argument(
        "--antenna-limit",
        type=int,
        default=0
    )
    return parser.parse_args(argv)

def main(argv):
    args = parse_args(argv)

    # pyuv = UVData()
    # pyuv.read_uvfits(args.uvfile)

    # antpairpols = pyuv.get_antpairpols()
    # for key in ['Nants_data', 'Nants_telescope', 'Nbls', 'Nfreqs', 'Npols', 'Ntimes']:
    #     print(f"-> {key}: {pyuv.__getattribute__(key)}")
    

    # for key, data in pyuv.antpairpol_iter():
    #     (ant1, ant2, pol) = key
    #     if args.antenna_limit:
    #         if ant1 >= args.antenna_limit or ant2 >= args.antenna_limit:
    #             continue
    #     print(f"-> ant1 {ant1}, ant2 {ant2}, pol {pol}")
    #     print(f" --> data \n{pd.DataFrame(data)}")
    #     flags = pyuv.get_flags(key)
    #     print(f" --> flags \n{pd.DataFrame(flags)}")


    hdus = fits.open(args.uvfile)
    print(f"-> hdus.info():")
    hdus.info()


    print("")
    print("VISIBILITIES")
    print("")

    print(repr(hdus[0].header))

    print(f"-> data.shape {hdus[0].data.shape}")

    baselines = np.unique(hdus[0].data.par("BASELINE"))
    num_baselines = len(baselines)
    # assert(pyuv.Nbls == num_baselines)
    baseline_limit = num_baselines
    if args.baseline_limit:
        if args.baseline_limit < baseline_limit:
            baseline_limit = args.baseline_limit
            show_baseline_ellipses = True
    timesteps = np.unique(hdus[0].data.par("DATE"))
    num_timesteps = len(timesteps)
    timestep_limit = num_timesteps
    if args.timestep_limit:
        if args.timestep_limit < timestep_limit:
            timestep_limit = args.timestep_limit
            show_timestep_ellipses = True

    # assert(pyuv.Ntimes == num_timesteps)
    assert(hdus[0].data.shape[0] == num_baselines * num_timesteps)
    # blt_data = np.reshape(np.array(hdus[0].data), (num_baselines, num_timesteps))
    # for bl in baselines[:baseline_limit]:
    #     for ts in timestep[:timestep_limit]:
    blt_data = hdus[0].data[::num_baselines]
    for bl_idx in range(baseline_limit):
        print(bl_idx)
        print(f"baseline idx {bl_idx}")
        bl_data = hdus[0].data[bl_idx::num_baselines]
        # vis_data = bl_data[:, :, :, :, :, 0] - 1j * bl_data[:, :, :, :, :, 1]
        # print(f"--> vis_data ({vis_data.shape}):")
        # print(vis_data)

        # flag_data = row_data[:, :, :, :, :, 2]
        # print(f"--> flag_data ({flag_data.shape}):")
        # print(flag_data)
        # pass
        vis_data = None
        flag_data = None
        for row in blt_data[:timestep_limit]:
            row_data = row['data']
            print(f"--> row_data {(row_data)}")

            # assert(pyuv.Nfreqs == row_data.shape[2])
            # assert(pyuv.Npols == row_data.shape[3])

            row_real_data = row_data[:, :, :, :, 0].reshape((row_data.shape[2], row_data.shape[3]))
            row_imag_data = row_data[:, :, :, :, 1].reshape((row_data.shape[2], row_data.shape[3]))
            row_flag_data = row_data[:, :, :, :, 2].reshape((row_data.shape[2], row_data.shape[3]))
            row_vis_data = row_real_data - 1j *  row_imag_data
            
            if vis_data is None:
                print(f"--> baseline {row['BASELINE']}")
                print(f"--> uvw {row[0:3]}")
                print(f"--> row_data shape {row_data.shape}")
                print(f"--> row_vis_data shape {row_vis_data.shape}")
                vis_data = row_vis_data.reshape(row_vis_data.shape + (1, ))
                flag_data = row_flag_data.reshape(row_flag_data.shape + (1,))
            else:
                vis_data = np.append(vis_data, row_vis_data.reshape(row_vis_data.shape + (1, )), axis=2)
                flag_data = np.append(flag_data, row_flag_data.reshape(row_flag_data.shape + (1,)), axis=2)

            
        for pol in range(vis_data.shape[1]):
            print(f"--> pol {pol}")
            print(f"---> vis_data ({vis_data.shape}):")
            print(pd.DataFrame(vis_data[:, pol, :].reshape(vis_data.shape[0], vis_data.shape[2]).transpose()))

            print(f"---> flag_data ({flag_data.shape}):")
            print(pd.DataFrame(flag_data[:, pol, :].reshape(flag_data.shape[0], flag_data.shape[2]).transpose()))
    if show_baseline_ellipses:
        print("...")


    # Antenna info
    print("")
    print("ANTENNA INFO")
    print("")
    
    print(repr(hdus[1].header))

    print(f"-> data ({hdus[1].data.shape})")
    print(tabulate(hdus[1].data, headers=hdus[1].data.names))


if __name__ == '__main__':
    main(sys.argv[1:])