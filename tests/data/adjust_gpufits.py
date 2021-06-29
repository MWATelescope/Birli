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

from common import get_gpufits_num_scans, chunk


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "--in-file",
        type=str
    )
    parser.add_argument(
        "--out-file",
        type=str
    )
    parser.add_argument(
        "--corr-type",
        help="the type of correlator",
        choices=["MWA_ORD", "MWAX"]
    )
    parser.add_argument(
        "--timestep-limit",
        type=int,
        default=0
    )
    parser.add_argument(
        "--timestep-offset",
        help="the number of scans to offset these visibilities",
        type=int,
        default=0
    )
    return parser.parse_args(argv)

def offset_hdu(hdu, offset):
    time = hdu.header['TIME'] + hdu.header['MILLITIM'] / 1000 + offset
    hdu.header['TIME'] = int(time)
    hdu.header['MILLITIM'] = int(1000 * (time % 1))


def main(argv):
    args = parse_args(argv)

    hdus = fits.open(args.in_file)
    print(f"-> hdus.info():")
    hdus.info()

    print("")
    print("HEADER")
    print("")

    primary_hdu = hdus[0]
    print(repr(primary_hdu.header))

    time_offset = 0

    if args.timestep_offset:
        int_time = primary_hdu.header['INTTIME']
        time_offset = int_time * args.timestep_offset
        offset_hdu(primary_hdu, time_offset)

    print("")
    print("VISIBILITIES")
    print("")

    timestep_limit = get_gpufits_num_scans(len(hdus), args.corr_type)
    if args.timestep_limit:
        timestep_limit = min(args.timestep_limit, timestep_limit)

    scan_hdus = []

    if args.corr_type == "MWAX":
        scan_hdu_chunks = chunk(hdus[1:][:timestep_limit*2], 2)
        for (_, (img_hdu, flag_hdu)) in enumerate(scan_hdu_chunks):

            if time_offset:
                offset_hdu(img_hdu, time_offset)
                offset_hdu(flag_hdu, time_offset)

            scan_hdus.append(img_hdu)
            scan_hdus.append(flag_hdu)

    elif args.corr_type == "MWA_ORD":
        for (_, img_hdu) in enumerate(hdus[1:][:timestep_limit]):
            if time_offset:
                offset_hdu(img_hdu, time_offset)

            scan_hdus.append(img_hdu)

    new_gpu_fits = fits.HDUList([primary_hdu] + scan_hdus)

    print(f'-> writing {len(scan_hdus)} scans to {args.out_file}')
    if not path_exists(dirname(args.out_file)):
        makedirs(dirname(args.out_file))
    new_gpu_fits.writeto(args.out_file, overwrite=True)

if __name__ == '__main__':
    main(sys.argv[1:])