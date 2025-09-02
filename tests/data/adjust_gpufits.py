#!/usr/bin/env python

from os.path import dirname
from os.path import exists as path_exists
from os import makedirs

from astropy.io import fits
import numpy as np
import sys
from argparse import ArgumentParser

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
        default=None
    )
    parser.add_argument(
        "--timestep-offset",
        help="the number of scans to offset these visibilities",
        type=int,
        default=None
    )
    parser.add_argument(
        "--empty-data",
        help="simulate empty data",
        action='store_true'
    )
    return parser.parse_args(argv)


def offset_hdu(hdu, offset):
    time = hdu.header['TIME'] + hdu.header['MILLITIM'] / 1000 + offset
    hdu.header['TIME'] = int(time)
    hdu.header['MILLITIM'] = int(1000 * (time % 1))


def empty_data(hdu):
    hdu.header['NAXIS1'] = 1
    hdu.header['NAXIS2'] = 1
    hdu.data = np.fromiter([0], dtype=np.float64).reshape((1, 1))


def main(argv):
    args = parse_args(argv)

    hdus = fits.open(args.in_file)
    print("-> hdus.info():")
    hdus.info()

    print("")
    print("HEADER")
    print("")

    primary_hdu = hdus[0]
    print(repr(primary_hdu.header))

    time_offset = 0

    if args.timestep_offset is not None:
        int_time = primary_hdu.header['INTTIME']
        time_offset = int_time * args.timestep_offset
        offset_hdu(primary_hdu, time_offset)

    timestep_limit = get_gpufits_num_scans(len(hdus), args.corr_type)
    if args.timestep_limit is not None:
        timestep_limit = min(args.timestep_limit, timestep_limit)

    print("")
    print(f"VISIBILITIES ({timestep_limit}):")
    print("")

    scan_hdus = []

    if args.corr_type == "MWAX":
        scan_hdu_chunks = chunk(hdus[1:][:timestep_limit * 2], 2)
        for (_, (img_hdu, flag_hdu)) in enumerate(scan_hdu_chunks):

            if time_offset:
                offset_hdu(img_hdu, time_offset)
                offset_hdu(flag_hdu, time_offset)

            if args.empty_data:
                empty_data(img_hdu)
                empty_data(flag_hdu)

            scan_hdus.append(img_hdu)
            scan_hdus.append(flag_hdu)

    elif args.corr_type == "MWA_ORD":
        for (_, img_hdu) in enumerate(hdus[1:][:timestep_limit]):
            if time_offset:
                offset_hdu(img_hdu, time_offset)

            if args.empty_data:
                empty_data(img_hdu)

            scan_hdus.append(img_hdu)

    new_gpu_fits = fits.HDUList([primary_hdu] + scan_hdus)

    print(f'-> writing {len(scan_hdus)} scans to {args.out_file}')
    if (parent := dirname(args.out_file)) and not path_exists(parent):
        makedirs(parent)
    new_gpu_fits.writeto(args.out_file, overwrite=True)


if __name__ == '__main__':
    main(sys.argv[1:])
