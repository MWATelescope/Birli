#!/usr/bin/env python

from astropy.io import fits
import sys
from argparse import ArgumentParser

from common import get_gpufits_num_scans, chunk


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "gpufile",
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
        "--baseline-limit",
        type=int,
        default=0
    )
    return parser.parse_args(argv)


def dump_img_hdu(scan_idx, img_hdu):

    print("")
    print(f" -> SCAN {scan_idx} IMAGE")
    print("")

    print(repr(img_hdu.header))


def dump_flag_hdu(scan_idx, flag_hdu):

    print("")
    print(f" -> SCAN {scan_idx} FLAG")
    print("")

    print(repr(flag_hdu.header))


def main(argv):
    args = parse_args(argv)

    hdus = fits.open(args.gpufile)
    print("-> hdus.info():")
    hdus.info()

    print("")
    print("HEADER")
    print("")

    print(repr(hdus[0].header))

    timestep_limit = num_scans = get_gpufits_num_scans(len(hdus), args.corr_type)
    if args.timestep_limit:
        timestep_limit = min(args.timestep_limit, timestep_limit)

    print("")
    print(f"VISIBILITIES ({num_scans}):")
    print("")

    if args.corr_type == "MWAX":
        scan_hdu_chunks = chunk(hdus[1:][:timestep_limit * 2], 2)
        for (scan_idx, (img_hdu, flag_hdu)) in enumerate(scan_hdu_chunks):

            dump_img_hdu(scan_idx, img_hdu)
            dump_flag_hdu(scan_idx, flag_hdu)

    elif args.corr_type == "MWA_ORD":
        for (scan_idx, img_hdu) in enumerate(hdus[1:][:timestep_limit]):

            dump_img_hdu(scan_idx, img_hdu)

    exit()


if __name__ == '__main__':
    main(sys.argv[1:])
