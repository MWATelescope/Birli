#!/usr/bin/env python

from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs

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
from tabulate import tabulate

import numpy as np


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "file",
        type=str
    )
    return parser.parse_args(argv)

def main(argv):
    args = parse_args(argv)
    hdus = fits.open(args.file)
    print(f"-> hdus.info():")
    hdus.info()

    print("")
    print("HEADER")
    print("")

    primary_hdu = hdus[0]
    print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")
    # hdus[0].info()

    ####
    # Tile Data
    ####

    tile_header = hdus[1].data.dtype.names
    tile_data = hdus[1].data
    print(tabulate(tile_data, headers=tile_header))

    # norths = sorted(list(tile_data['North']))
    # easts = sorted(list(tile_data['East']))
    # heights = sorted(list(tile_data['Height']))

    # print(f" -> norths range: {norths[0]}..{norths[-1]}")
    # print(f" -> easts range: {easts[0]}..{easts[-1]}")
    # print(f" -> heights range: {heights[0]}..{heights[-1]}")

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')

    # ax.plot3D((tile_data['North'], tile_data['East'], tile_data['Height'], 'gray'))

    # fig.show()

    # pass

if __name__ == '__main__':
    main(sys.argv[1:])
    # main(["tests/data/1196175296_mwa_ord/1196175296.metafits"])

