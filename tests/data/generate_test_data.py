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


DST_TEST_DIR_MWAX = "tests/data/1297526432_mwax"
SRC_TEST_DIR_MWAX = "/Volumes/T7/CIRA/1297526432_mwax"
TEST_METAFITS_NAME_MWAX = "1297526432.metafits"
TEST_GPUFITS_NAMES_MWAX = [
    "1297526432_20210216160014_ch117_000.fits",
    "1297526432_20210216160014_ch117_001.fits",
    "1297526432_20210216160014_ch118_000.fits",
    "1297526432_20210216160014_ch118_001.fits"
]
SRC_TEST_DIR_MWA_ORD = "/Volumes/T7/CIRA/1196175296_vis"
TEST_METAFITS_NAME_MWA_ORD = "1196175296.metafits"
TEST_GPUFITS_NAMES_MWA_ORD = [
    "1196175296_20171201145440_gpubox01_00.fits",
    "1196175296_20171201145540_gpubox01_01.fits",
    "1196175296_20171201145440_gpubox02_00.fits",
    "1196175296_20171201145540_gpubox02_01.fits"
]
RE_MWAX_NAME = (
    r"(?P<obsid>\d{10})_(?P<datetime>\d{8}(.)?\d{6})_ch(?P<rec_chan>\d{3})_(?P<batch>\d{3}).fits"
)
RE_MWA_ORD_NAME = (
    r"(?P<obsid>\d{10})_(?P<datetime>\d{14})_gpubox(?P<gpubox_num>\d{2})_(?P<batch>\d{2}).fits"
)

MAX_COARSE_CHANS = 2
MAX_BATCHES = 2
MAX_SCANS = 2
MAX_ANTENNAS = 2
MAX_FINE_CHANS = 2


def parse_filename(name, corr_type="MWAX", metafits_coarse_chans=[]):
    result = {}
    if corr_type == "MWAX":
        result = re.match(RE_MWAX_NAME, name).groupdict()
        result['corr_chan'] = metafits_coarse_chans.index(result['rec_chan'])
    elif corr_type == "MWA_ORD":
        result = re.match(RE_MWA_ORD_NAME, name).groupdict()
        result['corr_chan'] = int(result['gpubox_num']) - 1
        result['rec_chan'] = metafits_coarse_chans[result['corr_chan'] - 1]
    result['name'] = name
    return result


def chunk(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)


def get_input_df(tile_data):
    input_cols = ['Input', 'Antenna', 'Tile', 'TileName', 'Pol', 'Rx', 'Slot']
    return pd.DataFrame(dict([
        (col, list(tile_data.field(col)[:]))
        for col in input_cols
    ]))


def get_limited_antennas(input_df, max_antennas=None):
    antennas = sorted(list(np.unique(input_df.Antenna)))
    if max_antennas:
        antennas = antennas[:max_antennas]
    return antennas

def get_global_scan_index(channel_index, batch_index, max_batches, scan_index, max_scans):
    result = channel_index
    result = result << math.ceil( math.log2(max_batches)) | batch_index
    result = result << math.ceil( math.log2(max_scans)) | scan_index
    return result


def generate(args):
    def with_src_dir(n):
        return abspath(path_join(args['src_dir'], n))

    def with_dst_dir(n):
        return abspath(path_join(args['dst_dir'], n))

    ####
    # handle metafits
    ####

    metafits_path = with_src_dir(args['metafits_name'])
    dst_metafits_path = with_dst_dir(args['metafits_name'])
    print(metafits_path)
    assert exists(metafits_path), f"metafits_path {metafits_path} should exist"
    with fits.open(metafits_path) as meta_fits:
        primary_hdu = meta_fits[0]
        print(f" -> meta_fits.info()\n{pformat(meta_fits.info())}")
        print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")
        # print(f" -> meta_fits[1].header\n{repr(meta_fits[1].header)}")
        # print(f" -> meta_fits[1].columns\n{repr(meta_fits[1].columns)}")
        # print(f" -> meta_fits[1].data.shape\n{repr(meta_fits[1].data.shape)}")

        ####
        # Handle TILE_DATA: inputs, antennas, baselines
        ####

        tile_data = meta_fits[1].data

        input_df = get_input_df(tile_data)
        print(f" -> input_df:\n{input_df}")

        # Limit number of antennas
        antennas = get_limited_antennas(input_df, args.get('max_antennas'))
        num_antennas = len(antennas)
        print(f" -> antennas({num_antennas}):\n{antennas}")

        # Calculate baselines
        num_baselines = (num_antennas) * (num_antennas + 1) // 2
        print(f" -> num_baselines: {num_baselines}")

        filtered_input_df = input_df[np.isin(input_df.Antenna, antennas)]
        num_inputs = len(filtered_input_df)
        print(f" -> filtered_inputs({num_inputs}):\n{filtered_input_df}")


        tile_data = tile_data[np.isin(tile_data.Antenna, antennas)]

        primary_hdu.header['NINPUTS'] = len(tile_data)
        meta_fits[1].data = tile_data

        ####
        # Handle TILE_DATA: polarizations
        ####

        ant_pols = np.unique(input_df.Pol)
        print(f" -> ant_pols({len(ant_pols)}):\n{ant_pols}")
        num_corr_pols = len(ant_pols) ** 2
        print(f" -> num_corr_pols: {num_corr_pols}")

        ####
        # Handle TILE_DATA: receivers, delays
        ####

        receiver_delays = dict(zip(
            primary_hdu.header['RECVRS'].split(","),
            primary_hdu.header['DELAYS'].split(",")
        ))
        print(f" -> receiver_delays:\n{receiver_delays}")
        valid_receivers = list(np.unique(filtered_input_df.Rx))
        print(f" -> valid_receivers:\n{valid_receivers}")
        valid_reciever_delays = dict(filter(
            lambda rx_del: int(rx_del[0]) in valid_receivers,
            receiver_delays.items()))

        primary_hdu.header['RECVRS'] = ",".join(valid_reciever_delays.keys())
        primary_hdu.header['DELAYS'] = ",".join(valid_reciever_delays.values())


        ####
        # Handle Coarse Channels
        ####

        metafits_coarse_chans = primary_hdu.header['CHANNELS'].split(",")
        num_coarse_chans = len(metafits_coarse_chans)
        print(f" -> metafits_coarse_chans({num_coarse_chans}):\n{metafits_coarse_chans}")

        ####
        # Handle Coarse/Fine Channel Bandwidth
        ####

        fine_chan_bandwidth_hz = int(primary_hdu.header['FINECHAN'] * 1_000)  # Header in KHz
        print(f" -> fine_chan_bandwidth (Hz): {fine_chan_bandwidth_hz}")
        total_chan_bandwidth_hz = int(primary_hdu.header['BANDWDTH'] * 1_000_000)  # Header in MHz
        coarse_chan_bandwidth_hz = total_chan_bandwidth_hz // num_coarse_chans
        print(f" -> coarse_chan_bandwidth (Hz, derived): {coarse_chan_bandwidth_hz}")
        num_fine_chans = coarse_chan_bandwidth_hz // fine_chan_bandwidth_hz
        print(f" -> num_fine_chans (derived): {num_fine_chans}")
        num_fine_chans = min(args['max_fine_chans'], num_fine_chans)
        print(f" -> num_fine_chans (limited): {num_fine_chans}")
        fine_chan_bandwidth_hz = coarse_chan_bandwidth_hz // num_fine_chans
        print(f" -> fine_chan_bandwidth (Hz, limited): {fine_chan_bandwidth_hz}")
        coarse_chan_bandwidth_hz = fine_chan_bandwidth_hz * num_fine_chans
        print(f" -> coarse_chan_bandwidth (Hz, limited): {coarse_chan_bandwidth_hz}")
        primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000

        ####
        # Handle Fine Channel Selection
        ####

        coarse_channel_selection = list(filter(None, primary_hdu.header['CHANSEL'].split(',')))
        print(f" -> coarse_channel_selection (derived): {coarse_channel_selection}")
        coarse_channel_selection = coarse_channel_selection[:args['max_coarse_chans']]
        print(f" -> coarse_channel_selection (limited): {coarse_channel_selection}")
        primary_hdu.header['CHANSEL'] = ','.join(coarse_channel_selection)

        ####
        # handle NAXIS*
        ####

        floats_per_complex = 2  # [Real, Imag]

        if args['corr_type'] == "MWAX":
            naxis1 = num_fine_chans * num_corr_pols * floats_per_complex
            naxis2 = num_baselines
        elif args['corr_type'] == "MWA_ORD":
            naxis1 = num_corr_pols * num_baselines * floats_per_complex
            naxis2 = num_fine_chans

        print(f" -> naxis1: {naxis1}")
        print(f" -> naxis2: {naxis2}")
        floats_per_img = naxis1 * naxis2
        print(f" -> floats_per_img: {floats_per_img}")

        ####
        # Handle Scans
        ####

        primary_hdu.header['NSCANS'] = min(args['max_scans'], primary_hdu.header['NSCANS'])

        ####
        # Analyse GPUFits
        ####

        gpufits_df = pd.DataFrame([
            parse_filename(gpubox_name, args['corr_type'], metafits_coarse_chans)
            for gpubox_name in args['gpufits_names']
        ])
        print(f" -> gpufits_df:\n{gpufits_df}")

        rec_coarse_chans = np.unique(gpufits_df.rec_chan)
        print(f" -> rec_coarse_chans({len(rec_coarse_chans)}):\n{rec_coarse_chans}")
        valid_coarse_chans = sorted(rec_coarse_chans)
        if 'max_coarse_chans' in args:
            valid_coarse_chans = valid_coarse_chans[:args['max_coarse_chans']]
        num_coarse_chans = len(valid_coarse_chans)
        print(f" -> valid_coarse_chans({num_coarse_chans}):\n{valid_coarse_chans}")

        primary_hdu.header['CHANNELS'] = ','.join(f'{chan}' for chan in valid_coarse_chans)
        primary_hdu.header['CENTCHAN'] = valid_coarse_chans[num_coarse_chans//2]
        primary_hdu.header['NCHANS'] = num_coarse_chans * args['max_fine_chans']

        total_chan_bandwidth_hz = coarse_chan_bandwidth_hz * num_coarse_chans
        primary_hdu.header['BANDWDTH'] = total_chan_bandwidth_hz / 1_000_000

        gpufits_batches = sorted(list(np.unique(gpufits_df.batch)[:args['max_batches']]))
        print(f" -> gpufits_batches({len(gpufits_batches)}):\n{gpufits_batches}")

        filtered_gpufits = gpufits_df\
            [np.isin(gpufits_df.rec_chan, valid_coarse_chans)]\
            [np.isin(gpufits_df.batch, gpufits_batches)]\
            # .sort_values(['corr_chan', 'batch'])
        print(f" -> filtered_gpufits:\n{filtered_gpufits}")

        ####
        # Show modifications
        ####

        print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")

        ####
        # write MetaFits
        ####

        if not path_exists(dirname(dst_metafits_path)):
            makedirs(dirname(dst_metafits_path))

        meta_fits.writeto(dst_metafits_path, overwrite=True)

    ####
    # Handle GPUFits
    ####

    float_count = 0

    for row in filtered_gpufits.to_dict('records'):
        gpufits_name = row['name']
        print(f" -> gpufits_name: {gpufits_name}")
        # channel_index = row['corr_chan']
        channel_index = valid_coarse_chans.index(row['rec_chan'])
        print(f" -> channel_index: {channel_index:08b}")
        batch_index = gpufits_batches.index(row['batch'])
        gpufits_path = with_src_dir(gpufits_name)
        dst_gpufits_path = with_dst_dir(gpufits_name)
        with fits.open(gpufits_path) as gpu_fits:
            # print(f" -> gpu_fits[0].info()\n{pformat(gpu_fits.info())}")
            # print(f" -> gpu_fits[0].header\n{repr(gpu_fits[0].header)}")
            # print(f" -> gpu_fits[1].header\n{repr(gpu_fits[1].header)}")
            # print(f" -> gpu_fits[2].header\n{repr(gpu_fits[2].header)}")

            primary_hdu = gpu_fits[0]
            primary_hdu.header['CORR_VER'] = 2
            primary_hdu.header['NFINECHS'] = num_fine_chans
            primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000
            primary_hdu.header['NINPUTS'] = num_inputs
            scan_hdus = []

            if args['corr_type'] == "MWAX":

                scan_hdu_chunks = chunk(gpu_fits[1:][:args['max_scans']*2], 2)
                for (scan_index, (img_hdu, flag_hdu)) in enumerate(scan_hdu_chunks):
                    global_scan_index = get_global_scan_index(channel_index, batch_index, args['max_batches'], scan_index, args['max_scans'])
                    # prefix scan indices to start with the char 'A'
                    global_scan_index = (0x41 << 8) | global_scan_index
                    print(f" -> global_scan_index: {global_scan_index:08b}")

                    print(f" -> img_hdu[{scan_index}].data ({img_hdu.data.shape}, {img_hdu.data.dtype}): \n{img_hdu.data}")
                    print(f" -> flag_hdu[{scan_index}].data.shape({flag_hdu.data.shape}, {flag_hdu.data.dtype})")
                    img_hdu.header['NAXIS1'] = naxis1
                    img_hdu.header['NAXIS2'] = naxis2

                    float_start = global_scan_index << 8 * math.ceil(math.log2(floats_per_img) / 8)
                    print(f" -> float_start: {float_start} ({float_start:032b}, ~2^{math.log2(float_start)})")
                    float_end = float_start + floats_per_img
                    print(f" -> float_end: {float_end} ({float_end:032b}, ~2^{math.log2(float_end)})")
                    img_hdu.data = np.arange(float_start, float_end).reshape((naxis2, naxis1)).astype(np.int32)
                    flag_hdu.data = flag_hdu.data[:naxis2, :]

                    scan_hdus.append(img_hdu)
                    scan_hdus.append(flag_hdu)

            elif args['corr_type'] == "MWA_ORD":
                for (scan_index, img_hdu) in enumerate(gpu_fits[1:][:args['max_scans']]):
                    global_scan_index = get_global_scan_index(channel_index, batch_index, args['max_batches'], scan_index, args['max_scans'])
                    # global_scan_index = (0x1001 << 8) | global_scan_index
                    print(f" -> global_scan_index: {global_scan_index:08b}")

                    print(f" -> img_hdu[{scan_index}].data ({img_hdu.data.shape}, {img_hdu.data.dtype}): \n{img_hdu.data}")

                    img_hdu.header['NAXIS1'] = naxis1
                    img_hdu.header['NAXIS2'] = naxis2

                    if args.get('sequential'):
                        float_start = float_count
                    else:
                        float_start = global_scan_index << math.ceil(math.log2(floats_per_img))
                    print(f" -> float_start: {float_start} ({float_start:064b}, ")
                    float_end = float_start + floats_per_img
                    print(f" -> float_end: {float_end} ({float_end:064b}, ~2^{math.log2(float_end)})")
                    img_hdu.data = np.arange(float_start, float_end).reshape((naxis2, naxis1)).astype(np.float64)

                    float_count += floats_per_img
                    scan_hdus.append(img_hdu)

            new_gpu_fits = fits.HDUList([primary_hdu] + scan_hdus)

            if not path_exists(dirname(dst_gpufits_path)):
                makedirs(dirname(dst_gpufits_path))

            print(f'-> writing to {dst_gpufits_path}')
            new_gpu_fits.writeto(dst_gpufits_path, overwrite=True)


def main():
    generate({
        'corr_type': "MWA_ORD",
        'src_dir': '/Users/derwent/Documents/CIRA/code/mwalib/test_files/1101503312_1_timestep/',
        'dst_dir': "tests/data/1101503312_mwa_ord",
        'metafits_name': '1101503312.metafits',
        'gpufits_names': ['1101503312_20141201210818_gpubox01_00.fits'],
        'max_coarse_chans': MAX_COARSE_CHANS,
        'max_batches': MAX_BATCHES,
        'max_scans': MAX_SCANS,
        # 'max_antennas': MAX_ANTENNAS,
        'max_fine_chans': MAX_FINE_CHANS,
    })
    # # cargo run dump-all-data \
    # #     --dump-filename=../Birli/tests/data/1101503312_mwa_ord/1101503312_dump.csv \
    # #     --metafits=../Birli/tests/data/1101503312_mwa_ord/1101503312.metafits \
    # #     ../Birli/tests/data/1101503312_mwa_ord/1101503312_*.fits
    generate({
        'corr_type': "MWAX",
        'src_dir': SRC_TEST_DIR_MWAX,
        'dst_dir': "tests/data/1297526432_mwax",
        'metafits_name': TEST_METAFITS_NAME_MWAX,
        'gpufits_names': TEST_GPUFITS_NAMES_MWAX,
        'max_coarse_chans': MAX_COARSE_CHANS,
        'max_batches': MAX_BATCHES,
        'max_scans': MAX_SCANS,
        'max_antennas': MAX_ANTENNAS,
        'max_fine_chans': MAX_FINE_CHANS,
    })
    # cargo run dump-all-data \
    #   --dump-filename=../Birli/tests/data/1297526432_mwax/1297526432_dump.csv \
    #   --metafits=../Birli/tests/data/1297526432_mwax/1297526432.metafits \
    #   ../Birli/tests/data/1297526432_mwax/1297526432_20210216160014_ch*.fits
    generate({
        'corr_type': "MWA_ORD",
        'src_dir': SRC_TEST_DIR_MWA_ORD,
        'dst_dir': "tests/data/1196175296_mwa_ord",
        'metafits_name': TEST_METAFITS_NAME_MWA_ORD,
        'gpufits_names': TEST_GPUFITS_NAMES_MWA_ORD,
        'max_coarse_chans': MAX_COARSE_CHANS,
        'max_batches': MAX_BATCHES,
        'max_scans': MAX_SCANS,
        # 'max_antennas': MAX_ANTENNAS,
        'max_fine_chans': MAX_FINE_CHANS,
    })
    # cargo run dump-all-data --vis-radix=16 --absolute \
    #     --dump-filename=../Birli/tests/data/1196175296_mwa_ord/1196175296_dump_hex.csv \
    #     --metafits=../Birli/tests/data/1196175296_mwa_ord/1196175296.metafits \
    #     ../Birli/tests/data/1196175296_mwa_ord/1196175296_*.fits | tee dump_out.txt
    # cargo run dump-all-data \
    #   --dump-filename=../Birli/tests/data/1196175296_mwa_ord/1196175296_dump.csv \
    #   --metafits=../Birli/tests/data/1196175296_mwa_ord/1196175296.metafits \
    #   ../Birli/tests/data/1196175296_mwa_ord/1196175296_*.fits
    generate({
        'corr_type': "MWA_ORD",
        'src_dir': SRC_TEST_DIR_MWA_ORD,
        'dst_dir': "tests/data/1196175296_mwa_ord_seq",
        'metafits_name': TEST_METAFITS_NAME_MWA_ORD,
        'gpufits_names': TEST_GPUFITS_NAMES_MWA_ORD,
        'max_coarse_chans': MAX_COARSE_CHANS,
        'max_batches': MAX_BATCHES,
        'max_scans': MAX_SCANS,
        # 'max_antennas': MAX_ANTENNAS,
        'max_fine_chans': MAX_FINE_CHANS,
        'sequential': True,
    })
    # # cargo run dump-all-data \
    # #   --dump-filename=../Birli/tests/data/1196175296_mwa_ord_seq/1196175296_dump.csv \
    # #   --metafits=../Birli/tests/data/1196175296_mwa_ord_seq/1196175296.metafits \
    # #   ../Birli/tests/data/1196175296_mwa_ord_seq/1196175296_*.fits


if __name__ == '__main__':
    main()
