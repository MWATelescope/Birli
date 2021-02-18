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


SRC_TEST_DIR = "/Volumes/T7/CIRA/1297526432_mwax"
DST_TEST_DIR = "tests/data/1297526432_mwax"
TEST_METAFITS_NAME = "1297526432.metafits"
TEST_GPUFITS_NAMES = [
    "1297526432_20210216160014_ch117_000.fits",
    "1297526432_20210216160014_ch117_001.fits",
    "1297526432_20210216160014_ch118_000.fits",
    "1297526432_20210216160014_ch118_001.fits"
]
RE_MWAX_NAME = (
    r"(?P<obsid>\d{10})_(?P<datetime>\d{8}(.)?\d{6})_ch(?P<channel>\d{3})_(?P<batch>\d{3}).fits")

MAX_COARSE_CHANS = 2
MAX_BATCHES = 2
MAX_SCANS = 2
MAX_ANTENNAS = 2
MAX_FINE_CHANS = 2


def parse_filename(name):
    return re.match(RE_MWAX_NAME, name).groupdict()


def chunk(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)


def main():
    def with_src_dir(n):
        return abspath(path_join(SRC_TEST_DIR, n))

    def with_dst_dir(n):
        return abspath(path_join(DST_TEST_DIR, n))

    ####
    # handle metafits
    ####

    metafits_path = with_src_dir(TEST_METAFITS_NAME)
    dst_metafits_path = with_dst_dir(TEST_METAFITS_NAME)
    print(metafits_path)
    assert exists(metafits_path)
    with fits.open(metafits_path) as meta_fits:
        primary_hdu = meta_fits[0]
        print(f" -> meta_fits.info()\n{pformat(meta_fits.info())}")
        print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")
        # print(f" -> meta_fits[1].header\n{repr(meta_fits[1].header)}")
        # print(f" -> meta_fits[1].columns\n{repr(meta_fits[1].columns)}")
        # print(f" -> meta_fits[1].data.shape\n{repr(meta_fits[1].data.shape)}")

        ####
        # handle inputs
        ####

        inputs = meta_fits[1].data

        input_cols = ['Input', 'Antenna', 'Tile', 'TileName', 'Pol', 'Rx', 'Slot']
        input_df = pd.DataFrame(dict([
            (col, list(inputs.field(col)[:]))
            for col in input_cols
        ]))
        print(f" -> inputs:\n{input_df}")

        ant_pols = np.unique(inputs.field('Pol'))
        print(f" -> ant_pols({len(ant_pols)}):\n{ant_pols}")
        num_corr_pols = len(ant_pols) ** 2
        print(f" -> num_corr_pols: {num_corr_pols}")

        antennas = sorted(list(np.unique(inputs.field('Antenna'))))[:MAX_ANTENNAS]
        num_antennas = len(antennas)
        print(f" -> antennas({num_antennas}):\n{antennas}")
        num_baselines = (num_antennas) * (num_antennas + 1) // 2
        print(f" -> num_baselines: {num_baselines}")

        filtered_input_df = input_df[np.isin(input_df.Antenna, antennas)]
        print(f" -> filtered_inputs({len(filtered_input_df)}):\n{filtered_input_df}")
        valid_inputs = list(filtered_input_df.Input)

        filtered_inputs = inputs[np.isin(inputs['Input'], valid_inputs)]
        primary_hdu.header['NINPUTS'] = len(filtered_inputs)
        meta_fits[1].data = filtered_inputs

        ####
        # handle receivers, delays
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
        # Handle Channels
        ####

        fine_chan_bandwidth_hz = int(primary_hdu.header['FINECHAN'] * 1_000)  # Header in KHz
        print(f" -> fine_chan_bandwidth (Hz): {fine_chan_bandwidth_hz}")
        total_chan_bandwidth_hz = int(primary_hdu.header['BANDWDTH'] * 1_000_000)  # Header in MHz
        coarse_chans = primary_hdu.header['CHANNELS'].split(",")
        num_coarse_chans = len(coarse_chans)
        print(f" -> coarse_chans({num_coarse_chans}):\n{coarse_chans}")
        coarse_chan_bandwidth_hz = total_chan_bandwidth_hz // num_coarse_chans
        print(f" -> coarse_chan_bandwidth (Hz, derived): {coarse_chan_bandwidth_hz}")
        num_fine_chans = coarse_chan_bandwidth_hz // fine_chan_bandwidth_hz
        print(f" -> num_fine_chans (derived): {num_fine_chans}")
        num_fine_chans = min(MAX_FINE_CHANS, num_fine_chans)
        print(f" -> num_fine_chans (limited): {num_fine_chans}")
        fine_chan_bandwidth_hz = coarse_chan_bandwidth_hz // num_fine_chans
        print(f" -> fine_chan_bandwidth (Hz, limited): {fine_chan_bandwidth_hz}")
        primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000
        coarse_chan_bandwidth_hz = fine_chan_bandwidth_hz * num_fine_chans
        print(f" -> coarse_chan_bandwidth (Hz, limited): {coarse_chan_bandwidth_hz}")

        fine_channel_selection = list(filter(None, primary_hdu.header['CHANSEL'].split(',')))
        print(f" -> fine_channel_selection (derived): {fine_channel_selection}")
        fine_channel_selection = fine_channel_selection[:MAX_FINE_CHANS]
        print(f" -> fine_channel_selection (limited): {fine_channel_selection}")
        primary_hdu.header['CHANSEL'] = ','.join(fine_channel_selection)

        ####
        # handle NAXIS*
        ####

        floats_per_complex = 2  # [Real, Imag]

        # This is the legacy calculation.
        # naxis1 = num_corr_pols * num_baselines * floats_per_complex

        # This is the MWAX calculation
        naxis1 = num_fine_chans * num_corr_pols * floats_per_complex
        print(f" -> naxis1: {naxis1}")
        naxis2 = num_baselines
        print(f" -> naxis2: {naxis2}")
        floats_per_img = naxis1 * naxis2
        print(f" -> floats_per_img: {floats_per_img}")

        ####
        # Handle Scans
        ####

        primary_hdu.header['NSCANS'] = min(MAX_SCANS, primary_hdu.header['NSCANS'])

        ####
        # Analyse GPUFits
        ####

        gpufits_df = pd.DataFrame([
            dict(list(parse_filename(gpubox_name).items()) + [('name', gpubox_name)])
            for gpubox_name in TEST_GPUFITS_NAMES
        ])
        print(f" -> gpufits_df:\n{gpufits_df}")

        gpufits_coarse_chans = np.unique(gpufits_df.channel)
        print(f" -> gpufits_coarse_chans({len(gpufits_coarse_chans)}):\n{gpufits_coarse_chans}")
        valid_coarse_chans = sorted(list(
            set(gpufits_coarse_chans).intersection(set(coarse_chans))))[:MAX_COARSE_CHANS]
        num_coarse_chans = len(valid_coarse_chans)
        print(f" -> valid_coarse_chans({num_coarse_chans}):\n{valid_coarse_chans}")

        primary_hdu.header['CHANNELS'] = ','.join(f'{chan}' for chan in valid_coarse_chans)
        primary_hdu.header['CENTCHAN'] = valid_coarse_chans[num_coarse_chans//2]
        primary_hdu.header['NCHANS'] = num_coarse_chans * MAX_FINE_CHANS

        total_chan_bandwidth_hz = coarse_chan_bandwidth_hz * num_coarse_chans
        primary_hdu.header['BANDWDTH'] = total_chan_bandwidth_hz / 1_000_000

        gpufits_batches = sorted(list(np.unique(gpufits_df.batch)[:MAX_BATCHES]))
        print(f" -> gpufits_batches({len(gpufits_batches)}):\n{gpufits_batches}")

        filtered_gpufits = gpufits_df[np.isin(gpufits_df.channel, valid_coarse_chans)][np.isin(
            gpufits_df.batch, gpufits_batches)]
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

    for row in filtered_gpufits.to_dict('records'):
        gpufits_name = row['name']
        channel_index = valid_coarse_chans.index(row['channel'])
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
            primary_hdu.header['NINPUTS'] = len(filtered_inputs)
            scan_hdus = []

            scan_hdu_chunks = chunk(gpu_fits[1:1+(MAX_SCANS*2)], 2)
            for (scan_index, (img_hdu, flag_hdu)) in enumerate(scan_hdu_chunks):
                global_scan_index = channel_index
                global_scan_index = global_scan_index << math.ceil(
                    math.log2(MAX_BATCHES)) | batch_index
                global_scan_index = global_scan_index << math.ceil(
                    math.log2(MAX_SCANS)) | scan_index

                print(f" -> global_scan_index: {global_scan_index:08b}")
                # prefix scan indices to start with the char 'A'
                global_scan_index = (0x41 << 8) | global_scan_index
                print(f" -> img_hdu[{scan_index}].data.shape\n{img_hdu.data.shape}")
                print(f" -> flag_hdu[{scan_index}].data.shape\n{flag_hdu.data.shape}")
                img_hdu.header['NAXIS1'] = naxis1
                img_hdu.header['NAXIS2'] = naxis2

                # img_hdu.data = img_hdu.data[:naxis2, :naxis1]
                float_start = global_scan_index << 8 * math.ceil(math.log2(floats_per_img) / 8)
                img_hdu.data = np.arange(float_start, float_start +
                                         floats_per_img).reshape((naxis2, naxis1))
                flag_hdu.data = flag_hdu.data[:naxis2, :]

                scan_hdus.append(img_hdu)
                scan_hdus.append(flag_hdu)

            new_gpu_fits = fits.HDUList([primary_hdu] + scan_hdus)

            if not path_exists(dirname(dst_gpufits_path)):
                makedirs(dirname(dst_gpufits_path))

            print(f'-> writing to {dst_gpufits_path}')
            new_gpu_fits.writeto(dst_gpufits_path, overwrite=True)


if __name__ == '__main__':
    main()
