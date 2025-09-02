#!/usr/bin/env python3

from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs

from astropy.io import fits
from pprint import pformat
import numpy as np
import pandas as pd
import math
import re

from common import get_gpufits_num_scans, chunk

DST_TEST_DIR_MWAX = "tests/data/1297526432_mwax"
SRC_TEST_DIR_MWAX = "/mnt/data/1297526432_mwax"
TEST_METAFITS_NAME_MWAX = "1297526432.metafits"
TEST_GPUFITS_NAMES_MWAX = [
    "1297526432_20210216160014_ch117_000.fits",
    "1297526432_20210216160014_ch117_001.fits",
    "1297526432_20210216160014_ch118_000.fits",
    "1297526432_20210216160014_ch118_001.fits"
]
SRC_TEST_DIR_MWA_ORD = "/mnt/data/1196175296_vis"
TEST_METAFITS_NAME_MWA_ORD = "1196175296.metafits"
TEST_GPUFITS_NAMES_MWA_ORD = [
    "1196175296_20171201145440_gpubox01_00.fits",
    "1196175296_20171201145540_gpubox01_01.fits",
    "1196175296_20171201145440_gpubox02_00.fits",
    "1196175296_20171201145540_gpubox02_01.fits",
    # "1196175296_20171201145440_gpubox03_00.fits",
    # "1196175296_20171201145540_gpubox03_01.fits",
    # "1196175296_20171201145440_gpubox04_00.fits",
    # "1196175296_20171201145540_gpubox04_01.fits",
    # "1196175296_20171201145440_gpubox05_00.fits",
    # "1196175296_20171201145540_gpubox05_01.fits",
    # "1196175296_20171201145440_gpubox06_00.fits",
    # "1196175296_20171201145540_gpubox06_01.fits",
    # "1196175296_20171201145440_gpubox07_00.fits",
    # "1196175296_20171201145540_gpubox07_01.fits",
    # "1196175296_20171201145440_gpubox08_00.fits",
    # "1196175296_20171201145540_gpubox08_01.fits",
    # "1196175296_20171201145440_gpubox09_00.fits",
    # "1196175296_20171201145540_gpubox09_01.fits",
    # "1196175296_20171201145440_gpubox10_00.fits",
    # "1196175296_20171201145540_gpubox10_01.fits",
    # "1196175296_20171201145440_gpubox11_00.fits",
    # "1196175296_20171201145540_gpubox11_01.fits",
    # "1196175296_20171201145440_gpubox12_00.fits",
    # "1196175296_20171201145540_gpubox12_01.fits",
    # "1196175296_20171201145440_gpubox13_00.fits",
    # "1196175296_20171201145540_gpubox13_01.fits",
    # "1196175296_20171201145440_gpubox14_00.fits",
    # "1196175296_20171201145540_gpubox14_01.fits",
    # "1196175296_20171201145440_gpubox15_00.fits",
    # "1196175296_20171201145540_gpubox15_01.fits",
    # "1196175296_20171201145440_gpubox16_00.fits",
    # "1196175296_20171201145540_gpubox16_01.fits",
    # "1196175296_20171201145440_gpubox17_00.fits",
    # "1196175296_20171201145540_gpubox17_01.fits",
    # "1196175296_20171201145440_gpubox18_00.fits",
    # "1196175296_20171201145540_gpubox18_01.fits",
    # "1196175296_20171201145440_gpubox19_00.fits",
    # "1196175296_20171201145540_gpubox19_01.fits",
    # "1196175296_20171201145440_gpubox20_00.fits",
    # "1196175296_20171201145540_gpubox20_01.fits",
    # "1196175296_20171201145440_gpubox21_00.fits",
    # "1196175296_20171201145540_gpubox21_01.fits",
    # "1196175296_20171201145440_gpubox22_00.fits",
    # "1196175296_20171201145540_gpubox22_01.fits",
    # "1196175296_20171201145440_gpubox23_00.fits",
    # "1196175296_20171201145540_gpubox23_01.fits",
    # "1196175296_20171201145440_gpubox24_00.fits",
    # "1196175296_20171201145540_gpubox24_01.fits",
]
SRC_TEST_DIR_MWA_ORD_FLAGS = "/mnt/data/1247842824_vis"
TEST_METAFITS_NAME_MWA_ORD_FLAGS = "1247842824.metafits"
TEST_GPUFITS_NAMES_MWA_ORD_FLAGS = [
    "1247842824_20190722150008_gpubox01_00.fits",
]
RE_MWAX_NAME = (
    r"(?P<obsid>\d{10})_"
    r"(?P<datetime>\d{8}(.)?\d{6})_"
    r"ch(?P<rec_chan>\d{3})_"
    r"(?P<batch>\d{3}).fits"
)
RE_MWA_ORD_NAME = (
    r"(?P<obsid>\d{10})_"
    r"(?P<datetime>\d{14})_"
    r"gpubox(?P<gpubox_num>\d{2})_"
    r"(?P<batch>\d{2}).fits"
)

MAX_COARSE_CHANS = 2
MAX_BATCHES = 2
MAX_SCANS = 2
MAX_ANTENNAS = 2
MAX_FINE_CHANS = 2


def parse_filename(
    name,
    corr_type="MWAX",
    metafits_coarse_chans=None,
    with_src_dir=(lambda x: x),
):
    if metafits_coarse_chans is None:
        metafits_coarse_chans = []
    result = {}
    if corr_type == "MWAX":
        result = re.match(RE_MWAX_NAME, name).groupdict()
        result['corr_chan'] = metafits_coarse_chans.index(result['rec_chan'])
    elif corr_type == "MWA_ORD":
        result = re.match(RE_MWA_ORD_NAME, name).groupdict()
        result['corr_chan'] = int(result['gpubox_num']) - 1
        result['rec_chan'] = metafits_coarse_chans[result['corr_chan']]
    result['name'] = name
    with fits.open(with_src_dir(name)) as gpu_fits:
        num_hdus = len(gpu_fits)
        result['hdus'] = num_hdus
        result['time'] = (
            gpu_fits[1].header['TIME'] +
            gpu_fits[1].header['MILLITIM'] / 1000
        )

        result['nscans'] = get_gpufits_num_scans(num_hdus, corr_type)
    return result


def get_input_df(tile_data):
    input_cols = ['Input', 'Antenna', 'Tile', 'TileName', 'Pol', 'Rx', 'Slot']
    return pd.DataFrame(dict([
        (col, list(tile_data.field(col)[:]))
        for col in input_cols
    ]))


def get_limited_antennas(input_df, max_antennas=None):
    antennas = list(np.unique(input_df.Antenna))
    if max_antennas:
        antennas = antennas[:max_antennas]
    return antennas


def get_global_scan_index(
    channel_index,
    batch_index,
    max_batches,
    scan_index,
    max_scans,
):
    result = channel_index
    result = result << math.ceil(math.log2(max_batches)) | batch_index
    result = result << math.ceil(math.log2(max_scans)) | scan_index
    return result


def display_float(flt):
    return (
        f"{flt} (0x{flt:08x}, {flt:064b}, ~2^{math.log2(flt + 1)})"
    )


def split_strip_filter(str):
    return list(
        filter(None, map(lambda tok: tok.strip(), str.split(',')))
    )


"""
Notes and references (wrapped for linting):

for mwalib/mwax:
    - fine chan width = FINECHAN
    - coarse chan width = bandwidth / len(CHANNELS)

for mwalib/legacy:
    - coarse_chan_width_hz = BANDWIDTH / CHANNELS.len()
    - num_corr_fine_chans_per_coarse =
      metafits_coarse_chan_width_hz / FINECHAN
    - naxis2 = metafits_fine_chans_per_coarse

for cotter:

    double ChannelFrequencyHz(size_t channelIndex) const {
        const double channelWidthMHz =
            _header.bandwidthMHz / _header.nChannels;
        return (
            _header.centralFrequencyMHz +
            (channelIndex - _header.nChannels*0.5) * channelWidthMHz
        ) * 1000000.0;
    }

    double ChannelFrequencyHz(size_t coarseChannelNr,
                               size_t channelIndexInSubband) const {
        const double channelWidthHz =
            1000000.0*_header.bandwidthMHz / _header.nChannels;
        return (
            (double(coarseChannelNr) - 0.5) * 1280000.0 +
            double(channelIndexInSubband) * channelWidthHz
        );
    }

    std::cout << "Observation covers "
              << (ChannelFrequencyHz(0)/1000000.0)
              << '-'
              << (ChannelFrequencyHz(_header.nChannels-1)/1000000.0)
              << " MHz.\n";

    timeAvgFactor = round(timeRes_s / INTTIME);
    if (timeAvgFactor == 0)
        timeAvgFactor = 1;
    timeRes_s = timeAvgFactor * INTTIME;
    freqAvgFactor = round(freqRes_kHz / (1000.0 * BANDWIDTH / N_CHANS));
    if (freqAvgFactor == 0)
        freqAvgFactor = 1;
    freqRes_kHz = freqAvgFactor * (1000.0 * BANDWIDTH / N_CHANS);

    _subbandCount = `-sbcount`

    nChannelsPerSubband = N_CHANS / _subbandCount

"""


def generate(args):
    print(f"generating to {args['dst_dir']}")

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
        # print(
        #     f" -> meta_fits[1].data.shape\n{repr(meta_fits[1].data.shape)}"
        # )

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

        primary_hdu.header['NINPUTS'] = num_inputs
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
            split_strip_filter(primary_hdu.header['RECVRS']),
            split_strip_filter(primary_hdu.header['DELAYS'])
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

        metafits_coarse_chans = split_strip_filter(
            primary_hdu.header['CHANNELS']
        )
        num_metafits_coarse_chans = len(metafits_coarse_chans)
        print(
            f" -> metafits_coarse_chans({num_metafits_coarse_chans}):\n"
            f"{metafits_coarse_chans}"
        )

        ####
        # Analyse GPUFits
        ####

        print("anaylsing gpufits...")

        gpufits_df = pd.DataFrame([
            parse_filename(
                gpubox_name,
                args['corr_type'],
                metafits_coarse_chans,
                with_src_dir,
            )
            for gpubox_name in args['gpufits_names']
        ])
        print(f" -> gpufits_df:\n{gpufits_df}")

        rec_coarse_chans = np.unique(gpufits_df.rec_chan)
        print(
            f" -> rec_coarse_chans({len(rec_coarse_chans)}):\n"
            f"{rec_coarse_chans}"
        )
        # valid_coarse_chans = sorted(rec_coarse_chans)
        valid_coarse_chans = list(rec_coarse_chans)
        if 'max_coarse_chans' in args:
            valid_coarse_chans = valid_coarse_chans[:args['max_coarse_chans']]
        num_coarse_chans = len(valid_coarse_chans)
        print(
            f" -> valid_coarse_chans({num_coarse_chans}):\n"
            f"{valid_coarse_chans}"
        )

        # batches_df = sorted(
        #     list(np.unique(gpufits_df.batch)[:args['max_batches']])
        # )
        batches_df = pd.DataFrame(
            sorted(list(
                (
                    batch,
                    np.max(gpufits_df.nscans[gpufits_df.batch == batch])
                )
                for batch in np.unique(gpufits_df.batch)
            )),
            columns=['batch', 'max_scans']
        )
        num_batches = len(batches_df)
        max_scans_per_batch = np.max(batches_df.max_scans)
        print(
            f" -> batches_df (len: {num_batches},"
            f" max: {max_scans_per_batch}):\n"
            f"{batches_df}"
        )
        if args.get('max_batches'):
            batches_df = batches_df[:args['max_batches']]
        if args.get('max_scans'):
            batches_df.max_scans = batches_df.max_scans.apply(
                lambda x: min(x, args['max_scans'])
            )
        num_batches = len(batches_df)
        max_scans_per_batch = np.max(batches_df.max_scans)
        print(
            f" -> batches_df (limited, len: {num_batches},"
            f" max: {max_scans_per_batch}):\n{batches_df}"
        )

        filtered_gpufits = (
            gpufits_df[
                np.isin(gpufits_df.rec_chan, valid_coarse_chans) &
                np.isin(gpufits_df.batch, batches_df.batch)
            ]
            .sort_values(['corr_chan', 'batch'])
        )
        print(f" -> filtered_gpufits:\n{filtered_gpufits}")

        ####
        # Handle Scans
        ####

        num_scans_total = primary_hdu.header['NSCANS']
        print(f" -> num_scans_total:\n{num_scans_total}")
        num_scans_total = np.sum(batches_df.max_scans)
        print(f" -> num_scans_total (limited):\n{num_scans_total}")
        primary_hdu.header['NSCANS'] = num_scans_total
        int_time = primary_hdu.header['INTTIME']
        if int_time > 0:
            primary_hdu.header['EXPOSURE'] = int(num_scans_total * int_time)

        ####
        # Handle Coarse/Fine Channel Bandwidth
        ####

        # Header in KHz
        fine_chan_bandwidth_hz = int(primary_hdu.header['FINECHAN'] * 1_000)
        print(f" -> fine_chan_bandwidth (Hz): {fine_chan_bandwidth_hz}")
        # Header in MHz
        total_chan_bandwidth_hz = int(
            primary_hdu.header['BANDWDTH'] * 1_000_000
        )
        print(
            f" -> total_chan_bandwidth (Hz): {total_chan_bandwidth_hz}"
        )
        coarse_chan_bandwidth_hz = (
            total_chan_bandwidth_hz // num_metafits_coarse_chans
        )
        # coarse_chan_bandwidth_hz = (
        #     total_chan_bandwidth_hz // num_coarse_chans
        # )
        print(
            f" -> coarse_chan_bandwidth (Hz, derived): "
            f"{coarse_chan_bandwidth_hz}"
        )
        num_fine_chans = coarse_chan_bandwidth_hz // fine_chan_bandwidth_hz
        print(f" -> num_fine_chans (derived): {num_fine_chans}")
        if args.get('max_fine_chans'):
            num_fine_chans = min(args['max_fine_chans'], num_fine_chans)
            print(f" -> num_fine_chans (limited): {num_fine_chans}")
            fine_chan_bandwidth_hz = coarse_chan_bandwidth_hz // num_fine_chans
            print(
                f" -> fine_chan_bandwidth (Hz, limited): "
                f"{fine_chan_bandwidth_hz}"
            )
            coarse_chan_bandwidth_hz = fine_chan_bandwidth_hz * num_fine_chans
            print(
                f" -> coarse_chan_bandwidth (Hz, limited): "
                f"{coarse_chan_bandwidth_hz}"
            )
        # total_chan_bandwidth_hz = coarse_chan_bandwidth_hz * num_coarse_chans
        total_chan_bandwidth_hz = (
            coarse_chan_bandwidth_hz * num_metafits_coarse_chans
        )
        print(
            f" -> total_chan_bandwidth (Hz, limited): "
            f"{total_chan_bandwidth_hz}"
        )
        assert (
            total_chan_bandwidth_hz / num_metafits_coarse_chans /
            fine_chan_bandwidth_hz ==
            num_fine_chans
        )
        primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000
        primary_hdu.header['BANDWDTH'] = total_chan_bandwidth_hz / 1_000_000

        ####
        # Handle Fine Channel Selection
        ####

        # coarse_channel_selection = split_strip_filter(
        #     primary_hdu.header['CHANSEL']
        # )
        coarse_channel_selection = [
            # str(metafits_coarse_chans.index(i)) for i in rec_coarse_chans
            str(metafits_coarse_chans.index(i)) for i in valid_coarse_chans
        ]
        # print(
        #     f" -> coarse_channel_selection (derived):"
        #     f" {coarse_channel_selection}"
        # )
        # coarse_channel_selection = (
        #     coarse_channel_selection[:num_coarse_chans]
        # )
        print(
            f" -> coarse_channel_selection: {coarse_channel_selection}"
        )
        # primary_hdu.header['CENTCHAN'] = (
        #     coarse_channel_selection[num_coarse_chans//2]
        # )
        # primary_hdu.header['CHANNELS'] = ','.join(
        #     f'{chan}' for chan in valid_coarse_chans
        # )
        primary_hdu.header['CHANSEL'] = ','.join(coarse_channel_selection)

        # Observation covers 200.32-231.03 MHz.
        # Output resolution: 0.5 s / 10 kHz (time avg: 1x, freq avg: 1x).
        # Using a-priori subband passband with 1536 channels.

        print(
            f" -> num_metafits_coarse_chans: {num_metafits_coarse_chans} "
            f"num_fine_chans: {num_fine_chans}"
        )
        primary_hdu.header['NCHANS'] = (
            num_metafits_coarse_chans * num_fine_chans
        )
        # primary_hdu.header['CHANNELS'] = ','.join(
        #     f'{chan}' for chan in valid_coarse_chans
        # )

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
        # Show modifications
        ####

        print(f" -> primary_hdu.header\n{repr(primary_hdu.header)}")
        meta_fits[1].data

        ####
        # write MetaFits
        ####

        if not path_exists(dirname(dst_metafits_path)):
            makedirs(dirname(dst_metafits_path))

        meta_fits.writeto(dst_metafits_path, overwrite=True)

        # now write cotter-friendly
        if args['corr_type'] == "MWA_ORD":
            print(
                f" -> num_coarse_chans: {num_coarse_chans} "
                f"num_fine_chans: {num_fine_chans}"
            )
            total_chan_bandwidth_hz = (
                coarse_chan_bandwidth_hz * num_coarse_chans
            )
            print(
                f" -> total_chan_bandwidth (Hz, limited): "
                f"{total_chan_bandwidth_hz}"
            )
            assert (
                total_chan_bandwidth_hz / num_coarse_chans /
                fine_chan_bandwidth_hz ==
                num_fine_chans
            )
            primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000
            primary_hdu.header['BANDWDTH'] = (
                total_chan_bandwidth_hz / 1_000_000
            )
            num_chans = num_coarse_chans * num_fine_chans
            primary_hdu.header['NCHANS'] = num_chans

            assert (
                total_chan_bandwidth_hz / num_chans ==
                fine_chan_bandwidth_hz
            )

            meta_fits.writeto(dst_metafits_path.replace(
                ".metafits", ".cotter.metafits"), overwrite=True)

    ####
    # Handle GPUFits
    ####

    float_count = 0
    start_times = {}

    for file_idx, row in enumerate(filtered_gpufits.to_dict('records')):
        gpufits_name = row['name']
        print(f" -> gpufits_name: {gpufits_name}")
        # channel_index = row['corr_chan']
        channel_index = valid_coarse_chans.index(row['rec_chan'])
        print(f" -> channel_index: 0x{channel_index:08b}")
        batch_mask = (batches_df.batch == row['batch'])
        batch_index = batches_df.batch.index[batch_mask][0]
        gpufits_path = with_src_dir(gpufits_name)
        dst_gpufits_path = with_dst_dir(gpufits_name)
        with fits.open(gpufits_path) as gpu_fits:
            # print(
            #     f" -> gpu_fits[0].info()\n{pformat(gpu_fits.info())}"
            # )
            print(
                f" -> gpu_fits[0].header\n{repr(gpu_fits[0].header)}"
            )
            print(
                f" -> gpu_fits[1].header\n{repr(gpu_fits[1].header)}"
            )
            # print(f" -> gpu_fits[2].header\n{repr(gpu_fits[2].header)}")

            time = (
                gpu_fits[1].header['TIME'] +
                gpu_fits[1].header['MILLITIM'] / 1000
            )
            if batch_index == 0:
                start_times[channel_index] = time
            elif args.get('max_scans'):
                cum_scans = np.max(
                    batches_df.max_scans[batches_df.index < batch_index]
                )
                time = start_times[channel_index] + (int_time * cum_scans)
                gpu_fits[0].header['TIME'] = int(time)
                gpu_fits[0].header['MILLITIM'] = int(1000 * (time % 1))
                gpu_fits[1].header['TIME'] = int(time)
                gpu_fits[1].header['MILLITIM'] = int(1000 * (time % 1))

            primary_hdu = gpu_fits[0]
            primary_hdu.header['NFINECHS'] = num_fine_chans
            primary_hdu.header['FINECHAN'] = fine_chan_bandwidth_hz / 1_000
            primary_hdu.header['NINPUTS'] = num_inputs
            scan_hdus = []

            scan_mask = (batches_df.batch == row['batch'])
            num_scans = list(batches_df.max_scans[scan_mask])[0]

            if args['corr_type'] == "MWAX":
                primary_hdu.header['CORR_VER'] = 2

                scan_hdu_chunks = chunk(gpu_fits[1:][:num_scans * 2], 2)
                for (
                    scan_index,
                    (img_hdu, flag_hdu),
                ) in enumerate(scan_hdu_chunks):
                    # print(f" -> img_hdu.header\n{repr(img_hdu.header)}")
                    # print(f" -> flag_hdu.header\n{repr(flag_hdu.header)}")
                    hdu_time = time + int_time * scan_index

                    print(
                        f" -> img_hdu[{scan_index}].data ("
                        f"{img_hdu.data.shape}, {img_hdu.data.dtype}): \n"
                        f"{img_hdu.data}"
                    )
                    img_hdu.header['TIME'] = int(hdu_time)
                    img_hdu.header['MILLITIM'] = int(1000 * (hdu_time % 1))

                    if args.get('rewrite_viz'):
                        global_scan_index = get_global_scan_index(
                            channel_index,
                            batch_index,
                            num_batches,
                            scan_index,
                            num_scans,
                        )
                        # prefix scan indices to start with the char 'A'
                        global_scan_index = (0x41 << 8) | global_scan_index
                        print(
                            f" -> global_scan_index: {global_scan_index:08b}"
                        )
                        img_hdu.header['NAXIS1'] = naxis1
                        img_hdu.header['NAXIS2'] = naxis2
                        float_start = (
                            global_scan_index <<
                            8 * math.ceil(math.log2(floats_per_img) / 8)
                        )
                        print(
                            f" -> float_start: {display_float(float_start)}"
                        )
                        float_end = float_start + floats_per_img
                        print(
                            f" -> float_end: {display_float(float_end)}"
                        )
                        arr = np.arange(float_start, float_end)
                        img_hdu.data = arr.reshape(
                            (naxis2, naxis1)
                        ).astype(np.int32)

                    print(
                        f" -> flag_hdu[{scan_index}].data.shape("
                        f"{flag_hdu.data.shape}, {flag_hdu.data.dtype})"
                    )
                    flag_hdu.header['TIME'] = int(hdu_time)
                    flag_hdu.header['MILLITIM'] = int(1000 * (hdu_time % 1))
                    if args.get('rewrite_viz'):
                        flag_hdu.data = flag_hdu.data[:naxis2, :]
                    scan_hdus.append(img_hdu)
                    scan_hdus.append(flag_hdu)

            elif args['corr_type'] == "MWA_ORD":
                for (
                    scan_index,
                    img_hdu,
                ) in enumerate(gpu_fits[1:][:num_scans]):
                    # print(f" -> img_hdu.header\n{repr(img_hdu.header)}")
                    hdu_time = time + int_time * scan_index

                    print(
                        f" -> img_hdu[{scan_index}].data ("
                        f"{img_hdu.data.shape}, {img_hdu.data.dtype}): \n"
                        f"{img_hdu.data}"
                    )
                    img_hdu.header['TIME'] = int(hdu_time)
                    img_hdu.header['MILLITIM'] = int(1000 * (hdu_time % 1))

                    if args.get('rewrite_viz'):
                        img_hdu.header['NAXIS1'] = naxis1
                        img_hdu.header['NAXIS2'] = naxis2
                        global_scan_index = get_global_scan_index(
                            channel_index,
                            batch_index,
                            num_batches,
                            scan_index,
                            num_scans,
                        )
                        # global_scan_index = (0x1001 << 8) | global_scan_index
                        print(
                            f" -> global_scan_index: {global_scan_index:08b}"
                        )
                        if args.get('sequential'):
                            float_start = float_count
                        else:
                            float_start = (
                                global_scan_index <<
                                math.ceil(math.log2(floats_per_img))
                            )
                        print(f" -> float_start: {display_float(float_start)}")
                        float_end = float_start + floats_per_img
                        print(f" -> float_end: {display_float(float_end)}")
                        arr = np.arange(float_start, float_end)
                        img_hdu.data = arr.reshape(
                            (naxis2, naxis1)
                        ).astype(np.float64)

                    float_count += floats_per_img

                    # print(
                    #     f" -> modified img_hdu.header\n"
                    #     f"{repr(img_hdu.header)}"
                    # )
                    scan_hdus.append(img_hdu)

            new_gpu_fits = fits.HDUList([primary_hdu] + scan_hdus)
            # print(
            #     f" -> new_gpu_fits[0].header\n{repr(new_gpu_fits[0].header)}"
            # )
            # print(
            #     f" -> new_gpu_fits[1].header\n{repr(new_gpu_fits[1].header)}"
            # )

            print(
                f"-> writing {len(scan_hdus)} scans to {dst_gpufits_path}"
            )
            if not path_exists(dirname(dst_gpufits_path)):
                makedirs(dirname(dst_gpufits_path))
            new_gpu_fits.writeto(dst_gpufits_path, overwrite=True)


def main():
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
        'rewrite_viz': True,
    })
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
        #   limiting ord antennas breaks mwalib and cotter.
        'max_fine_chans': MAX_FINE_CHANS,
        'rewrite_viz': True,
    })
    generate({
        'corr_type': "MWA_ORD",
        'src_dir': SRC_TEST_DIR_MWA_ORD_FLAGS,
        'dst_dir': "tests/data/1247842824_flags",
        'metafits_name': TEST_METAFITS_NAME_MWA_ORD_FLAGS,
        'gpufits_names': TEST_GPUFITS_NAMES_MWA_ORD_FLAGS,
        'max_coarse_chans': 1,
        'max_batches': 1,
        'max_scans': 2,
        # 'max_antennas': MAX_ANTENNAS,
        #   limiting ord antennas breaks mwalib and cotter.
        # 'max_fine_chans': MAX_FINE_CHANS,
        'rewrite_viz': False,
    })


if __name__ == '__main__':
    main()
