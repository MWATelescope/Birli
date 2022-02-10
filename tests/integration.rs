use birli::cli::main_with_args;
use float_cmp::F32Margin;
use std::path::PathBuf;
use tempfile::tempdir;

mod common;
use common::{compare_ms_with_csv, compare_uvfits_with_csv, get_1254670392_avg_paths};

#[test]
fn test_1254670392_avg_uvfits_none_chunked() {
    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1254670392.none.chunked.uvfits");
    let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

    let expected_csv_path =
        PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

    let mut args = vec![
        "birli",
        "-m",
        metafits_path,
        "-u",
        uvfits_path.to_str().unwrap(),
        "--no-digital-gains",
        "--no-draw-progress",
        "--no-cable-delay",
        "--no-geometric-delay",
        "--no-rfi",
        "--emulate-cotter",
        "--time-chunk",
        "1",
    ];
    args.extend_from_slice(&gpufits_paths);

    main_with_args(&args);

    compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), true);
}

/// Data generated with
///
/// ```bash
/// cotter \
///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms \
///  -allowmissing \
///  -edgewidth 0 \
///  -endflag 0 \
///  -initflag 0 \
///  -noantennapruning \
///  -nocablelength \
///  -norfi \
///  -nogeom \
///  -noflagautos \
///  -noflagdcchannels \
///  -nosbgains \
///  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
///  -nostats \
///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
///  -full-apply tests/data/1254670392_avg/1254690096.bin \
///  tests/data/1254670392_avg/*gpubox*.fits
/// ```
///
/// then casa
///
/// ```bash
/// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms')
/// exec(open('tests/data/casa_dump_ms.py').read())
/// ```
#[test]
fn test_1254670392_avg_ms_none_norfi_cal() {
    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.norfi.cal.ms");
    let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

    let expected_csv_path =
        PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms.csv");

    let mut args = vec![
        "birli",
        "-m",
        metafits_path,
        "-M",
        ms_path.to_str().unwrap(),
        "--no-digital-gains",
        "--no-draw-progress",
        "--no-cable-delay",
        "--no-geometric-delay",
        "--no-rfi",
        "--emulate-cotter",
        "--apply-di-cal",
        "tests/data/1254670392_avg/1254690096.bin",
    ];
    args.extend_from_slice(&gpufits_paths);

    main_with_args(&args);

    // ignoring weights because Cotter doesn't flag NaNs
    compare_ms_with_csv(ms_path, expected_csv_path, F32Margin::default(), true);
}

/// Data generated with
///
/// ```bash
/// cotter \
///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms \
///  -allowmissing \
///  -edgewidth 0 \
///  -endflag 0 \
///  -initflag 0 \
///  -noantennapruning \
///  -nocablelength \
///  -norfi \
///  -nogeom \
///  -noflagautos \
///  -noflagdcchannels \
///  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
///  -nostats \
///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
///  tests/data/1254670392_avg/*gpubox*.fits
/// ```
///
/// then casa
///
/// ```bash
/// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms')
/// exec(open('tests/data/casa_dump_ms.py').read())
/// ```
#[test]
fn test_1254670392_avg_ms_none_norfi_nodigital() {
    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.norfi.nodigital.ms");
    let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

    let expected_csv_path =
        PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms.csv");

    let mut args = vec![
        "birli",
        "-m",
        metafits_path,
        "-M",
        ms_path.to_str().unwrap(),
        "--no-draw-progress",
        "--no-cable-delay",
        "--no-geometric-delay",
        "--no-rfi",
        "--emulate-cotter",
    ];
    args.extend_from_slice(&gpufits_paths);

    main_with_args(&args);
    compare_ms_with_csv(
        ms_path,
        expected_csv_path,
        F32Margin::default().epsilon(7e-5),
        false,
    );
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use super::*;
    use birli::{CorrelatorContext, FlagFileSet};
    use marlu::rubbl_casatables::{Table, TableOpenMode};
    use itertools::izip;

    macro_rules! assert_flagsets_eq {
        ($context:expr, $left_flagset:expr, $right_flagset:expr, $gpubox_ids:expr) => {
            let num_baselines = $context.metafits_context.num_baselines;
            let num_flags_per_row = $context.metafits_context.num_corr_fine_chans_per_coarse;
            let num_common_timesteps = $context.num_common_timesteps;
            let num_rows = num_common_timesteps * num_baselines;
            let num_flags_per_timestep = num_baselines * num_flags_per_row;

            assert!(num_baselines > 0);
            assert!(num_rows > 0);
            assert!(num_flags_per_row > 0);

            let right_chan_header_flags_raw =
            $right_flagset.read_chan_header_flags_raw().unwrap();

            let left_chan_header_flags_raw = $left_flagset.read_chan_header_flags_raw().unwrap();

            for gpubox_id in $gpubox_ids {
                let (left_header, left_flags) = left_chan_header_flags_raw.get(&gpubox_id).unwrap();
                let (right_header, right_flags) =
                    right_chan_header_flags_raw.get(&gpubox_id).unwrap();
                assert_eq!(left_header.obs_id, right_header.obs_id);
                assert_eq!(left_header.num_channels, right_header.num_channels);
                assert_eq!(left_header.num_ants, right_header.num_ants);
                // assert_eq!(left_header.num_common_timesteps, right_header.num_common_timesteps);
                assert_eq!(left_header.num_timesteps, num_common_timesteps);
                assert_eq!(left_header.num_pols, right_header.num_pols);
                assert_eq!(left_header.gpubox_id, right_header.gpubox_id);
                assert_eq!(left_header.bytes_per_row, right_header.bytes_per_row);
                // assert_eq!(left_header.num_rows, right_header.num_rows);
                assert_eq!(left_header.num_rows, num_rows);

                // assert_eq!(left_flags.len(), right_flags.len());
                assert_eq!(
                    left_flags.len(),
                    num_common_timesteps * num_baselines * num_flags_per_row
                );

                izip!(
                    left_flags.chunks(num_flags_per_timestep),
                    right_flags.chunks(num_flags_per_timestep)
                ).enumerate().for_each(|(common_timestep_idx, (left_timestep_chunk, right_timestep_chunk))| {
                    izip!(
                        $context.metafits_context.baselines.iter(),
                        left_timestep_chunk.chunks(num_flags_per_row),
                        right_timestep_chunk.chunks(num_flags_per_row)
                    ).enumerate().for_each(|(baseline_idx, (baseline, left_baseline_chunk, right_baseline_chunk))| {
                        if baseline.ant1_index == baseline.ant2_index {
                            return
                        }

                        assert_eq!(
                            left_baseline_chunk, right_baseline_chunk,
                            "flag chunks for common timestep {}, baseline {} (ants {}, {}) do not match! \nbirli:\n{:?}\ncotter:\n{:?}",
                            common_timestep_idx, baseline_idx, baseline.ant1_index, baseline.ant2_index, left_baseline_chunk, right_baseline_chunk
                        )
                    });
                });
            }
        };
    }

    pub(crate) use assert_flagsets_eq;

    #[test]
    fn aoflagger_outputs_flags() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--no-draw-progress",
            "-f",
            mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let mut birli_flag_file_set = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        let mut cotter_flag_file_set = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        assert_flagsets_eq!(
            context,
            birli_flag_file_set,
            cotter_flag_file_set,
            gpubox_ids
        );
    }

    #[test]
    #[ignore = "chunks not supported for .mwaf"]
    fn aoflagger_outputs_flags_chunked() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--no-draw-progress",
            "--time-chunk", "1",
            "-f",
            mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let mut birli_flag_file_set = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        let mut cotter_flag_file_set = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        assert_flagsets_eq!(
            context,
            birli_flag_file_set,
            cotter_flag_file_set,
            gpubox_ids
        );
    }

    #[test]
    fn aoflagger_outputs_uvfits() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1247842824.uvfits");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        assert!(uvfits_path.exists());

        assert!(uvfits_path.metadata().unwrap().len() > 0);
    }

    #[test]
    fn test_1254670392_avg_uvfits_no_corrections() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), false);
    }

    #[test]
    #[ignore = "slow"]
    fn test_1254670392_avg_uvfits_cable_only() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.cable.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.cable.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
            false,
        );
    }

    #[test]
    #[ignore = "slow"]
    fn test_1254670392_avg_uvfits_geom_only() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.geom.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--no-cable-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.geom.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
            false,
        );
    }

    #[test]
    fn test_1254670392_avg_uvfits_both() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path =
        //     PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.corrected.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
        );
    }

    /// Test corrections using arbitrary phase centre.
    /// data generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -centre 00h00m00.0s 00d00m00.0s \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// ```bash
    /// python tests/data/dump_uvfits.py \
    ///     tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits \
    ///     --timestep-limit=12 --baseline-limit=12 --dump-mode=vis-only \
    ///     --dump-csv=tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits.csv
    /// ```
    #[test]
    fn test_1254670392_avg_uvfits_both_phase0() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits.csv",
        );

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--phase-centre",
            "0.0",
            "0.0",
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
        );
    }

    /// Test corrections using pointing centre as phase centre.
    /// data generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -usepcentre \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// Then the following CASA commands
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn test_1254670392_avg_ms_phase_pointing() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms.csv",
        );

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--pointing-centre",
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
            false,
        );
    }

    /// Test generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// then the following casa commands:
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.corrected.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn test_1254670392_avg_ms_corrected() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.ms.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
        );
    }

    /// Test generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -nocablelength \
    ///   -nogeom \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///   -timeres 4 \
    ///   -freqres 160 \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// then the following casa commands:
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn test_1254670392_avg_ms_none_avg_4s_160khz() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-res",
            "4",
            "--avg-freq-res",
            "160",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
        );
    }

    /// Same as above but with factors instead of resolution
    #[test]
    fn test_1254670392_avg_ms_none_avg_4s_160khz_factors() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor",
            "2",
            "--avg-freq-factor",
            "4",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
        );
    }

    /// Same as above but forcing chunks by using a small --max-memory
    #[test]
    fn test_1254670392_avg_ms_none_avg_4s_160khz_max_mem() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--sel-time", "0", "2",
            "--time-chunk", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        let main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        assert_eq!(main_table.n_rows(), 2 * 8256);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
        );
    }

    /// test when time_chunk is not a multiple of avg_time
    #[test]
    #[should_panic]
    fn test_1254670392_avg_ms_none_avg_4s_160khz_tiny_chunk() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--sel-time", "0", "2",
            "--time-chunk", "1",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
    }
}
