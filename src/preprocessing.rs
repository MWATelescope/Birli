//! Crate for preprocessing visibilities
use crate::{
    calibration::apply_di_calsol,
    correct_cable_lengths, correct_geometry,
    corrections::{correct_coarse_passband_gains, correct_digital_gains, ScrunchType},
    marlu::{mwalib::CorrelatorContext, ndarray::prelude::*, Jones, LatLngHeight, RADec},
    van_vleck::{correct_van_vleck, get_vv_sample_scale},
    with_increment_duration, BirliError, VisSelection,
};
use cfg_if::cfg_if;
use derive_builder::Builder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::trace;
use std::{
    fmt::{Debug, Display},
    time::Duration,
};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use crate::{
            flags::flag_jones_array_existing,
        };
        use aoflagger_sys::{cxx_aoflagger_new};
    }
}

/// Options for preprocessing a chunk of correlator data
#[derive(Builder, Debug, Default, Clone)]
pub struct PreprocessContext<'a> {
    /// The array position used for geometric corrections
    pub array_pos: LatLngHeight,
    /// The phase centre used for geometric corrections
    pub phase_centre: RADec,

    /// Whether Van Vleck corrections are enabled
    #[builder(default = "false")]
    pub correct_van_vleck: bool,
    /// Whether cable length corrections are enabled
    #[builder(default = "true")]
    pub correct_cable_lengths: bool,
    /// Whether digital gain corrections are enabled
    #[builder(default = "true")]
    pub correct_digital_gains: bool,
    /// the pfb passband gains to use for corrections
    pub passband_gains: Option<&'a [f64]>,
    /// The calibration solutions to apply
    pub calsols: Option<Array2<Jones<f64>>>,
    /// Whether geometric corrections are enabled
    #[builder(default = "true")]
    pub correct_geometry: bool,

    /// `AOFlagger` strategy path for flagging
    #[builder(default)]
    #[cfg(feature = "aoflagger")]
    pub aoflagger_strategy: Option<String>,

    /// Whether to draw progress bars
    #[builder(default = "true")]
    pub draw_progress: bool,
}

impl Display for PreprocessContext<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{} correct Van Vleck.",
            if self.correct_van_vleck {
                "Will"
            } else {
                "Will not"
            }
        )?;
        writeln!(
            f,
            "{} correct cable lengths.",
            if self.correct_cable_lengths {
                "Will"
            } else {
                "Will not"
            }
        )?;
        writeln!(
            f,
            "{} correct digital gains.",
            if self.correct_digital_gains {
                "Will"
            } else {
                "Will not"
            }
        )?;
        writeln!(
            f,
            "{} correct coarse pfb passband gains.",
            if self.passband_gains.is_some() {
                "Will"
            } else {
                "Will not"
            }
        )?;
        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                if let Some(strategy) = &self.aoflagger_strategy {
                    writeln!(f, "Will flag with aoflagger strategy {strategy}")?;
                } else {
                    writeln!(f, "Will not flag with aoflagger")?;
                }
            }
        }
        writeln!(
            f,
            "{} correct geometry.",
            if self.correct_geometry {
                "Will"
            } else {
                "Will not"
            }
        )?;
        Ok(())
    }
}

impl PreprocessContext<'_> {
    /// A one line description of the tasks preprocessing will do.
    pub fn as_comment(&self) -> String {
        [
            if self.correct_van_vleck {
                Some("Van Vleck corrections".to_string())
            } else {
                None
            },
            if self.correct_cable_lengths {
                Some("cable length corrections".to_string())
            } else {
                None
            },
            if self.correct_digital_gains {
                Some("digital gains".to_string())
            } else {
                None
            },
            if self.passband_gains.is_some() {
                Some("pfb gains".to_string())
            } else {
                None
            },
            #[cfg(feature = "aoflagger")]
            self.aoflagger_strategy
                .as_ref()
                .map(|strategy| format!("aoflagging with {strategy}")),
            if self.correct_geometry {
                Some("geometric corrections".to_string())
            } else {
                None
            },
        ]
        .into_iter()
        .flatten()
        .collect::<Vec<String>>()
        .join(", ")
    }

    /// Preprocess visibilities for a chunk of correlator data
    ///
    /// # Arguments
    /// * `corr_ctx` - [`marlu::mwalib::CorrelatorContext`]
    /// * `jones_array` - Array of Jones visibilties
    /// * `weight_array` - Array of weights associated with Jones visibilities
    /// * `flag_array` - Array of flags associated with Jones visibilities
    ///
    /// # Errors
    /// will wrap errors from `correct_van_vleck`, `correct_cable_lengths`, `correct_digital_gains`,
    /// `correct_coarse_passband_gains`, `flag_jones_array_existing`, `correct_geometry`,
    /// `apply_di_calsol`
    #[allow(clippy::too_many_arguments)]
    pub fn preprocess(
        &self,
        corr_ctx: &CorrelatorContext,
        mut jones_array: ArrayViewMut3<Jones<f32>>,
        mut weight_array: ArrayViewMut3<f32>,
        mut flag_array: ArrayViewMut3<bool>,
        vis_sel: &VisSelection,
    ) -> Result<(), BirliError> {
        let sel_ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

        if self.correct_van_vleck
            || self.correct_cable_lengths
            || self.correct_digital_gains
            || self.passband_gains.is_some()
        {
            // get selected antenna pairs, and flagged antenna
            let flagged_ants = corr_ctx
                .metafits_context
                .antennas
                .iter()
                .enumerate()
                .filter_map(|(idx, ant)| {
                    if ant.rfinput_x.flagged || ant.rfinput_y.flagged {
                        Some(idx)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            let sample_scale = get_vv_sample_scale(corr_ctx)?;
            // mwalib reads data one timestep and coarse channel at a time.
            // each of these reads can fail, in which case the whole chunk is flagged.

            // For each timestep and coarse channel pair
            let (num_timesteps, num_channels, _) = jones_array.dim();
            let num_coarse_chans = num_channels / fine_chans_per_coarse;
            let scrunch_type =
                ScrunchType::from_mwa_version(corr_ctx.metafits_context.mwa_version.unwrap())?;
            // combined progress bar for vv, cable, digital
            let draw_target = if self.draw_progress {
                ProgressDrawTarget::stderr()
            } else {
                ProgressDrawTarget::hidden()
            };

            // Create a progress bar to show the status of the correction
            let preflag_progress = ProgressBar::with_draw_target(
                Some(num_coarse_chans as u64 * num_timesteps as u64),
                draw_target,
            );
            preflag_progress.set_style(
                ProgressStyle::default_bar()
                    .template(
                        "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                    )
                    .unwrap()
                    .progress_chars("=> "),
            );
            preflag_progress.set_message("preflag");

            // iterate over coarse channels
            for (mut jones_chunk, mut weight_chunk, flag_chunk, coarse_chan_idx) in izip!(
                jones_array.axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse),
                weight_array.axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse),
                flag_array.axis_chunks_iter(Axis(1), fine_chans_per_coarse),
                vis_sel.coarse_chan_range.clone(), // todo: .chunks(size)?,
            ) {
                let coarse_chan_range = coarse_chan_idx..(coarse_chan_idx + 1);

                // determine which timestep ranges are not completely flagged
                let unflagged_timestep_ranges = flag_chunk
                    .axis_iter(Axis(0))
                    .enumerate()
                    .filter_map(|(timestep_idx, flag_chunk)| {
                        if flag_chunk.iter().any(|&x| !x) {
                            Some(timestep_idx)
                        } else {
                            None
                        }
                    })
                    // now convert this iterator of unflagged timestep ranges to an iterator of contiguous ranges
                    .fold(Vec::new(), |mut acc, timestep_idx| {
                        if acc.is_empty() {
                            acc.push(timestep_idx..(timestep_idx + 1));
                        } else {
                            let mut last_acc = acc.pop().unwrap();
                            if last_acc.end == timestep_idx {
                                last_acc.end = timestep_idx + 1;
                                acc.push(last_acc);
                            } else {
                                acc.push(last_acc);
                                acc.push(timestep_idx..(timestep_idx + 1));
                            }
                        }
                        acc
                    });

                for timestep_range in unflagged_timestep_ranges {
                    let mut jones_chunk = jones_chunk.slice_mut(s![timestep_range.clone(), .., ..]);
                    let mut weight_chunk = weight_chunk.slice_mut(s![timestep_range, .., ..]);
                    if self.correct_van_vleck {
                        trace!("correcting van vleck");
                        // Only correct successfully read data
                        with_increment_duration!(
                            "correct_van_vleck",
                            correct_van_vleck(
                                jones_chunk.view_mut(),
                                &sel_ant_pairs,
                                &flagged_ants,
                                sample_scale,
                            )?
                        );
                    }
                    if self.correct_cable_lengths {
                        trace!("correcting cable lengths");
                        with_increment_duration!(
                            "correct_cable",
                            correct_cable_lengths(
                                corr_ctx,
                                jones_chunk.view_mut(),
                                &coarse_chan_range,
                                &sel_ant_pairs,
                            )?
                        );
                    }
                    if self.correct_digital_gains {
                        trace!("correcting digital gains");
                        with_increment_duration!(
                            "correct_digital",
                            correct_digital_gains(
                                corr_ctx,
                                jones_chunk.view_mut(),
                                &coarse_chan_range,
                                &sel_ant_pairs,
                            )?
                        );
                    }
                    // perform pfb passband gain corrections
                    if let Some(passband_gains) = self.passband_gains {
                        trace!("correcting pfb gains");
                        with_increment_duration!(
                            "correct_passband",
                            correct_coarse_passband_gains(
                                jones_chunk.view_mut(),
                                weight_chunk.view_mut(),
                                passband_gains,
                                fine_chans_per_coarse,
                                &scrunch_type,
                            )?
                        );
                    }
                    preflag_progress.inc(1);
                }
            }
            preflag_progress.finish();
        }

        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                if let Some(strategy) = self.aoflagger_strategy.as_ref() {
                    trace!("using aoflagger");
                    let aoflagger = unsafe { cxx_aoflagger_new() };
                    with_increment_duration!(
                        "flag",
                        flag_jones_array_existing(
                            &aoflagger,
                            strategy,
                            jones_array.view(),
                            flag_array.view_mut(),
                            true,
                            self.draw_progress,
                        )
                    );
                }
            }
        }

        if self.correct_geometry {
            trace!("correcting geometric delays");
            with_increment_duration!(
                "correct_geom",
                correct_geometry(
                    corr_ctx,
                    jones_array.view_mut(),
                    vis_sel,
                    Some(self.array_pos),
                    Some(self.phase_centre),
                    self.draw_progress,
                )
            );
        }

        if let Some(ref calsols) = self.calsols {
            trace!("applying calibration solutions");
            with_increment_duration!(
                "calibrate",
                apply_di_calsol(
                    calsols.view(),
                    jones_array.view_mut(),
                    weight_array.view_mut(),
                    flag_array.view_mut(),
                    &sel_ant_pairs,
                )?
            );
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use float_cmp::F32Margin;
    use marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib::MWAVersion,
    };
    use tempfile::tempdir;

    use crate::{
        flag_to_weight_array,
        flags::get_weight_factor,
        io::{read_mwalib, write_uvfits},
        passband_gains::PFB_JAKE_2022_200HZ,
        test_common::{compare_uvfits_with_csv, get_1254670392_avg_paths},
        FlagContext, VisSelection,
    };

    use super::*;

    #[test]
    #[allow(clippy::field_reassign_with_default)]
    fn test_1254670392_avg_uvfits_no_correct_norfi() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths)
            .expect("unable to get mwalib context");

        let mut prep_ctx = PreprocessContext::default();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        prep_ctx.array_pos = LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        };
        prep_ctx.phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);

        prep_ctx.correct_cable_lengths = false;
        prep_ctx.correct_digital_gains = false;
        prep_ctx.correct_geometry = false;
        prep_ctx.draw_progress = false;

        let (avg_time, avg_freq) = (1, 1);

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        flag_ctx
            .set_flags(
                flag_array.view_mut(),
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )
            .unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        // generate weights
        let weight_factor = get_weight_factor(&corr_ctx);
        let mut weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        prep_ctx
            .preprocess(
                &corr_ctx,
                jones_array.view_mut(),
                weight_array.view_mut(),
                flag_array.view_mut(),
                &vis_sel,
            )
            .unwrap();

        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        write_uvfits(
            &uvfits_path,
            &corr_ctx,
            jones_array.view(),
            weight_array.view(),
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            &vis_sel.baseline_idxs,
            Some(prep_ctx.array_pos),
            Some(prep_ctx.phase_centre),
            avg_time,
            avg_freq,
        )
        .unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
            true,
            false,
        );
    }

    #[test]
    pub fn test_handle_weird_mwa_versions() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let mut corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        let mut prep_ctx = PreprocessContext::default();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        prep_ctx.array_pos = LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        };
        prep_ctx.phase_centre = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);

        prep_ctx.correct_cable_lengths = false;
        prep_ctx.correct_digital_gains = false;
        prep_ctx.correct_geometry = false;
        prep_ctx.draw_progress = false;
        prep_ctx.passband_gains = Some(PFB_JAKE_2022_200HZ);

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();
        let mut weight_array = vis_sel.allocate_weights(fine_chans_per_coarse).unwrap();
        weight_array.fill(get_weight_factor(&corr_ctx) as _);

        corr_ctx.metafits_context.mwa_version = Some(MWAVersion::CorrOldLegacy);

        assert!(prep_ctx
            .preprocess(
                &corr_ctx,
                jones_array.view_mut(),
                weight_array.view_mut(),
                flag_array.view_mut(),
                &vis_sel,
            )
            .is_ok());

        corr_ctx.metafits_context.mwa_version = Some(MWAVersion::VCSLegacyRecombined);

        let result = prep_ctx.preprocess(
            &corr_ctx,
            jones_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &vis_sel,
        );

        assert!(matches!(result, Err(BirliError::BadMWAVersion { .. })));
    }

    #[test]
    #[rustfmt::skip]
    fn test_preprocess_runs_correct_van_vleck() {
        use float_cmp::assert_approx_eq;
        use marlu::mwalib::CorrelatorContext;
        use ndarray::Array3;

        let vis_dims = (1, 1, 3);
        let mut corr_ctx = CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &["tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits"],
        )
        .unwrap();
        corr_ctx.metafits_context.corr_fine_chan_width_hz = 1;
        corr_ctx.metafits_context.corr_int_time_ms = 500;
        corr_ctx.metafits_context.corr_raw_scale_factor = 1.0;
        let sample_scale = get_vv_sample_scale(&corr_ctx).unwrap();
        assert_eq!(sample_scale, 1.0);

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[0, 1]);
        assert_eq!(vis_sel.get_ant_pairs(&corr_ctx.metafits_context), vec![(0, 0), (0, 1), (1, 1)]);

        let mut jones_array = Array3::<Jones<f32>>::zeros(vis_dims);

        // $\hat σ$
        let sighats1: [f64; 2] = [1.453730115943, 1.373255711803];
        let sighats2: [f64; 2] = [1.495798281855, 1.441745903410];
        // $\hat κ$
        let khats1: [f64; 2] = [-0.046568750000, 0.012337500000];
        let khats2: [f64; 2] = [-0.042031250000, -0.011650000000];

        let khatsxx: [f64; 2] = [-0.000037500000, 0.001550000000];
        let khatsyy: [f64; 2] = [-0.008881250000, -0.004287500000];
        let khatsxy: [f64; 2] = [-0.001587500000, 0.009337500000];
        let khatsyx: [f64; 2] = [-0.002425000000, -0.004268750000];
        // $σ$
        let sigmas1: [f64; 2] = [1.424780710577, 1.342571513473];
        let sigmas2: [f64; 2] = [1.467679841384, 1.412550740902];
        // $κ$
        let kappas1: [f64; 2] = [-0.046568781137, 0.012337508611];
        let kappas2: [f64; 2] = [-0.042031317949, -0.011650019325];

        let kappasxx: [f64; 2] = [-0.000037500067, 0.001550002722];
        let kappasyy: [f64; 2] = [-0.008881254122, -0.004287502263];
        let kappasxy: [f64; 2] = [-0.001587501562, 0.009337509186];
        let kappasyx: [f64; 2] = [-0.002425003098, -0.004268755205];

        jones_array[(0, 0, 0)][0].re = sighats1[0].powi(2) as f32;
        jones_array[(0, 0, 0)][1].re = khats1[0] as f32;
        jones_array[(0, 0, 0)][1].im = khats1[1] as f32;
        jones_array[(0, 0, 0)][2].re = khats1[0] as f32;
        jones_array[(0, 0, 0)][2].im = -khats1[1] as f32;
        jones_array[(0, 0, 0)][3].re = sighats1[1].powi(2) as f32;

        jones_array[(0, 0, 1)][0].re = khatsxx[0] as f32;
        jones_array[(0, 0, 1)][0].im = khatsxx[1] as f32;
        jones_array[(0, 0, 1)][1].re = khatsxy[0] as f32;
        jones_array[(0, 0, 1)][1].im = khatsxy[1] as f32;
        jones_array[(0, 0, 1)][2].re = khatsyx[0] as f32;
        jones_array[(0, 0, 1)][2].im = khatsyx[1] as f32;
        jones_array[(0, 0, 1)][3].re = khatsyy[0] as f32;
        jones_array[(0, 0, 1)][3].im = khatsyy[1] as f32;

        jones_array[(0, 0, 2)][0].re = sighats2[0].powi(2) as f32;
        jones_array[(0, 0, 2)][1].re = khats2[0] as f32;
        jones_array[(0, 0, 2)][1].im = khats2[1] as f32;
        jones_array[(0, 0, 2)][2].re = khats2[0] as f32;
        jones_array[(0, 0, 2)][2].im = -khats2[1] as f32;
        jones_array[(0, 0, 2)][3].re = sighats2[1].powi(2) as f32;

        // Set up PreprocessContext to only run Van Vleck
        let mut prep_ctx = PreprocessContext::default();
        prep_ctx.correct_van_vleck = true;
        prep_ctx.correct_cable_lengths = false;
        prep_ctx.correct_digital_gains = false;
        prep_ctx.correct_geometry = false;
        prep_ctx.draw_progress = false;
        prep_ctx.passband_gains = None;
        // phase_centre and array_pos are required but not used for van vleck
        prep_ctx.array_pos = LatLngHeight {
            longitude_rad: 0.0,
            latitude_rad: 0.0,
            height_metres: 0.0,
        };
        prep_ctx.phase_centre = RADec { ra: 0.0, dec: 0.0 };
        println!("prep_ctx: {}", prep_ctx.as_comment());

        // Dummy arrays for weights and flags
        let mut weight_array = Array3::<f32>::ones(vis_dims);
        let mut flag_array = Array3::<bool>::from_elem(vis_dims, false);

        // Call preprocess
        prep_ctx
            .preprocess(
                &corr_ctx,
                jones_array.view_mut(),
                weight_array.view_mut(),
                flag_array.view_mut(),
                &vis_sel,
            )
            .unwrap();

        // Check that the output matches the expected values after correction

        assert_approx_eq!(f32, jones_array[(0, 0, 0)][0].re, sigmas1[0].powi(2) as f32);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][1].re, kappas1[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][1].im, kappas1[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][2].re, kappas1[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][2].im, -kappas1[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][3].re, sigmas1[1].powi(2) as f32);

        assert_approx_eq!(f32, jones_array[(0, 0, 1)][0].re, kappasxx[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][0].im, kappasxx[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][1].re, kappasxy[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][1].im, kappasxy[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][2].re, kappasyx[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][2].im, kappasyx[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][3].re, kappasyy[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 1)][3].im, kappasyy[1] as f32, epsilon=1e-9);

        assert_approx_eq!(f32, jones_array[(0, 0, 2)][0].re, sigmas2[0].powi(2) as f32);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][1].re, kappas2[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][1].im, kappas2[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][2].re, kappas2[0] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][2].im, -kappas2[1] as f32, epsilon=1e-9);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][3].re, sigmas2[1].powi(2) as f32);

    }
}
