//! Crate for preprocessing visibilities
use crate::{
    calibration::apply_di_calsol,
    correct_cable_lengths, correct_geometry,
    corrections::{correct_coarse_passband_gains, correct_digital_gains, ScrunchType},
    marlu::{
        mwalib::CorrelatorContext,
        ndarray::{Array2, Array3},
        Jones, LatLngHeight, RADec,
    },
    with_increment_duration, BirliError, VisSelection,
};
use cfg_if::cfg_if;
use derive_builder::Builder;
use log::info;
use std::{
    collections::HashMap,
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
#[derive(Builder, Debug, Default)]
pub struct PreprocessContext {
    /// The array position used for geometric corrections
    pub array_pos: LatLngHeight,
    /// The phase centre used for geometric corrections
    pub phase_centre: RADec,

    /// Whether cable length corrections are enabled
    #[builder(default = "true")]
    pub correct_cable_lengths: bool,
    /// Whether digital gain corrections are enabled
    #[builder(default = "true")]
    pub correct_digital_gains: bool,
    /// the pfb passband gains to use for corrections
    pub passband_gains: Option<Vec<f64>>,
    /// The calibration solutions to apply
    pub calsols: Option<Array2<Jones<f64>>>,
    /// Whether geometric corrections are enabled
    #[builder(default = "true")]
    pub correct_geometry: bool,

    /// AOFlagger strategy path for flagging
    #[builder(default)]
    #[cfg(feature = "aoflagger")]
    pub aoflagger_strategy: Option<String>,

    /// Whether to draw progress bars
    #[builder(default = "true")]
    pub draw_progress: bool,
}

impl Display for PreprocessContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
                    writeln!(f, "Will flag with aoflagger strategy {}", strategy)?;
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

impl PreprocessContext {
    /// Preprocess visibilities for a chunk of correlator data
    ///
    /// # Arguments
    /// * `corr_ctx` - [`marlu::mwalib::CorrelatorContext`]
    /// * `jones_array` - Array of Jones visibilties
    /// * `weight_array` - Array of weights associated with Jones visibilities
    /// * `flag_array` - Array of flags associated with Jones visibilities
    /// * `durations` - Hashmap used to record timing info
    ///
    /// # Errors
    /// will wrap errors from `correct_digital_gains`, `correct_coarse_passband_gains`
    ///
    /// TODO: more granular error types: `PreprocessingError` -> {`DigitalGainsError`, etc.}
    #[allow(clippy::too_many_arguments)]
    pub fn preprocess(
        &self,
        corr_ctx: &CorrelatorContext,
        jones_array: &mut Array3<Jones<f32>>,
        weight_array: &mut Array3<f32>,
        flag_array: &mut Array3<bool>,
        durations: &mut HashMap<String, Duration>,
        vis_sel: &VisSelection,
    ) -> Result<(), BirliError> {
        if self.correct_cable_lengths {
            with_increment_duration!(
                durations,
                "correct",
                correct_cable_lengths(corr_ctx, jones_array, &vis_sel.coarse_chan_range, false)
            );
        }

        let sel_ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);

        if self.correct_digital_gains {
            with_increment_duration!(
                durations,
                "correct",
                correct_digital_gains(
                    corr_ctx,
                    jones_array,
                    &vis_sel.coarse_chan_range,
                    &sel_ant_pairs,
                )?
            );
        }

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

        // perform pfb passband gain corrections
        if let Some(passband_gains) = self.passband_gains.as_ref() {
            info!("correcting pfb gains");
            with_increment_duration!(
                durations,
                "correct",
                correct_coarse_passband_gains(
                    jones_array,
                    weight_array,
                    passband_gains,
                    fine_chans_per_coarse,
                    &ScrunchType::from_mwa_version(corr_ctx.metafits_context.mwa_version.unwrap())?,
                )?
            );
        }

        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                if let Some(strategy) = self.aoflagger_strategy.as_ref() {
                    let aoflagger = unsafe { cxx_aoflagger_new() };
                    with_increment_duration!(durations,
                        "flag",
                        flag_jones_array_existing(
                            &aoflagger,
                            strategy,
                            jones_array,
                            flag_array,
                            true,
                            self.draw_progress,
                        )
                    );
                }
            }
        }

        if self.correct_geometry {
            // perform geometric delay corrections
            with_increment_duration!(
                durations,
                "correct",
                correct_geometry(
                    corr_ctx,
                    jones_array,
                    &vis_sel.timestep_range,
                    &vis_sel.coarse_chan_range,
                    Some(self.array_pos),
                    Some(self.phase_centre),
                    self.draw_progress,
                )
            );
        }

        if let Some(ref calsols) = self.calsols {
            with_increment_duration!(
                durations,
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
        passband_gains::PFB_JAKE_2022_200HZ,
        test_common::{compare_uvfits_with_csv, get_1254670392_avg_paths},
        write_uvfits, FlagContext, VisSelection,
    };

    use super::*;

    #[test]
    #[allow(clippy::field_reassign_with_default)]
    fn test_1254670392_avg_uvfits_no_correct_norfi() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths)
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
                &mut flag_array,
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )
            .unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        // generate weights
        let weight_factor = get_weight_factor(&corr_ctx);
        let mut weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

        let mut durations = HashMap::new();

        prep_ctx
            .preprocess(
                &corr_ctx,
                &mut jones_array,
                &mut weight_array,
                &mut flag_array,
                &mut durations,
                &vis_sel,
            )
            .unwrap();

        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        write_uvfits(
            uvfits_path.clone(),
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
            prep_ctx.draw_progress,
        )
        .unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    #[test]
    pub fn test_handle_weird_mwa_versions() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let mut corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

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
        prep_ctx.passband_gains = Some(PFB_JAKE_2022_200HZ.to_vec());

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();
        let mut weight_array = vis_sel.allocate_weights(fine_chans_per_coarse).unwrap();
        weight_array.fill(get_weight_factor(&corr_ctx) as _);

        corr_ctx.metafits_context.mwa_version = Some(MWAVersion::CorrOldLegacy);

        assert!(matches!(
            prep_ctx.preprocess(
                &corr_ctx,
                &mut jones_array,
                &mut weight_array,
                &mut flag_array,
                &mut HashMap::new(),
                &vis_sel
            ),
            Ok(_)
        ));

        corr_ctx.metafits_context.mwa_version = Some(MWAVersion::VCSLegacyRecombined);

        let result = prep_ctx.preprocess(
            &corr_ctx,
            &mut jones_array,
            &mut weight_array,
            &mut flag_array,
            &mut HashMap::new(),
            &vis_sel,
        );

        assert!(matches!(result, Err(BirliError::BadMWAVersion { .. })));
    }
}
