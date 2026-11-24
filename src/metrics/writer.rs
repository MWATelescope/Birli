//! Writing all of the metrics of an observation to one FITS file.

use std::path::{Path, PathBuf};

use log::{info, warn};

#[cfg(feature = "aoflagger")]
use super::AOFlagMetricsBuilder;
use super::{AutoMetricsBuilder, CrossMetricsBuilder, EAVILSBuilder, MetricsContext, SSINSBuilder};
use crate::{
    io::error::IOError,
    marlu::{
        fitsio::FitsFile, mwalib::CorrelatorContext, ndarray::ArrayView3, Jones, VisSelection,
    },
    FlagContext,
};

/// Computes all of the quality-assurance metrics (see [`crate::metrics`]) of
/// preprocessed visibilities, one chunk of timesteps at a time, and writes
/// them to one FITS file at the end. This is what `birli --metrics-out` uses.
///
/// The metrics are the same however the timesteps are chunked. The
/// visibilities should not be averaged.
///
/// With the `aoflagger` feature, this includes the flag occupancy, and the
/// SSINS z-scores are flagged with `AOFlagger`'s minimal strategy.
pub struct MetricsWriter {
    path: PathBuf,
    fptr: FitsFile,
    auto: AutoMetricsBuilder,
    /// SSINS is the difference between timesteps, so it needs two of them.
    ssins: Option<SSINSBuilder>,
    eavils: EAVILSBuilder,
    #[cfg(feature = "aoflagger")]
    aoflag: AOFlagMetricsBuilder,
    cross: CrossMetricsBuilder,
}

impl MetricsWriter {
    /// Create a new FITS file at `path`, replacing it if it exists, for the
    /// metrics of the selection `vis_sel` of the observation `corr_ctx`.
    ///
    /// The delay spectra of the cross-correlations use baselines up to
    /// `cross_metrics_cutoff_m` metres long.
    ///
    /// Returns `None`, and does not create the file, if all the timesteps or
    /// all the channels of the selection are flagged. The SSINS metrics are
    /// left out if the selection has fewer than two timesteps.
    ///
    /// # Errors
    ///
    /// Returns an error if the file can't be replaced or created.
    pub fn new(
        path: &Path,
        corr_ctx: &CorrelatorContext,
        vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
        cross_metrics_cutoff_m: f32,
    ) -> Result<Option<Self>, IOError> {
        let timestep_range = &vis_sel.timestep_range;
        let timestep_flags = &flag_ctx.timestep_flags[timestep_range.clone()];
        let chan_flags = &flag_ctx.get_raw_chan_flags(&vis_sel.coarse_chan_range);
        for (what, flags) in [("timesteps", timestep_flags), ("channels", chan_flags)] {
            if flags.iter().all(|&flag| flag) {
                warn!(
                    "Not writing metrics to {}: all {what} are flagged",
                    path.display()
                );
                return Ok(None);
            }
        }
        if timestep_range.len() < 2 {
            warn!(
                "Not writing SSINS metrics to {}: timesteps {timestep_range:?} are fewer than 2",
                path.display(),
            );
        }

        if path.exists() {
            std::fs::remove_file(path)?;
        }
        let fptr = FitsFile::create(path).open()?;

        let metadata = &MetricsContext::from_mwalib(corr_ctx, vis_sel);
        Ok(Some(Self {
            path: path.to_path_buf(),
            fptr,
            auto: AutoMetricsBuilder::new(metadata, timestep_flags, chan_flags),
            ssins: (timestep_range.len() >= 2)
                .then(|| SSINSBuilder::new(metadata, timestep_flags, chan_flags)),
            eavils: EAVILSBuilder::new(metadata, timestep_flags, chan_flags),
            #[cfg(feature = "aoflagger")]
            aoflag: AOFlagMetricsBuilder::new(metadata, timestep_flags, chan_flags),
            cross: CrossMetricsBuilder::new(
                metadata,
                timestep_flags,
                chan_flags,
                cross_metrics_cutoff_m,
            ),
        }))
    }

    /// Add the preprocessed visibilities and flags
    /// `[timestep][channel][baseline]` of the next chunk of timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the channels or baselines are not those of the selection, or
    /// if there are more timesteps than are left in it.
    pub fn push(&mut self, jones_array: ArrayView3<Jones<f32>>, flag_array: ArrayView3<bool>) {
        self.auto.push(jones_array);
        if let Some(ssins) = self.ssins.as_mut() {
            ssins.push(jones_array, flag_array);
        }
        self.eavils.push(jones_array);
        #[cfg(feature = "aoflagger")]
        self.aoflag.push(flag_array);
        self.cross.push(jones_array);
    }

    /// Compute the metrics of all the chunks, and write them to the file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file can't be written.
    ///
    /// # Panics
    ///
    /// Panics if not all of the timesteps of the selection have been pushed.
    pub fn finish(self) -> Result<(), IOError> {
        let mut fptr = self.fptr;
        self.auto.finish().save_to_fits(&mut fptr)?;
        #[cfg(feature = "aoflagger")]
        self.aoflag.finish().save_to_fits(&mut fptr)?;
        self.cross.finish().save_to_fits(&mut fptr)?;
        if let Some(ssins) = self.ssins {
            #[cfg_attr(not(feature = "aoflagger"), allow(unused_mut))]
            let mut ssins = ssins.finish();
            #[cfg(feature = "aoflagger")]
            ssins.flag(None);
            ssins.save_to_fits(&mut fptr)?;
        }
        self.eavils.finish().save_to_fits(&mut fptr)?;
        info!("Saved metrics to {}", self.path.display());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::{io::read_mwalib, marlu::ndarray::s, test_common::get_mwax_data_paths};

    /// The names of the HDUs in a FITS file, after the primary HDU.
    fn hdu_names(path: &Path) -> Vec<String> {
        let mut fptr = FitsFile::open(path).unwrap();
        let num_hdus = fptr.iter().count();
        (1..num_hdus)
            .map(|i| {
                let hdu = fptr.hdu(i).unwrap();
                hdu.read_key::<String>(&mut fptr, "EXTNAME").unwrap()
            })
            .collect()
    }

    /// The 1297526432 mwax dataset, which has 2 tiles, 4 timesteps and 2
    /// coarse channels of 2 fine channels, with nothing flagged.
    fn mwax_context() -> (CorrelatorContext, VisSelection, FlagContext) {
        let (metafits_path, gpufits_paths) = get_mwax_data_paths();
        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        flag_ctx.timestep_flags.fill(false);
        flag_ctx.fine_chan_flags.fill(false);
        (corr_ctx, vis_sel, flag_ctx)
    }

    /// Chunks of timesteps can be pushed one at a time, and the file is
    /// written when the writer is finished.
    #[test]
    fn test_metrics_writer_chunks() {
        let (corr_ctx, vis_sel, flag_ctx) = mwax_context();
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();
        assert_eq!(jones_array.dim().0, 4);

        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("metrics.fits");
        // an existing file is replaced
        std::fs::write(&path, "not a fits file").unwrap();

        let mut writer = MetricsWriter::new(&path, &corr_ctx, &vis_sel, &flag_ctx, 1e4)
            .unwrap()
            .expect("not everything is flagged");
        for chunk in [s![0..1, .., ..], s![1..4, .., ..]] {
            writer.push(jones_array.slice(chunk), flag_array.slice(chunk));
        }
        writer.finish().unwrap();

        let hdus = hdu_names(&path);
        for prefix in [
            "AUTO_POL",
            "AUTO_DELAY_POL",
            "CROSS_DELAY_POL",
            "SSINS_POL",
            "EAVILS_POL",
        ] {
            let extname = format!("{prefix}=XX");
            assert!(hdus.contains(&extname), "{extname} is not in {hdus:?}");
        }
        assert_eq!(
            hdus.contains(&"AO_FLAG_METRICS".to_string()),
            cfg!(feature = "aoflagger")
        );
        // 3 differences between the 4 timesteps, 4 fine channels
        let mut fptr = FitsFile::open(&path).unwrap();
        let hdu = fptr.hdu("SSINS_POL=XX").unwrap();
        let zscores: Vec<f32> = hdu.read_image(&mut fptr).unwrap();
        assert_eq!(zscores.len(), 12);
        assert!(zscores.iter().all(|zscore| zscore.is_finite()));
    }

    /// There are no SSINS metrics for a single timestep.
    #[test]
    fn test_metrics_writer_single_timestep() {
        let (corr_ctx, mut vis_sel, flag_ctx) = mwax_context();
        vis_sel.timestep_range = vis_sel.timestep_range.start..vis_sel.timestep_range.start + 1;
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("metrics.fits");
        let mut writer = MetricsWriter::new(&path, &corr_ctx, &vis_sel, &flag_ctx, 30.0)
            .unwrap()
            .expect("not everything is flagged");
        writer.push(jones_array.view(), flag_array.view());
        writer.finish().unwrap();

        let hdus = hdu_names(&path);
        assert!(hdus.contains(&"EAVILS_POL=XX".to_string()));
        assert!(!hdus.iter().any(|name| name.starts_with("SSINS")));
    }

    /// No file is created if every timestep, or every channel is flagged.
    #[test]
    fn test_metrics_writer_all_flagged() {
        let (corr_ctx, vis_sel, flag_ctx) = mwax_context();
        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("metrics.fits");

        let mut timesteps_flagged = flag_ctx.clone();
        timesteps_flagged.timestep_flags.fill(true);
        let mut chans_flagged = flag_ctx;
        chans_flagged.fine_chan_flags.fill(true);

        for flag_ctx in [timesteps_flagged, chans_flagged] {
            let writer = MetricsWriter::new(&path, &corr_ctx, &vis_sel, &flag_ctx, 30.0).unwrap();
            assert!(writer.is_none());
            assert!(!path.exists());
        }
    }

    /// It is an error if the file can't be created.
    #[test]
    fn test_metrics_writer_unwritable() {
        let (corr_ctx, vis_sel, flag_ctx) = mwax_context();
        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("does_not_exist").join("metrics.fits");

        let result = MetricsWriter::new(&path, &corr_ctx, &vis_sel, &flag_ctx, 30.0);
        assert!(result.is_err());
    }
}
