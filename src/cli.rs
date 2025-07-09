//! Command Line Interface helpers for Birli

#![allow(clippy::doc_markdown)]

use std::{
    convert::Into,
    ffi::OsString,
    fmt::{Debug, Display},
    path::Path,
    time::Duration,
};

use cfg_if::cfg_if;
use clap::{arg, command, ErrorKind::ArgumentNotFound, PossibleValue, ValueHint::FilePath};
use indicatif::{ProgressDrawTarget, ProgressStyle};
use itertools::{izip, Itertools};
use log::{debug, info, trace};
use prettytable::{format as prettyformat, row, table};

use crate::error::{
    BirliError::{self, BadMWAVersion, DryRun},
    CLIError::{InvalidCommandLineArgument, InvalidRangeSpecifier},
};
use crate::flags::FlagContext;
use crate::io::{aocal::AOCalSols, read_mwalib, IOContext};
use crate::marlu::{
    built_info::PKG_VERSION as MARLU_PKG_VERSION,
    constants::{
        COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
    },
    fitsio::FitsFile,
    hifitime::{self, Epoch},
    io::{error::BadArrayShape, ms::MeasurementSetWriter, uvfits::UvfitsWriter, VisWrite},
    mwalib::{
        built_info::PKG_VERSION as MWALIB_PKG_VERSION, CableDelaysApplied, CorrelatorContext,
        GeometricDelaysApplied, MWAVersion, MetafitsContext,
    },
    ndarray::s,
    precession::{precess_time, PrecessionInfo},
    History, Jones, LatLngHeight, MwaObsContext, ObsContext, RADec, VisContext, ENH,
};
use crate::passband_gains::{OSPFB_JAKE_2025_200HZ, PFB_COTTER_2014_10KHZ, PFB_JAKE_2022_200HZ};
use crate::ssins::{EAVILS, SSINS};
use crate::{with_increment_duration, Axis, Complex, FlagFileSet, PreprocessContext, VisSelection};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use aoflagger_sys::{cxx_aoflagger_new};
    }
}

/// Args for preprocessing a correlator context.
pub struct BirliContext<'a> {
    /// `mwalib::CorrelatorContext`
    pub corr_ctx: CorrelatorContext,
    /// Preprocessing parameters
    pub prep_ctx: PreprocessContext<'a>,
    /// selected visibility indices
    pub vis_sel: VisSelection,
    /// Flagging Parameters
    pub flag_ctx: FlagContext,
    /// Input / output paths
    pub io_ctx: IOContext,
    /// temporal averaging factor
    pub avg_time: usize,
    /// spectral averaging factor
    pub avg_freq: usize,
    /// temporal chunking factor
    pub num_timesteps_per_chunk: Option<usize>,
    /// channel selections for picket-fencing
    pub channel_range_sel: ChannelRanges,
}

// Add build-time information from the "built" crate.
include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// stolen from hyperdrive
/// Write many info-level log lines of how this executable was compiled.
///
/// # Errors
///
/// propagates writeln! fails
pub fn fmt_build_info(f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match GIT_HEAD_REF {
        Some(hr) => {
            let dirty = GIT_DIRTY.unwrap_or(false);
            writeln!(
                f,
                "Compiled on git commit hash: {}{}",
                GIT_COMMIT_HASH.unwrap(),
                if dirty { " (dirty)" } else { "" }
            )?;
            writeln!(f, "            git head ref: {hr}")?;
        }
        None => writeln!(f, "Compiled on git commit hash: <no git info>")?,
    }
    writeln!(f, "            {BUILT_TIME_UTC}")?;
    writeln!(f, "         with compiler {RUSTC_VERSION}")?;
    writeln!(f, "libraries:")?;
    writeln!(f, "- marlu v{MARLU_PKG_VERSION}")?;
    writeln!(f, "- mwalib v{MWALIB_PKG_VERSION}")?;

    cfg_if! {
        if #[cfg(feature = "aoflagger")] {
            use std::os::raw::c_short;
            let mut major: c_short = -1;
            let mut minor: c_short = -1;
            let mut sub_minor: c_short = -1;
            let aoflagger = unsafe { cxx_aoflagger_new() };
            aoflagger.GetVersion(&mut major, &mut minor, &mut sub_minor);
            assert!(major >= 3);
            assert!(minor >= 0);
            assert!(sub_minor >= 0);
            writeln!(f, "- aoflagger v{major}.{minor}.{sub_minor}")?;
        }
    }
    writeln!(f)?;
    Ok(())
}

fn time_details(
    gps_time_ms: u64,
    dut1: hifitime::Duration,
    phase_centre: RADec,
    array_pos: LatLngHeight,
) -> (String, String, f64, PrecessionInfo) {
    let epoch = Epoch::from_gpst_seconds(gps_time_ms as f64 / 1e3);
    let (y, mo, d, h, mi, s, ms) = epoch.to_gregorian_utc();
    let precession_info = precess_time(
        array_pos.longitude_rad,
        array_pos.latitude_rad,
        phase_centre,
        epoch,
        dut1,
    );
    (
        format!("{y:02}-{mo:02}-{d:02}"),
        format!(
            "{:02}:{:02}:{:02}.{:03}",
            h,
            mi,
            s,
            (ms as f64 / 1e6).round()
        ),
        epoch.to_mjd_utc_seconds(),
        precession_info,
    )
}

/// Structure representing a series of channel ranges
/// e.g. 1-10, 20-30, 40-50
#[derive(Debug)]
pub struct ChannelRanges {
    /// Vector of channel ranges
    pub ranges: Vec<(usize, usize)>,
}

impl ChannelRanges {
    /// Create a new `ChannelRanges` object from a string
    /// e.g. "1-10, 20-30, 40-50"
    ///
    /// # Errors
    /// `BirliError::CLIError` if the string does not conform to the expected format
    pub fn new(s: &str) -> Result<Self, BirliError> {
        let mut ranges = Vec::new();
        for range in s.split(',') {
            match range.split('-').collect::<Vec<_>>().as_slice() {
                [start, end] => {
                    let start_result = start.trim().parse::<usize>();
                    let end_result = end.trim().parse::<usize>();
                    match (start_result, end_result) {
                        (Ok(start), Ok(end)) => {
                            ranges.push((start, end));
                        }
                        _ => {
                            return Err(BirliError::CLIError(InvalidRangeSpecifier {
                                reason: format!("invalid channel range: {}", range),
                            }));
                        }
                    }
                }
                [start] => {
                    let start_result = start.trim().parse::<usize>();
                    match start_result {
                        Ok(start) => {
                            ranges.push((start, start));
                        }
                        _ => {
                            return Err(BirliError::CLIError(InvalidRangeSpecifier {
                                reason: format!("invalid channel range: {}", range),
                            }));
                        }
                    }
                }
                _ => {
                    return Err(BirliError::CLIError(InvalidRangeSpecifier {
                        reason: format!("invalid channel range: {}", range),
                    }));
                }
            }
        }
        Ok(Self { ranges })
    }

    /// create coarse channel ranges from a slice of indices
    pub fn from_idxs(context: &CorrelatorContext, coarse_chan_indices: &[usize]) -> Self {
        // get the first value in the vector
        let mut range_start_idx = coarse_chan_indices[0];
        let mut range_end_idx = range_start_idx;
        let mut range_end_ch = context.coarse_chans[range_start_idx].rec_chan_number;

        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for &idx in &coarse_chan_indices[1..] {
            let ch = context.coarse_chans[idx].rec_chan_number;
            if idx == range_end_idx + 1 && ch == range_end_ch + 1 {
                range_end_idx = idx;
                range_end_ch = ch;
            } else {
                ranges.push((range_start_idx, range_end_idx));
                range_start_idx = idx;
                range_end_idx = range_start_idx;
                range_end_ch = ch;
            }
        }
        ranges.push((range_start_idx, range_end_idx));
        Self { ranges }
    }

    /// Create a new `ChannelRanges` object spanning all available channels in the metafits context
    pub fn all(context: &CorrelatorContext) -> Self {
        let coarse_chan_indices = &(0..context.num_coarse_chans).collect_vec();
        Self::from_idxs(context, coarse_chan_indices)
    }

    /// Create a new `ChannelRanges` object spanning all available channels in the metafits context
    pub fn provided(context: &CorrelatorContext) -> Self {
        Self::from_idxs(context, &context.provided_coarse_chan_indices)
    }
}

impl Display for BirliContext<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{PKG_NAME} version {PKG_VERSION}")?;

        fmt_build_info(f)?;

        writeln!(
            f,
            "observation name:     {}",
            self.corr_ctx.metafits_context.obs_name
        )?;

        writeln!(f, "Array position:       {}", &self.prep_ctx.array_pos)?;
        writeln!(f, "Phase centre:         {}", &self.prep_ctx.phase_centre)?;
        let pointing_centre = RADec::from_mwalib_tile_pointing(&self.corr_ctx.metafits_context);
        if pointing_centre != self.prep_ctx.phase_centre {
            writeln!(f, "Pointing centre:      {}", &pointing_centre)?;
        }

        let coarse_chan_flag_idxs: Vec<usize> = self
            .flag_ctx
            .coarse_chan_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        // TODO: actually display this.
        let _fine_chan_flag_idxs: Vec<usize> = self
            .flag_ctx
            .fine_chan_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        let timestep_flag_idxs: Vec<usize> = self
            .flag_ctx
            .timestep_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        let ant_pairs = self.vis_sel.get_ant_pairs(&self.corr_ctx.metafits_context);
        #[allow(clippy::needless_collect)]
        let baseline_flag_idxs: Vec<usize> = self
            .flag_ctx
            .get_baseline_flags(&ant_pairs)
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();

        let dut1 =
            hifitime::Duration::from_seconds(self.corr_ctx.metafits_context.dut1.unwrap_or(0.0));
        let (sched_start_date, sched_start_time, sched_start_mjd_s, sched_start_prec) =
            time_details(
                self.corr_ctx.metafits_context.sched_start_gps_time_ms,
                dut1,
                self.prep_ctx.phase_centre,
                self.prep_ctx.array_pos,
            );
        writeln!(
            f,
            "Scheduled start:      {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
            sched_start_date, sched_start_time,
            self.corr_ctx.metafits_context.sched_start_unix_time_ms as f64 / 1e3,
            self.corr_ctx.metafits_context.sched_start_gps_time_ms as f64 / 1e3,
            sched_start_mjd_s,
            sched_start_prec.lmst.to_degrees(),
            sched_start_prec.lmst_j2000.to_degrees(),
            sched_start_prec.array_latitude_j2000.to_degrees(),
        )?;
        let (sched_end_date, sched_end_time, sched_end_mjd_s, sched_end_prec) = time_details(
            self.corr_ctx.metafits_context.sched_end_gps_time_ms,
            dut1,
            self.prep_ctx.phase_centre,
            self.prep_ctx.array_pos,
        );
        writeln!(
            f,

            "Scheduled end:        {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
            sched_end_date, sched_end_time,
            self.corr_ctx.metafits_context.sched_end_unix_time_ms as f64 / 1e3,
            self.corr_ctx.metafits_context.sched_end_gps_time_ms as f64 / 1e3,
            sched_end_mjd_s,
            sched_end_prec.lmst.to_degrees(),
            sched_end_prec.lmst_j2000.to_degrees(),
            sched_end_prec.array_latitude_j2000.to_degrees(),
        )?;
        let int_time_s = self.corr_ctx.metafits_context.corr_int_time_ms as f64 / 1e3;
        let sched_duration_s = self.corr_ctx.metafits_context.sched_duration_ms as f64 / 1e3;
        writeln!(
            f,
            "Scheduled duration:   {:.3}s = {:3} * {:.3}s",
            sched_duration_s,
            (sched_duration_s / int_time_s).ceil(),
            int_time_s
        )?;
        let quack_duration_s = self.corr_ctx.metafits_context.quack_time_duration_ms as f64 / 1e3;
        writeln!(
            f,
            "Quack duration:       {:.3}s = {:3} * {:.3}s",
            quack_duration_s,
            (quack_duration_s / int_time_s).ceil(),
            int_time_s
        )?;
        let num_avg_timesteps =
            (self.vis_sel.timestep_range.len() as f64 / self.avg_time as f64).ceil() as usize;
        let avg_int_time_s = int_time_s * self.avg_time as f64;
        writeln!(
            f,
            "Output duration:      {:.3}s = {:3} * {:.3}s{}",
            num_avg_timesteps as f64 * avg_int_time_s,
            num_avg_timesteps,
            avg_int_time_s,
            if self.avg_time == 1 {
                "".into()
            } else {
                format!(" ({}x)", self.avg_time)
            }
        )?;

        let total_bandwidth_mhz = self.corr_ctx.metafits_context.obs_bandwidth_hz as f64 / 1e6;
        let fine_chan_width_khz =
            self.corr_ctx.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;
        let fine_chans_per_coarse = self
            .corr_ctx
            .metafits_context
            .num_corr_fine_chans_per_coarse;

        writeln!(
            f,
            "Scheduled Bandwidth:  {:.3}MHz = {:3} * {:3} * {:.3}kHz",
            total_bandwidth_mhz,
            self.corr_ctx.metafits_context.num_metafits_coarse_chans,
            fine_chans_per_coarse,
            fine_chan_width_khz
        )?;

        let out_bandwidth_mhz = self.vis_sel.coarse_chan_range.len() as f64
            * fine_chans_per_coarse as f64
            * fine_chan_width_khz
            / 1e3;
        let num_avg_chans = (self.vis_sel.coarse_chan_range.len() as f64
            * fine_chans_per_coarse as f64
            / self.avg_freq as f64)
            .ceil() as usize;
        let avg_fine_chan_width_khz = fine_chan_width_khz * self.avg_freq as f64;
        writeln!(
            f,
            "Output Bandwidth:     {:.3}MHz = {:9} * {:.3}kHz{}",
            out_bandwidth_mhz,
            num_avg_chans,
            avg_fine_chan_width_khz,
            if self.avg_freq == 1 {
                "".into()
            } else {
                format!(" ({}x)", self.avg_freq)
            }
        )?;

        let first_epoch =
            Epoch::from_gpst_seconds(self.corr_ctx.timesteps[0].gps_time_ms as f64 / 1e3);
        let (y, mo, d, ..) = first_epoch.to_gregorian_utc();

        let mut timestep_table = table!([
            "",
            format!("{y:02}-{mo:02}-{d:02} UTC +"),
            "unix [s]",
            "gps [s]",
            "p",
            "c",
            "g",
            "s",
            "f"
        ]);
        timestep_table.set_format(*prettyformat::consts::FORMAT_CLEAN);

        let provided_timestep_indices = &self.corr_ctx.provided_timestep_indices;
        let common_timestep_indices = &self.corr_ctx.common_timestep_indices;
        let common_good_timestep_indices = &self.corr_ctx.common_good_timestep_indices;
        for (timestep_idx, timestep) in self.corr_ctx.timesteps.iter().enumerate() {
            let provided = provided_timestep_indices.contains(&timestep_idx);
            let selected = self.vis_sel.timestep_range.contains(&timestep_idx);
            let common = common_timestep_indices.contains(&timestep_idx);
            let good = common_good_timestep_indices.contains(&timestep_idx);
            let flagged = timestep_flag_idxs.contains(&timestep_idx);

            let (_, time, ..) = time_details(
                timestep.gps_time_ms,
                dut1,
                self.prep_ctx.phase_centre,
                self.prep_ctx.array_pos,
            );
            let row = row![r =>
                format!("ts{timestep_idx}:"),
                time,
                format!("{:.3}", timestep.unix_time_ms as f64 / 1e3),
                format!("{:.3}", timestep.gps_time_ms as f64 / 1e3),
                if provided {"p"} else {""},
                if common {"c"} else {""},
                if good {"g"} else {""},
                if selected {"s"} else {""},
                if flagged {"f"} else {""}
            ];
            timestep_table.add_row(row);
        }

        writeln!(
            f,
            "Timestep details (all={}, provided={}, common={}, good={}, select={}, flag={}):\n{}",
            self.corr_ctx.num_timesteps,
            self.corr_ctx.num_provided_timesteps,
            self.corr_ctx.num_common_timesteps,
            self.corr_ctx.num_common_good_timesteps,
            self.vis_sel.timestep_range.len(),
            timestep_flag_idxs.len(),
            timestep_table
        )?;

        let mut coarse_chan_table = table!([
            "",
            "gpu",
            "corr",
            "rec",
            "cen [MHz]",
            "p",
            "c",
            "g",
            "s",
            "f",
            "range",
        ]);
        coarse_chan_table.set_format(*prettyformat::consts::FORMAT_CLEAN);
        // coarse_chan_table
        let provided_coarse_chan_indices = &self.corr_ctx.provided_coarse_chan_indices;
        let common_coarse_chan_indices = &self.corr_ctx.common_coarse_chan_indices;
        let common_good_coarse_chan_indices = &self.corr_ctx.common_good_coarse_chan_indices;
        for (chan_idx, chan) in self.corr_ctx.coarse_chans.iter().enumerate() {
            let provided = provided_coarse_chan_indices.contains(&chan_idx);
            let mut selected = self.vis_sel.coarse_chan_range.contains(&chan_idx);
            let common = common_coarse_chan_indices.contains(&chan_idx);
            let good = common_good_coarse_chan_indices.contains(&chan_idx);
            let mut flagged = coarse_chan_flag_idxs.contains(&chan_idx);
            let range = self
                .channel_range_sel
                .ranges
                .iter()
                .enumerate()
                .find(|(_idx, (start, end))| chan_idx >= *start && chan_idx <= *end);
            selected &= range.is_some();
            flagged &= range.is_some();
            let row = row![r =>
                format!("cc{chan_idx}:"),
                chan.gpubox_number,
                chan.corr_chan_number,
                chan.rec_chan_number,
                format!("{:.4}", chan.chan_centre_hz as f64 / 1e6),
                if provided {"p"} else {""},
                if common {"c"} else {""},
                if good {"g"} else {""},
                if selected {"s"} else {""},
                if flagged {"f"} else {""},
                if let Some((idx, (start,end))) = range { format!("{} ({}-{})", idx, start, end) } else { "".into() }
            ];
            coarse_chan_table.add_row(row);
        }

        writeln!(
            f,
            "Coarse channel details (metafits={}, provided={}, common={}, good={}, select={}, flag={}, ranges={}):\n{}",
            self.corr_ctx.num_coarse_chans,
            self.corr_ctx.num_provided_coarse_chans,
            self.corr_ctx.num_common_coarse_chans,
            self.corr_ctx.num_common_good_coarse_chans,
            self.vis_sel.coarse_chan_range.len(),
            coarse_chan_flag_idxs.len(),
            self.channel_range_sel.ranges.len(),
            coarse_chan_table
        )?;

        writeln!(
            f,
            "Antenna details (all={}, flag={})",
            self.corr_ctx.metafits_context.num_ants,
            self.flag_ctx
                .antenna_flags
                .iter()
                .enumerate()
                .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
                .count(),
            // format!("\n{}", ant_table)
        )?;

        writeln!(
            f,
            "Baseline Details (all={}, auto={}, select={}, flag={}):",
            self.corr_ctx.metafits_context.num_baselines,
            self.corr_ctx.metafits_context.num_ants,
            self.vis_sel.baseline_idxs.len(),
            baseline_flag_idxs.len(),
        )?;

        // TODO: show free memory with https://docs.rs/sys-info/latest/sys_info/fn.mem_info.html

        let num_sel_timesteps = self.vis_sel.timestep_range.len();
        let num_sel_chans = self.vis_sel.coarse_chan_range.len() * fine_chans_per_coarse;
        let num_sel_baselines = self.vis_sel.baseline_idxs.len();
        let num_sel_pols = self.corr_ctx.metafits_context.num_visibility_pols;
        let mem_selected_bytes = self.vis_sel.estimate_bytes_best(fine_chans_per_coarse);
        let mem_per_timestep_gib =
            mem_selected_bytes as f64 / num_sel_timesteps as f64 / 1024.0_f64.powi(3);

        writeln!(
            f,
            "Estimated memory usage per timestep =           {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_chans,
            num_sel_baselines,
            std::mem::size_of::<Jones<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib,
        )?;

        if let Some(num_timesteps) = self.num_timesteps_per_chunk {
            writeln!(
                f,
                "Estimated memory per chunk          = {:5}ts * {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
                num_timesteps,
                num_sel_chans,
                num_sel_baselines,
                std::mem::size_of::<Jones<f32>>(),
                std::mem::size_of::<f32>(),
                std::mem::size_of::<bool>(),
                mem_per_timestep_gib * num_timesteps as f64,
            )?;
        }

        writeln!(
            f,
            "Estimated memory selected           = {:5}ts * {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_timesteps,
            num_sel_chans,
            num_sel_baselines,
            std::mem::size_of::<Jones<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib * num_sel_timesteps as f64,
        )?;

        let avg_mem_per_timestep_gib = (num_avg_chans
            * num_sel_baselines
            * num_sel_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())) as f64
            / 1024.0_f64.powi(3);

        writeln!(
            f,
            "Estimated output size               = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_avg_timesteps,
            num_avg_chans,
            num_sel_baselines,
            num_sel_pols,
            std::mem::size_of::<Complex<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            avg_mem_per_timestep_gib * num_avg_timesteps as f64,
        )?;

        writeln!(f, "Preprocessing Context: \n{}", &self.prep_ctx)?;

        Ok(())
    }
}

impl<'a> BirliContext<'a> {
    // TODO: try struct instead of builder
    #[allow(clippy::cognitive_complexity)]
    fn get_matches<I, T>(args: I) -> Result<clap::ArgMatches, BirliError>
    where
        I: IntoIterator<Item = T> + Debug,
        T: Into<OsString> + Clone,
    {
        #[allow(unused_mut)]
        let mut app = command!()
            .subcommand_precedence_over_arg(true)
            .arg_required_else_help(true)
            .next_line_help(false)
            .about("Preprocess Murchison Widefield Array MetaFITS and GPUFITS data \
                    into usable astronomy formats.")
            .args(&[
                // input options
                arg!(-m --metafits <PATH> "Metadata file for the observation")
                    .required(true)
                    .value_hint(FilePath)
                    .help_heading("INPUT"),
                arg!(fits_paths: <PATHS>... "GPUBox files to process")
                    .help_heading("INPUT")
                    .value_hint(FilePath)
                    .required(true),

                // processing options
                arg!(--"phase-centre" "Override Phase centre from metafits (degrees)")
                    .value_names(&["RA", "DEC"])
                    .allow_hyphen_values(true)
                    .required(false),
                arg!(--"pointing-centre" "Use pointing instead phase centre")
                    .conflicts_with("phase-centre"),
                arg!(--"emulate-cotter" "Use Cotter's array position, not MWAlib's"),
                arg!(--"dry-run" "Just print the summary and exit"),
                arg!(--"no-draw-progress" "do not show progress bars"),

                // selection options
                arg!(--"sel-time" "Timestep index range (inclusive) to select")
                    .help_heading("SELECTION")
                    .value_names(&["MIN", "MAX"])
                    .required(false),
                arg!(--"sel-ants" <ANTS>... "Antenna indices to select")
                    .help_heading("SELECTION")
                    .multiple_values(true)
                    .required(false),
                arg!(--"no-sel-flagged-ants" "Deselect flagged antennas")
                    .help_heading("SELECTION"),
                arg!(--"no-sel-autos" "Deselect autocorrelations")
                    .help_heading("SELECTION"),

                arg!(--"sel-chan-ranges" <RANGES> "Select separate channel ranges")
                    .help_heading("SELECTION")
                    .required(false),
                arg!(--"provided-chan-ranges" "Only consider provided channels")
                    .help_heading("SELECTION")
                    .required(false),

                // resource limit options
                arg!(--"time-chunk" <STEPS> "Process observation in chunks of <STEPS> timesteps.")
                    .help_heading("RESOURCE LIMITS")
                    .required(false)
                    .conflicts_with("max-memory"),
                arg!(--"max-memory" <GIBIBYTES> "Estimate --time-chunk with <GIBIBYTES> GiB each chunk.")
                    .help_heading("RESOURCE LIMITS")
                    .required(false),

                // flagging options
                // -> timesteps
                arg!(--"flag-init" <SECONDS> "Flag <SECONDS> after first common time (quack time)")
                    .alias("--quack-time")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-init-steps" <COUNT> "Flag <COUNT> steps after first common time")
                    .help_heading("FLAGGING")
                    .required(false)
                    .conflicts_with("flag-init"),
                arg!(--"flag-end" <SECONDS> "Flag seconds before the last provided time")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-end-steps" <COUNT> "Flag <COUNT> steps before the last provided")
                    .help_heading("FLAGGING")
                    .required(false)
                    .conflicts_with("flag-end"),
                arg!(--"flag-times" <STEPS>... "Flag additional time steps")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                // -> channels
                arg!(--"flag-coarse-chans" <CHANS> ... "Flag additional coarse chan indices")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                arg!(--"flag-edge-width" <KHZ> "Flag bandwidth [kHz] at the ends of each coarse chan")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-edge-chans" <COUNT> "Flag <COUNT> fine chans on the ends of each coarse")
                    .help_heading("FLAGGING")
                    .conflicts_with("flag-edge-width")
                    .required(false),
                arg!(--"flag-fine-chans" <CHANS>... "Flag fine chan indices in each coarse chan")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                arg!(--"flag-dc" "Force flagging of DC centre chans")
                    .help_heading("FLAGGING")
                    .conflicts_with("no-flag-dc"),
                arg!(--"no-flag-dc" "Do not flag DC centre chans")
                    .help_heading("FLAGGING")
                    .conflicts_with("flag-dc"),
                // -> antennas
                arg!(--"no-flag-metafits" "Ignore antenna flags in metafits")
                    .help_heading("FLAGGING"),
                arg!(--"flag-antennas" <ANTS>... "Flag antenna indices")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                // -> baselines
                arg!(--"flag-autos" "Flag auto correlations")
                    .help_heading("FLAGGING"),

                // corrections
                arg!(--"van-vleck" "Apply Van Vleck corrections")
                    .help_heading("CORRECTION"),
                arg!(--"no-cable-delay" "Do not perform cable length corrections")
                    .help_heading("CORRECTION"),
                arg!(--"no-geometric-delay" "Do not perform geometric corrections")
                    .help_heading("CORRECTION")
                    .alias("no-geom")
                    .conflicts_with("pointing-centre")
                    .conflicts_with("phase-centre"),
                arg!(--"no-digital-gains" "Do not perform digital gains corrections")
                    .help_heading("CORRECTION"),
                arg!(--"passband-gains" <TYPE> "Type of PFB passband filter gains correction to apply")
                    .required(false)
                    .possible_values([
                        PossibleValue::new("none").help("No passband gains correction (unitary)"),
                        PossibleValue::new("cotter")
                            .help(
                                "_sb128ChannelSubbandValue2014FromMemo from
                                    subbandpassband.cpp in Cotter. Can only be used with resolutions of
                                    n * 10kHz"
                            ),
                        PossibleValue::new("jake")
                            .help("see: PFB_JAKE_2022_200HZ in src/passband_gains.rs"),
                        PossibleValue::new("jake_oversampled")
                            .help("see: OSPFB_JAKE_2025_200HZ in src/passband_gains.rs"),
                        PossibleValue::new("auto")
                            .help("MWAX => jake or jake_oversampled, legacy => cotter"),
                    ])
                    .default_value("auto")
                    .alias("pfb-gains")
                    .help_heading("CORRECTION"),

                // calibration
                arg!(--"apply-di-cal" <PATH> "Apply DI calibration solutions before averaging")
                    .required(false)
                    .value_hint(FilePath),

                // averaging
                arg!(--"avg-time-res" <SECONDS> "Time resolution of averaged data")
                    .help_heading("AVERAGING")
                    .required(false),
                arg!(--"avg-time-factor" <FACTOR> "Average <FACTOR> timesteps per averaged timestep")
                    .help_heading("AVERAGING")
                    .required(false)
                    .conflicts_with("avg-time-res"),
                arg!(--"avg-freq-res" <KHZ> "Frequency resolution of averaged data")
                    .help_heading("AVERAGING")
                    .required(false),
                arg!(--"avg-freq-factor" <FACTOR> "Average <FACTOR> channels per averaged channel")
                    .help_heading("AVERAGING")
                    .required(false)
                    .conflicts_with("avg-freq-res"),

                // output options
                arg!(-f --"flag-template" <TEMPLATE> "The template used to name flag files. \
                        Percents are substituted for the zero-prefixed GPUBox ID, which can be up to \
                        3 characters long. Example: FlagFile%%%.mwaf")
                    .help_heading("OUTPUT")
                    .required(false),
                arg!(-u --"uvfits-out" <PATH> "Path for uvfits output")
                    .help_heading("OUTPUT")
                    .required(false),
                arg!(-M --"ms-out" <PATH> "Path for measurement set output")
                    .help_heading("OUTPUT")
                    .required(false),
                arg!(--"metrics-out" <PATH> "Path for csv metrics output")
                    .help_heading("OUTPUT")
                    .required(false),
            ]);
        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                app = app.args(&[
                    arg!(--"no-rfi" "Do not perform RFI Flagging with aoflagger")
                        .help_heading("AOFLAGGER"),
                    arg!(--"aoflagger-strategy" <PATH> "Strategy to use for RFI Flagging")
                        .value_hint(FilePath)
                        .help_heading("AOFLAGGER")
                        .required(false)
                ]);
            }
        };
        let matches = app.try_get_matches_from_mut(args)?;
        Ok(matches)
    }

    fn parse_io_matches(matches: &clap::ArgMatches) -> IOContext {
        IOContext {
            metafits_in: matches
                .value_of_t("metafits")
                .unwrap_or_else(|_| panic!("--metafits <PATH> is required, enforced by clap")),
            gpufits_in: matches
                .values_of_t("fits_paths")
                .unwrap_or_else(|_| panic!("<PATHS> is required, enforced by clap")),
            aocalsols_in: matches.value_of("apply-di-cal").map(Into::into),
            uvfits_out: matches.value_of("uvfits-out").map(Into::into),
            ms_out: matches.value_of("ms-out").map(Into::into),
            flag_template: matches.value_of("flag-template").map(Into::into),
            metrics_out: matches.value_of("metrics-out").map(Into::into),
        }
    }

    fn parse_vis_sel_matches(
        corr_ctx: &CorrelatorContext,
        flag_ctx: &FlagContext,
        matches: &clap::ArgMatches,
    ) -> Result<VisSelection, BirliError> {
        let mut vis_sel = VisSelection::from_mwalib(corr_ctx).unwrap();
        let CorrelatorContext {
            metafits_context: meta_ctx,
            num_timesteps,
            ..
        } = corr_ctx;
        let MetafitsContext { num_ants, .. } = meta_ctx;
        match matches
            .values_of_t::<usize>("sel-time")
            .map(|v| (v[0], v[1]))
        {
            Ok((from, to)) => {
                if from > to || to >= *num_timesteps {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--sel-time <FROM> <TO>".into(),
                        expected: format!("from <= to < num_timesteps={}", num_timesteps),
                        received: format!("from={from} to={to}"),
                    }));
                }
                vis_sel.timestep_range = from..(to + 1);
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        }
        match matches.values_of_t::<usize>("sel-ants") {
            Ok(antenna_idxs) => {
                for (value_idx, &antenna_idx) in antenna_idxs.iter().enumerate() {
                    if antenna_idx >= *num_ants {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--sel-ants <ANTS>...".into(),
                            expected: format!("antenna_idx < num_ants={}", num_ants),
                            received: format!(
                                "antenna_idxs[{value_idx}]={antenna_idx}. all:{antenna_idxs:?}"
                            ),
                        }));
                    }
                }
                // filter vis_sel.baseline_idxs that correspond with antennas not in antenna_idxs
                vis_sel.retain_antennas(meta_ctx, &antenna_idxs);
                // if no baselines are selected, return an error
                if vis_sel.baseline_idxs.is_empty() {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--sel-ants <ANTS>...".into(),
                        expected: "at least one baseline matched".into(),
                        received: format!("antenna_idxs={antenna_idxs:?}"),
                    }));
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        }
        if matches.is_present("no-sel-flagged-ants") {
            let flagged_antenna_idxs = flag_ctx.get_flagged_antenna_idxs();
            vis_sel.filter_antennas(meta_ctx, &flagged_antenna_idxs);
            if vis_sel.baseline_idxs.is_empty() {
                return Err(BirliError::CLIError(InvalidCommandLineArgument {
                    option: "--no-sel-flagged-ants".into(),
                    expected: "at least one baseline matched".into(),
                    received: "".into(),
                }));
            }
        }
        if matches.is_present("no-sel-autos") {
            vis_sel.filter_autos(meta_ctx);
            if vis_sel.baseline_idxs.is_empty() {
                return Err(BirliError::CLIError(InvalidCommandLineArgument {
                    option: "--no-sel-autos".into(),
                    expected: "at least one baseline matched".into(),
                    received: "".into(),
                }));
            }
        }
        Ok(vis_sel)
    }

    fn parse_sel_chan_ranges(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
    ) -> Result<ChannelRanges, BirliError> {
        let all_chan_ranges = ChannelRanges::all(corr_ctx);
        let provided_chan_ranges = ChannelRanges::provided(corr_ctx);
        match (matches.is_present("sel-chan-ranges"), matches.is_present("provided-chan-ranges")) {
            (true, true) => panic!("can't use both --provided-chan-ranges and --sel-chan-ranges"),
            (true, _) => {
                #[allow(clippy::option_if_let_else)]
                match matches.value_of("sel-chan-ranges") {
                    Some(range_str) => ChannelRanges::new(range_str),
                    None => Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--sel-chan-ranges <RANGES>".into(),
                        expected: "comma-separated ranges indexing the metafits coarse channels, e.g. 0-10,20-30".into(),
                        received: "no value".into(),
                    })),
                }
            },
            (_, true) => Ok(provided_chan_ranges),
            (_, _) => Ok(all_chan_ranges),
        }
    }

    fn parse_flag_matches(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
    ) -> Result<FlagContext, BirliError> {
        let mut flag_ctx = FlagContext::from_mwalib(corr_ctx);
        let CorrelatorContext {
            metafits_context: meta_ctx,
            num_timesteps,
            num_coarse_chans,
            ..
        } = corr_ctx;
        let MetafitsContext {
            num_corr_fine_chans_per_coarse: fine_chans_per_coarse,
            corr_fine_chan_width_hz,
            corr_int_time_ms,
            num_ants,
            ..
        } = meta_ctx;
        match matches.values_of_t::<usize>("flag-times") {
            Ok(timestep_idxs) => {
                for (value_idx, &timestep_idx) in timestep_idxs.iter().enumerate() {
                    if timestep_idx >= *num_timesteps {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-times <TIMESTEPS>...".into(),
                            expected: format!("timestep_idx < num_timesteps={}", *num_timesteps),
                            received: format!(
                                "timestep_idxs[{value_idx}]={timestep_idx}. all:{timestep_idxs:?}"
                            ),
                        }));
                    }
                    flag_ctx.timestep_flags[timestep_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.values_of_t::<usize>("flag-coarse-chans") {
            Ok(coarse_chan_idxs) => {
                for (value_idx, &coarse_chan_idx) in coarse_chan_idxs.iter().enumerate() {
                    if coarse_chan_idx >= *num_coarse_chans {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-coarse-chans <CHANS>...".into(),
                            expected: format!(
                                "coarse_chan_idx < num_coarse_chans={}",
                                num_coarse_chans
                            ),
                            received: format!(
                                "coarse_chan_idxs[{value_idx}]={coarse_chan_idx}. all:{coarse_chan_idxs:?}"
                            ),
                        }));
                    }
                    flag_ctx.coarse_chan_flags[coarse_chan_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.values_of_t::<usize>("flag-fine-chans") {
            Ok(fine_chan_idxs) => {
                for (value_idx, &fine_chan_idx) in fine_chan_idxs.iter().enumerate() {
                    if fine_chan_idx >= *fine_chans_per_coarse {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-fine-chans <CHANS>...".into(),
                            expected: format!(
                                "fine_chan_idx < num_fine_chans={fine_chans_per_coarse}"
                            ),
                            received: format!(
                                "fine_chan_idxs[{value_idx}]={fine_chan_idx}. all:{fine_chan_idxs:?}"
                            ),
                        }));
                    }
                    flag_ctx.fine_chan_flags[fine_chan_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        if matches.is_present("no-flag-metafits") {
            info!("Ignoring antenna flags from metafits.");
            // set antenna flags to all false
            flag_ctx.antenna_flags = vec![false; flag_ctx.antenna_flags.len()];
        }
        match matches.values_of_t::<usize>("flag-antennas") {
            Ok(antenna_idxs) => {
                for (value_idx, &antenna_idx) in antenna_idxs.iter().enumerate() {
                    if antenna_idx >= *num_ants {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-antennas <ANTS>...".into(),
                            expected: format!("antenna_idx < num_ants={}", *num_ants),
                            received: format!(
                                "antenna_idxs[{value_idx}]={antenna_idx}. all:{antenna_idxs:?}"
                            ),
                        }));
                    }
                    flag_ctx.antenna_flags[antenna_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        if matches.is_present("flag-autos") {
            flag_ctx.autos = true;
        }
        if matches.is_present("flag-dc") {
            flag_ctx.flag_dc = true;
        }
        if matches.is_present("no-flag-dc") {
            flag_ctx.flag_dc = false;
        }
        match matches.value_of_t::<usize>("flag-edge-chans") {
            Ok(n) => {
                if n >= flag_ctx.fine_chan_flags.len() / 2 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--flag-edge-chans <COUNT>".into(),
                        expected: "fewer than N/2-1 fine channels".into(),
                        received: format!("{n}"),
                    }));
                }
                Self::flag_edge_channels(n, &mut flag_ctx.fine_chan_flags);
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.value_of_t::<usize>("flag-edge-width") {
            Ok(width) => {
                let fine_chan_width = *corr_fine_chan_width_hz / 1000;
                let n = width as f32 / fine_chan_width as f32;
                if (n - n.floor()).abs() > 0.00001 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--flag-edge-width <COUNT>".into(),
                        expected: format!("multiple of fine channel width ({fine_chan_width})"),
                        received: format!("{width}"),
                    }));
                }
                if n as usize >= flag_ctx.fine_chan_flags.len() / 2 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--flag-edge-width <COUNT>".into(),
                        expected: "width equal to fewer than N/2-1 fine channels".into(),
                        received: format!("{n}"),
                    }));
                }
                Self::flag_edge_channels(n as usize, &mut flag_ctx.fine_chan_flags);
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.value_of_t::<f32>("flag-init") {
            Ok(init_time) => {
                let d = *corr_int_time_ms as f32 / 1000.0;
                if init_time % d < 0.000001 {
                    flag_ctx.flag_init = init_time;
                } else {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "flag-init".into(),
                        expected: format!("A multiple of the timestep length ({d})"),
                        received: format!("{init_time}"),
                    }));
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.value_of_t::<f32>("flag-end") {
            Ok(end_time) => {
                let d = *corr_int_time_ms as f32 / 1000.0;
                if end_time % d < 0.000001 {
                    flag_ctx.flag_end = end_time;
                } else {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "flag-end".into(),
                        expected: format!("A multiple of the timestep length ({d})"),
                        received: format!("{end_time}"),
                    }));
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.value_of_t::<u32>("flag-init-steps") {
            Ok(init_steps) => {
                flag_ctx.flag_init = init_steps as f32 * *corr_int_time_ms as f32 / 1000.0;
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        match matches.value_of_t::<u32>("flag-end-steps") {
            Ok(end_steps) => {
                flag_ctx.flag_end = end_steps as f32 * *corr_int_time_ms as f32 / 1000.0;
            }
            Err(err) => match err.kind() {
                ArgumentNotFound => {}
                _ => return Err(err.into()),
            },
        };
        flag_ctx.finalise_flag_settings(corr_ctx);
        Ok(flag_ctx)
    }

    fn flag_edge_channels(n: usize, channels: &mut [bool]) {
        channels.iter_mut().take(n).for_each(|x| {
            *x = true;
        });
        channels.iter_mut().rev().take(n).for_each(|x| {
            *x = true;
        });
    }

    fn parse_avg_matches(
        matches: &clap::ArgMatches,
        corr_ctx: &CorrelatorContext,
    ) -> Result<(usize, usize), BirliError> {
        let CorrelatorContext {
            metafits_context: meta_ctx,
            ..
        } = corr_ctx;
        let MetafitsContext {
            corr_int_time_ms,
            corr_fine_chan_width_hz,
            ..
        } = meta_ctx;

        let avg_time: usize = match (
            matches.value_of_t::<usize>("avg-time-factor"),
            matches.value_of_t::<f64>("avg-time-res"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) | (_, Err(err)) if err.kind() != ArgumentNotFound => {
                return Err(err.into())
            }
            (Ok(_), Ok(_)) => {
                unreachable!("--avg-time-res conflicts with --avg-time-factor, enforced by clap")
            }
            (Ok(factor), _) => {
                if factor == 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-time-factor <FACTOR>".into(),
                        expected: "a positive, non-zero integer".into(),
                        received: format!("{factor}"),
                    }));
                }
                factor
            }
            (_, Ok(res)) => {
                let int_time_s = *corr_int_time_ms as f64 / 1e3;
                let ratio = res / int_time_s;
                if ratio.is_infinite() || ratio.fract() > 1e-6 || ratio < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-time-res <RES>".into(),
                        expected: format!("a multiple of the integration time, {int_time_s} [s]"),
                        received: format!("{res}"),
                    }));
                }
                ratio.round() as _
            }
            _ => 1,
        };
        let avg_freq: usize = match (
            matches.value_of_t::<usize>("avg-freq-factor"),
            matches.value_of_t::<f64>("avg-freq-res"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) | (_, Err(err)) if err.kind() != ArgumentNotFound => {
                return Err(err.into())
            }
            (Ok(_), Ok(_)) => {
                unreachable!("--avg-freq-res conflicts with --avg-freq-factor, enforced by clap")
            }
            (Ok(factor), _) => {
                if factor == 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-freq-factor <FACTOR>".into(),
                        expected: "a positive, non-zero integer".into(),
                        received: format!("{factor}"),
                    }));
                }
                factor
            }
            (_, Ok(res)) => {
                let fine_chan_width_khz = *corr_fine_chan_width_hz as f64 / 1e3;
                let ratio = res / fine_chan_width_khz;
                if ratio.is_infinite() || ratio.fract() > 1e-6 || ratio < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-freq-res <RES>".into(),
                        expected: format!(
                            "a multiple of the fine channel width, {fine_chan_width_khz} [KHz]"
                        ),
                        received: format!("{res}"),
                    }));
                }
                ratio.round() as _
            }
            _ => 1,
        };
        Ok((avg_time, avg_freq))
    }

    fn parse_chunk_matches(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
        avg_time: usize,
        vis_sel: &VisSelection,
    ) -> Result<Option<usize>, BirliError> {
        let CorrelatorContext {
            metafits_context: meta_ctx,
            ..
        } = corr_ctx;
        let MetafitsContext {
            num_corr_fine_chans_per_coarse: fine_chans_per_coarse,
            ..
        } = meta_ctx;
        let num_timesteps_per_chunk: Option<usize> = match (
            matches.value_of_t::<usize>("time-chunk"),
            matches.value_of_t::<f64>("max-memory"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) | (_, Err(err)) if err.kind() != ArgumentNotFound => {
                return Err(err.into())
            }
            (Ok(_), Ok(_)) => {
                unreachable!("--time-chunk conflicts with --max-memory, enforced by clap")
            }
            (Ok(steps), _) => {
                if steps % avg_time != 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--time-chunk <STEPS>".into(),
                        expected: format!(
                            "a multiple of the temporal averaging factor, {avg_time}"
                        ),
                        received: format!("{steps}"),
                    }));
                }
                Some(steps)
            }
            (_, Ok(mem_gib)) => {
                let max_mem_bytes = mem_gib * 1024.0_f64.powi(3);
                if max_mem_bytes < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--max-memory <GIBIBYTES>".into(),
                        expected: "at least one Byte".into(),
                        received: format!("{max_mem_bytes}B"),
                    }));
                }
                let bytes_selected = vis_sel.estimate_bytes_best(*fine_chans_per_coarse);
                let bytes_per_timestep = bytes_selected / vis_sel.timestep_range.len();
                let bytes_per_avg_time = bytes_per_timestep * avg_time;
                if max_mem_bytes < bytes_selected as f64 {
                    if max_mem_bytes < bytes_per_avg_time as f64 {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--max-memory <GIBIBYTES>".into(),
                            expected: format!("at least enough memory for an averaged timestep ({} * {:.02} = {:.02} GiB)", avg_time, bytes_per_timestep as f64 / 1024.0_f64.powi(3), bytes_per_avg_time as f64 / 1024.0_f64.powi(3)),
                            received: format!("{}GiB", max_mem_bytes / 1024.0_f64.powi(3)),
                        }));
                    }
                    Some((max_mem_bytes / bytes_per_avg_time as f64).floor() as usize * avg_time)
                } else {
                    None
                }
            }
            _ => None,
        };

        // validate chunk size
        if let Some(chunk_size) = num_timesteps_per_chunk {
            info!("chunking output to {} timesteps per chunk", chunk_size);
        }

        Ok(num_timesteps_per_chunk)
    }

    fn parse_prep_matches(
        matches: &clap::ArgMatches,
        corr_ctx: &CorrelatorContext,
    ) -> Result<PreprocessContext<'a>, BirliError> {
        let mut prep_ctx = PreprocessContext {
            draw_progress: !matches.is_present("no-draw-progress"),
            ..PreprocessContext::default()
        };
        let CorrelatorContext {
            metafits_context: meta_ctx,
            mwa_version,
            ..
        } = corr_ctx;
        let MetafitsContext {
            oversampled,
            cable_delays_applied,
            geometric_delays_applied,
            deripple_applied,
            ..
        } = meta_ctx;
        prep_ctx.array_pos = if matches.is_present("emulate-cotter") {
            info!("Using array position from Cotter.");
            LatLngHeight {
                longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
                latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
                height_metres: COTTER_MWA_HEIGHT_METRES,
            }
        } else {
            info!("Using default MWA array position.");
            LatLngHeight::mwa()
        };
        prep_ctx.phase_centre = match (
            matches
                .values_of_t::<f64>("phase-centre")
                .map(|v| (v[0], v[1])),
            matches.is_present("pointing-centre"),
        ) {
            (Err(err), _) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (Ok(_), true) => {
                unreachable!("--phase-centre conflicts with --pointing-centre, enforced by clap");
            }
            (Ok((ra, dec)), _) => RADec::from_degrees(ra, dec),
            (_, true) => RADec::from_mwalib_tile_pointing(meta_ctx),
            _ => RADec::from_mwalib_phase_or_pointing(meta_ctx),
        };
        prep_ctx.correct_van_vleck = match (matches.is_present("van-vleck"), mwa_version) {
            (true, MWAVersion::CorrLegacy) => true,
            (true, _) => {
                return Err(BirliError::CLIError(InvalidCommandLineArgument {
                    option: "--van-vleck".into(),
                    expected: "legacy correlator files".into(),
                    received: mwa_version.to_string(),
                }))
            }
            _ => false,
        };
        prep_ctx.correct_cable_lengths = {
            let cable_delays_disabled = matches.is_present("no-cable-delay");
            info!(
                "cable corrections: applied={:?}, disabled={}",
                cable_delays_applied, cable_delays_disabled
            );
            matches!(
                cable_delays_applied,
                CableDelaysApplied::NoCableDelaysApplied
            ) && !cable_delays_disabled
        };
        prep_ctx.correct_digital_gains = !matches.is_present("no-digital-gains");
        prep_ctx.passband_gains = match matches.value_of("passband-gains") {
            None | Some("none") => None,
            Some(g) if g == "jake" => {
                info!("passband gains: {} (mwax, not oversampled)", g);
                Some(PFB_JAKE_2022_200HZ)
            }
            Some(g) if g == "jake_oversampled" => {
                info!("passband gains: {} (mwax, oversampled)", g);
                Some(OSPFB_JAKE_2025_200HZ)
            }
            Some(g) if g == "cotter" => {
                info!("passband gains: {} (legacy)", g);
                Some(PFB_COTTER_2014_10KHZ)
            }
            Some(g) if g == "auto" => {
                if *deripple_applied {
                    info!("passband gains: {} (disabled, deripple already applied)", g);
                    None
                } else {
                    match (mwa_version, oversampled) {
                        (MWAVersion::CorrMWAXv2, false) => {
                            info!("passband gains: {} (mwax, not oversampled)", g);
                            Some(PFB_JAKE_2022_200HZ)
                        }
                        (MWAVersion::CorrMWAXv2, true) => {
                            info!("passband gains: {} (mwax, oversampled)", g);
                            Some(OSPFB_JAKE_2025_200HZ)
                        }
                        (MWAVersion::CorrLegacy | MWAVersion::CorrOldLegacy, _) => {
                            info!("passband gains: {} (legacy)", g);
                            Some(PFB_COTTER_2014_10KHZ)
                        }
                        (ver, _) => {
                            return Err(BadMWAVersion {
                                message: "unknown mwa version".into(),
                                version: ver.to_string(),
                            })
                        }
                    }
                }
            }
            Some(option) => panic!("unknown option for --passband-gains: {option}"),
        };
        prep_ctx.correct_geometry = {
            let geometric_delays_disabled = matches.is_present("no-geometric-delay");
            info!(
                "geometric corrections: applied={:?}, disabled={}",
                geometric_delays_applied, !geometric_delays_disabled
            );
            matches!(geometric_delays_applied, GeometricDelaysApplied::No)
                && !geometric_delays_disabled
        };
        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                prep_ctx.aoflagger_strategy = if matches.is_present("no-rfi") {
                    None
                } else {
                    match matches.value_of_t("aoflagger-strategy") {
                        Err(err) if err.kind() != ArgumentNotFound => return Err(err.into()),
                        Ok(strategy) => Some(strategy),
                        Err(_) => Some(unsafe {
                            cxx_aoflagger_new().FindStrategyFileMWA()
                        }),
                    }
                };
            }
        }
        Ok(prep_ctx)
    }

    /// Parse an iterator of arguments, `args` into a `BirliContext`.
    ///
    /// # Errors
    ///
    /// Can raise:
    /// - `clap::Error` if clap cannot parse `args`
    /// - `mwalib::MwalibError` if mwalib can't open the input files.
    /// - `BirliError::CLIError` if the arguments are invalid.
    pub fn from_args<I, T>(args: I) -> Result<Self, BirliError>
    where
        I: IntoIterator<Item = T> + Debug,
        T: Into<OsString> + Clone,
    {
        debug!("args:\n{:?}", &args);

        let matches = Self::get_matches(args)?;
        trace!("arg matches:\n{:?}", &matches);

        let io_ctx = Self::parse_io_matches(&matches);
        let corr_ctx = io_ctx.get_corr_ctx()?;
        debug!("mwalib correlator context:\n{}", &corr_ctx);
        let flag_ctx = Self::parse_flag_matches(&corr_ctx, &matches)?;
        let vis_sel = Self::parse_vis_sel_matches(&corr_ctx, &flag_ctx, &matches)?;
        let prep_ctx = Self::parse_prep_matches(&matches, &corr_ctx)?;
        let (avg_time, avg_freq) = Self::parse_avg_matches(&matches, &corr_ctx)?;
        let num_timesteps_per_chunk =
            Self::parse_chunk_matches(&corr_ctx, &matches, avg_time, &vis_sel)?;
        let channel_range_sel = Self::parse_sel_chan_ranges(&corr_ctx, &matches)?;
        let result = Self {
            corr_ctx,
            prep_ctx,
            vis_sel,
            flag_ctx,
            io_ctx,
            avg_time,
            avg_freq,
            num_timesteps_per_chunk,
            channel_range_sel,
        };

        info!("{}", &result);

        if matches.is_present("dry-run") {
            return Err(DryRun {});
        }

        Ok(result)
    }

    /// Call `run()` for every channel range in `self.channel_range_sel`.
    ///
    /// # Errors
    /// see: `run()`
    pub fn run_ranges(self) -> Result<(), BirliError> {
        let ranges = self.channel_range_sel.ranges.clone();
        let coarse_chans = self.corr_ctx.coarse_chans.clone();
        let original_io_ctx = self.io_ctx.clone();
        let mut ranged_context: BirliContext = BirliContext {
            corr_ctx: self.corr_ctx,
            prep_ctx: self.prep_ctx,
            vis_sel: self.vis_sel,
            flag_ctx: self.flag_ctx,
            io_ctx: self.io_ctx,
            avg_time: self.avg_time,
            avg_freq: self.avg_freq,
            num_timesteps_per_chunk: self.num_timesteps_per_chunk,
            channel_range_sel: self.channel_range_sel,
        };
        for &(range_start, range_end) in &ranges {
            ranged_context.vis_sel.coarse_chan_range = range_start..range_end + 1;

            ranged_context.io_ctx = original_io_ctx.clone();
            let suffix = match (range_start, range_end) {
                // single channel range
                (m, n) if m == n => format!("_ch{}", coarse_chans[range_start].rec_chan_number),
                // multi channel range
                (m, n) if m < n => format!(
                    "_ch{}-{}",
                    coarse_chans[range_start].rec_chan_number,
                    coarse_chans[range_end].rec_chan_number
                ),
                _ => unreachable!(),
            };

            if let Some(path) = ranged_context.io_ctx.ms_out.as_mut() {
                path.set_file_name(format!(
                    "{}{}.{}",
                    path.file_stem().unwrap().to_str().unwrap(),
                    suffix,
                    path.extension().unwrap().to_str().unwrap()
                ));
            }
            if let Some(path) = ranged_context.io_ctx.uvfits_out.as_mut() {
                path.set_file_name(format!(
                    "{}{}.{}",
                    path.file_stem().unwrap().to_str().unwrap(),
                    suffix,
                    path.extension().unwrap().to_str().unwrap()
                ));
            }
            ranged_context.run()?;
        }
        Ok(())
    }

    /// Read, Preprocess and write corrected visibilities chunks.
    ///
    /// # Errors
    ///
    /// can raise:
    /// - `BadArrayShape` if the shape of the calibration solutions
    ///   is incompatible with the visibility shape.
    /// - preprocessing errors
    pub fn run(&self) -> Result<(), BirliError> {
        let Self {
            corr_ctx,
            vis_sel,
            flag_ctx,
            io_ctx,
            avg_time,
            avg_freq,
            num_timesteps_per_chunk,
            ..
        } = self;
        let mut prep_ctx = self.prep_ctx.clone();

        // ////////// //
        // Prepare IO //
        // ////////// //

        let vis_ctx = VisContext::from_mwalib(
            corr_ctx,
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            &vis_sel.baseline_idxs,
            *avg_time,
            *avg_freq,
        );

        let num_avg_timesteps = vis_ctx.num_avg_timesteps();
        let draw_target = if prep_ctx.draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        };
        let write_progress =
            indicatif::ProgressBar::with_draw_target(Some(num_avg_timesteps as u64), draw_target);
        write_progress.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .unwrap()
                .progress_chars("=> "),
        );
        write_progress.set_message("write vis");

        // TODO: move phase_centre, array_pos out of prep_ctx
        let obs_ctx = ObsContext {
            phase_centre: prep_ctx.phase_centre,
            array_pos: prep_ctx.array_pos,
            ..ObsContext::from_mwalib(&corr_ctx.metafits_context)
        };

        let mwa_ctx = MwaObsContext::from_mwalib(&corr_ctx.metafits_context);

        prep_ctx.calsols = if let Some(ref calsol_file) = io_ctx.aocalsols_in {
            let calsols = AOCalSols::read_andre_binary(calsol_file).unwrap();
            assert!(
                calsols.di_jones.dim().0 == 1,
                "only 1 timeblock must be supplied for calsols. \
                Instead found {} timeblocks. dimensions {:?}",
                calsols.di_jones.dim().1,
                calsols.di_jones.dim()
            );
            let calsol_chans = calsols.di_jones.dim().2;
            if calsol_chans % corr_ctx.num_coarse_chans != 0 {
                return Err(BirliError::BadArrayShape(BadArrayShape {
                    argument: "input AO calibration solutions",
                    function: "BirliContext::run",
                    expected: format!(
                        "a multiple of metafits_num_coarse_chans={}",
                        corr_ctx.metafits_context.num_metafits_coarse_chans
                    ),
                    received: format!("{calsol_chans}"),
                }));
            }
            let num_calsol_fine_chans_per_coarse = calsol_chans / corr_ctx.num_coarse_chans;
            Some(
                calsols
                    .di_jones
                    .index_axis(Axis(0), 0)
                    .slice(s![
                        ..,
                        (vis_sel.coarse_chan_range.start * num_calsol_fine_chans_per_coarse)
                            ..(vis_sel.coarse_chan_range.end * num_calsol_fine_chans_per_coarse)
                    ])
                    .to_owned(),
            )
        } else {
            None
        };

        let args_strings = std::env::args().collect_vec();
        let cmd_line = shlex::try_join(args_strings.iter().map(String::as_str)).unwrap();
        let application = format!("{PKG_NAME} {PKG_VERSION}");
        let message = prep_ctx.as_comment();
        let history = History {
            cmd_line: Some(&cmd_line),
            application: Some(&application),
            message: Some(&message),
        };
        let (s_lat, c_lat) = obs_ctx.array_pos.latitude_rad.sin_cos();
        let (antenna_names, antenna_positions): (Vec<String>, Vec<marlu::XyzGeodetic>) = corr_ctx
            .metafits_context
            .antennas
            .iter()
            .map(|a| {
                let enh = ENH {
                    e: a.east_m,
                    n: a.north_m,
                    h: a.height_m,
                };
                let xyz = enh.to_xyz_inner(s_lat, c_lat);
                (a.tile_name.clone(), xyz)
            })
            .unzip();
        let dut1 = hifitime::Duration::from_seconds(corr_ctx.metafits_context.dut1.unwrap_or(0.0));
        let mut uvfits_writer = io_ctx.uvfits_out.as_ref().map(|uvfits_out| {
            with_increment_duration!("init", {
                UvfitsWriter::from_marlu(
                    uvfits_out,
                    &vis_ctx,
                    obs_ctx.array_pos,
                    obs_ctx.phase_centre,
                    dut1,
                    obs_ctx.name.as_deref(),
                    antenna_names,
                    antenna_positions.clone(),
                    true,
                    Some(&history),
                )
                .expect("unable to initialize uvfits writer")
            })
        });
        let mut ms_writer = io_ctx.ms_out.as_ref().map(|ms_out| {
            let writer = MeasurementSetWriter::new(
                ms_out,
                obs_ctx.phase_centre,
                obs_ctx.array_pos,
                antenna_positions,
                dut1,
                true,
            );
            println!(
                "Writing to MS: {} with {} chans selected",
                ms_out.display(),
                vis_ctx.num_sel_chans
            );
            with_increment_duration!("init", {
                writer
                    .initialize_mwa(
                        &vis_ctx,
                        &obs_ctx,
                        &mwa_ctx,
                        Some(&history),
                        &vis_sel.coarse_chan_range,
                    )
                    .expect("unable to initialize ms writer");
            });
            writer
        });

        #[cfg(feature = "aoflagger")]
        let (aoflagger_version, aoflagger_strategy) = {
            let mut major = 0;
            let mut minor = 0;
            let mut subminor = 0;
            unsafe {
                aoflagger_sys::cxx_aoflagger_new().GetVersion(
                    &mut major,
                    &mut minor,
                    &mut subminor,
                );
            }
            (
                Some(format!("v{major}.{minor}.{subminor}")),
                prep_ctx.aoflagger_strategy.clone(),
            )
        };
        #[cfg(not(feature = "aoflagger"))]
        let (aoflagger_version, aoflagger_strategy) = (None, None);

        let mut flag_file_set = io_ctx.flag_template.as_ref().map(|flag_template| {
            FlagFileSet::new(
                flag_template,
                corr_ctx,
                vis_sel,
                aoflagger_version,
                aoflagger_strategy,
            )
            .expect("cannot create flag file writer")
        });

        // //////// //
        // Chunking //
        // //////// //

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let chunk_size = num_timesteps_per_chunk.map_or_else(
            || {
                let num_timesteps: usize = vis_sel.timestep_range.len();
                num_timesteps
            },
            |steps| steps,
        );

        // Allocate our big arrays once, reuse them for each chunk unless the chunk shape changes
        let chunk_vis_sel = VisSelection {
            timestep_range: (vis_sel.timestep_range.start
                ..vis_sel.timestep_range.start + chunk_size),
            ..vis_sel.clone()
        };
        let mut jones_array = chunk_vis_sel.allocate_jones(fine_chans_per_coarse)?;
        let mut flag_array = chunk_vis_sel.allocate_flags(fine_chans_per_coarse)?;
        let mut weight_array = chunk_vis_sel.allocate_weights(fine_chans_per_coarse)?;

        for mut timestep_chunk in &vis_sel.timestep_range.clone().chunks(chunk_size) {
            let chunk_first_timestep = timestep_chunk.next().expect("zero-sized chunk");
            let chunk_vis_sel = VisSelection {
                timestep_range: (chunk_first_timestep
                    ..(timestep_chunk.last().unwrap_or(chunk_first_timestep) + 1)),
                ..vis_sel.clone()
            };
            if num_timesteps_per_chunk.is_some() {
                info!(
                    "processing timestep chunk {:?} of {:?} % {}",
                    chunk_vis_sel.timestep_range,
                    vis_sel.timestep_range.clone(),
                    chunk_size
                );
            }

            // only reallocate arrays if the chunk dimensions have changed.
            let chunk_dims = chunk_vis_sel.get_shape(fine_chans_per_coarse);
            let (mut jones_array, mut flag_array, mut weight_array) = if jones_array.dim()
                == chunk_dims
            {
                (
                    jones_array.view_mut(),
                    flag_array.view_mut(),
                    weight_array.view_mut(),
                )
            } else {
                (
                    jones_array.slice_mut(s![0..chunk_dims.0, 0..chunk_dims.1, 0..chunk_dims.2]),
                    flag_array.slice_mut(s![0..chunk_dims.0, 0..chunk_dims.1, 0..chunk_dims.2]),
                    weight_array.slice_mut(s![0..chunk_dims.0, 0..chunk_dims.1, 0..chunk_dims.2]),
                )
            };

            // populate flags
            flag_ctx.set_flags(
                flag_array.view_mut(),
                &chunk_vis_sel.timestep_range,
                &chunk_vis_sel.coarse_chan_range,
                &chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )?;

            // populate visibilities
            with_increment_duration!(
                "read",
                read_mwalib(
                    &chunk_vis_sel,
                    corr_ctx,
                    jones_array.view_mut(),
                    flag_array.view_mut(),
                    prep_ctx.draw_progress
                )?
            );

            // populate weights
            weight_array.fill(vis_ctx.weight_factor() as f32);

            prep_ctx.preprocess(
                corr_ctx,
                jones_array.view_mut(),
                weight_array.view_mut(),
                flag_array.view_mut(),
                &chunk_vis_sel,
            )?;

            let mut ssins = SSINS::new(jones_array.view(), &corr_ctx, &chunk_vis_sel);
            // this requires the aoflagger feature flag.
            #[cfg(feature = "aoflagger")]
            {
                ssins.flag(
                    None,
                    // Some( "/Users/dev/Code/aoflagger/data/strategies/generic-minimal.lua".to_string(), ),
                );
                println!("SSINS flags:");
                let flag_dims = ssins.flag_array.dim();
                for t in 0..flag_dims.0 {
                    for f in 0..flag_dims.1 {
                        let c = if ssins.flag_array[[t, f]] { "#" } else { " " };
                        print!("{}", c);
                    }
                    println!();
                }
            }

            let eavils = EAVILS::new(jones_array.view(), &corr_ctx, &chunk_vis_sel);

            // Save SSINS metrics to CSV if output path is specified
            if let Some(metrics_path) = &io_ctx.metrics_out {
                let path = Path::new(metrics_path.to_str().unwrap());
                if path.exists() {
                    std::fs::remove_file(path).unwrap();
                }
                let mut fptr = FitsFile::create(path).open().unwrap();
                if let Err(e) = ssins.save_to_fits(&mut fptr) {
                    log::warn!("Failed to save SSINS metrics to FITS: {}", e);
                }
                if let Err(e) = eavils.save_to_fits(&mut fptr) {
                    log::warn!("Failed to save EAVILS metrics to FITS: {}", e);
                }
                log::info!(
                    "Saved SSINS and EAVILS metrics to {}",
                    metrics_path.display()
                );
            }

            if let Some(flag_file_set) = flag_file_set.as_mut() {
                with_increment_duration!(
                    "write",
                    flag_file_set
                        .write_flag_array(flag_array.view(), prep_ctx.draw_progress)
                        .expect("unable to write flags")
                );
            }

            // bake flags into weights
            for (weight, flag) in izip!(weight_array.iter_mut(), flag_array.iter()) {
                *weight = if *flag {
                    -(*weight).abs()
                } else {
                    (*weight).abs()
                };
            }

            let sel_timestep_chunks = chunk_vis_sel.timestep_range.clone().chunks(*avg_time);
            izip!(
                &sel_timestep_chunks,
                jones_array.axis_chunks_iter(Axis(0), *avg_time),
                weight_array.axis_chunks_iter(Axis(0), *avg_time),
            )
            .for_each(|(mut out_times, jones_array, weight_array)| {
                let first_avg_timestep = out_times.next().expect("zero-sized chunk");
                let last_avg_timestep = out_times.last().unwrap_or(first_avg_timestep) + 1;
                let avg_timestep_range = first_avg_timestep..last_avg_timestep;
                let avg_vis_sel = VisSelection {
                    timestep_range: avg_timestep_range,
                    ..chunk_vis_sel.clone()
                };

                let avg_vis_ctx = VisContext::from_mwalib(
                    corr_ctx,
                    &avg_vis_sel.timestep_range,
                    &avg_vis_sel.coarse_chan_range,
                    &avg_vis_sel.baseline_idxs,
                    *avg_time,
                    *avg_freq,
                );

                // output uvfits
                if let Some(uvfits_writer) = uvfits_writer.as_mut() {
                    with_increment_duration!(
                        "write",
                        uvfits_writer
                            .write_vis(jones_array.view(), weight_array.view(), &avg_vis_ctx)
                            .expect("unable to write uvfits")
                    );
                }

                // output ms
                if let Some(ms_writer) = ms_writer.as_mut() {
                    with_increment_duration!(
                        "write",
                        ms_writer
                            .write_vis(jones_array.view(), weight_array.view(), &avg_vis_ctx)
                            .expect("unable to write ms")
                    );
                }

                write_progress.inc(1);
            });
        }
        // Finalise the uvfits writer.
        if let Some(uvfits_writer) = uvfits_writer.as_mut() {
            with_increment_duration!(
                "write",
                uvfits_writer
                    .finalise()
                    .expect("couldn't write antenna table to uvfits")
            );
        };

        // Finalise the MS writer.
        if let Some(ms_writer) = ms_writer.as_mut() {
            with_increment_duration!("write", ms_writer.finalise().expect("couldn't finalise MS"));
        };

        write_progress.finish();

        // Finalise the mwaf files.
        if let Some(flag_file_set) = flag_file_set {
            flag_file_set
                .finalise()
                .expect("couldn't finalise mwaf files");
        }

        Ok(())
    }
}

#[cfg(test)]
mod argparse_tests {
    use marlu::{mwalib::MWAVersion, RADec};

    use approx::{assert_abs_diff_eq, assert_abs_diff_ne};

    use crate::{
        passband_gains::{OSPFB_JAKE_2025_200HZ, PFB_COTTER_2014_10KHZ, PFB_JAKE_2022_200HZ},
        test_common::{get_1254670392_avg_paths, get_mwax_data_paths},
        BirliContext, BirliError, VisSelection,
    };
    use tempfile::tempdir;

    /// Middle channel rounded-down is DC flagged when `flag_dc` is set.
    #[test]
    fn test_dc_flagging_odd_channel_count() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--flag-init", "0",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext {
            flag_ctx,
            vis_sel,
            num_timesteps_per_chunk,
            corr_ctx,
            ..
        } = BirliContext::from_args(&args).unwrap();

        let chunk_size = if let Some(steps) = num_timesteps_per_chunk {
            steps
        } else {
            vis_sel.timestep_range.len()
        };
        let chunk_vis_sel = VisSelection {
            timestep_range: (vis_sel.timestep_range.start
                ..vis_sel.timestep_range.start + chunk_size),
            ..vis_sel
        };
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = chunk_vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        flag_ctx
            .set_flags(
                flag_array.view_mut(),
                &chunk_vis_sel.timestep_range,
                &chunk_vis_sel.coarse_chan_range,
                &chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )
            .unwrap();

        assert!(flag_array[[0, fine_chans_per_coarse / 2, 0]]);
        assert!(!flag_array[[0, fine_chans_per_coarse / 2 + 1, 0]]);
    }

    /// DC flagging is on by default for old correlator inputs.
    #[test]
    fn test_automatic_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.flag_dc);
    }

    /// pfb gains is cotter by default for legacy correlator.
    #[test]
    fn test_auto_pfb_legacy() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { prep_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(prep_ctx.passband_gains.is_some());
        assert_abs_diff_eq!(
            prep_ctx.passband_gains.unwrap()[0],
            PFB_COTTER_2014_10KHZ[0]
        );
        // santiy check
        assert_abs_diff_ne!(prep_ctx.passband_gains.unwrap()[0], PFB_JAKE_2022_200HZ[0]);
    }

    /// invalid pfb gains are rejected.
    #[test]
    fn test_invalid_pfb() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "invalid",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let birli_ctx = BirliContext::from_args(&args);

        assert!(birli_ctx.is_err());
    }

    /// correlator version throws an error when pfb gains are auto,
    /// selects oversampled pfb for oversampled ovservations.
    #[test]
    fn test_corr_type_auto_pfb() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "auto",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        // create an initial birli context so we can give it an unfamiliar correlator version and parse again
        let mut birli_ctx = BirliContext::from_args(&args).unwrap();
        let matches = BirliContext::get_matches(&args).unwrap();

        // case 1: invalid corr type
        birli_ctx.corr_ctx.mwa_version = MWAVersion::VCSLegacyRecombined;
        let result = BirliContext::<'_>::parse_prep_matches(&matches, &birli_ctx.corr_ctx);
        assert!(matches!(result, Err(BirliError::BadMWAVersion { .. })));

        // case 2: mwax, oversampled
        // TODO: make this its own test case, need new metafits, ideally with sigchain corrections
        birli_ctx.corr_ctx.mwa_version = MWAVersion::CorrMWAXv2;
        birli_ctx.corr_ctx.metafits_context.oversampled = true;
        let result = BirliContext::<'_>::parse_prep_matches(&matches, &birli_ctx.corr_ctx);
        assert!(result.is_ok());
        assert_abs_diff_eq!(
            result.unwrap().passband_gains.unwrap()[0..10],
            OSPFB_JAKE_2025_200HZ[0..10]
        );

        // case 3: mwax, not oversampled
        birli_ctx.corr_ctx.metafits_context.oversampled = false;
        let result = BirliContext::<'_>::parse_prep_matches(&matches, &birli_ctx.corr_ctx);
        assert!(result.is_ok());
        assert_abs_diff_eq!(
            result.unwrap().passband_gains.unwrap()[0..10],
            PFB_JAKE_2022_200HZ[0..10]
        );
    }

    /// pfb gains is jake by default for mwax correlator.
    #[test]
    fn test_auto_pfb_mwax() {
        let (metafits_path, gpufits_paths) = get_mwax_data_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { prep_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(prep_ctx.passband_gains.is_some());
        assert_abs_diff_eq!(prep_ctx.passband_gains.unwrap()[0], PFB_JAKE_2022_200HZ[0]);
        // santiy check
        assert_abs_diff_ne!(
            prep_ctx.passband_gains.unwrap()[0],
            PFB_COTTER_2014_10KHZ[0]
        );
    }

    /// pfb gains is disabled when `deripple_applied` is true
    #[test]
    fn test_no_pfb_when_deripple_applied() {
        #[rustfmt::skip]
        let args_no_dry_run = vec![
            "birli",
            "-m", "tests/data/1439922144_deripple/1439922144.metafits",
            "--no-draw-progress",
            "tests/data/1439922144_deripple/1439922144_20250822182206_ch137_000.fits",
        ];

        let BirliContext { prep_ctx, .. } = BirliContext::from_args(&args_no_dry_run).unwrap();

        assert!(prep_ctx.passband_gains.is_none());
    }

    /// DC flagging is off by default for MWAX inputs.
    #[test]
    fn test_no_automatic_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_mwax_data_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(!flag_ctx.flag_dc);
    }

    /// DC flagging for old correlator inputs is disabled by `no-flag-dc`.
    #[test]
    fn test_suppress_automatic_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-flag-dc",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(!flag_ctx.flag_dc);
    }

    /// DC flagging for MWAX inputs is enabled by `flag-dc`.
    #[test]
    fn test_force_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_mwax_data_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--flag-dc",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.flag_dc);
    }

    /// Implicit DC flagging for old correlator inputs can be made explicit with `flag-dc`.
    #[test]
    fn test_explicit_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--flag-dc",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.flag_dc);
    }

    /// DC flagging for MWAX can be explicitly disabled with `flag-dc`.
    #[test]
    fn test_explicit_disable_dc_flagging() {
        let (metafits_path, gpufits_paths) = get_mwax_data_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-flag-dc",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(!flag_ctx.flag_dc);
    }

    /// Flag 3 fine channels on the edge of a coarse channel with `flag-edge-chans`
    #[test]
    fn test_flag_edge_chans_ms_output() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("test_flags.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--flag-edge-chans", "3",
            "--flag-init", "0",
            "--no-draw-progress",
            gpufits_paths[0],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        // Verify edge channels are flagged in the flag context
        assert!(birli_ctx.flag_ctx.fine_chan_flags[0]); // first channel
        assert!(birli_ctx.flag_ctx.fine_chan_flags[1]); // second channel
        assert!(birli_ctx.flag_ctx.fine_chan_flags[2]); // third channel
        assert!(!birli_ctx.flag_ctx.fine_chan_flags[3]); // fourth should not be flagged
        assert!(birli_ctx.flag_ctx.fine_chan_flags[31]); // last channel
        assert!(birli_ctx.flag_ctx.fine_chan_flags[30]); // second to last channel
        assert!(birli_ctx.flag_ctx.fine_chan_flags[29]); // third to last channel
        assert!(!birli_ctx.flag_ctx.fine_chan_flags[28]); // fourth from last should not be flagged

        // Run birli to create the measurement set
        birli_ctx.run().unwrap();

        // Check that the measurement set was created
        assert!(ms_path.exists());

        // Test that FLAGS are properly written by reading the measurement set
        use marlu::rubbl_casatables::{Table, TableOpenMode};
        let mut main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        let num_rows = main_table.n_rows();
        assert!(num_rows > 0);

        // Check the FLAG column exists and has the correct flags set
        // Get flags for the first row
        let obs_flags = main_table.get_cell_as_vec::<bool>("FLAG", 0).unwrap();

        // The flags should reflect our edge channel flagging
        // With 32 channels per coarse channel and 4 pols, we expect:
        // Channels 0,1,2 and 29,30,31 to be flagged for each pol
        let num_chans = obs_flags.len() / 4; // 4 polarizations
        assert_eq!(num_chans, 32);

        // Check first few channels are flagged for first pol
        assert!(obs_flags[0]); // pol 0, chan 0
        assert!(obs_flags[4]); // pol 0, chan 1
        assert!(obs_flags[8]); // pol 0, chan 2
        assert!(!obs_flags[12]); // pol 0, chan 3 (should not be flagged)

        // Check last few channels are flagged for first pol
        assert!(obs_flags[116]); // pol 0, chan 29 (29*4 = 116)
        assert!(obs_flags[120]); // pol 0, chan 30
        assert!(obs_flags[124]); // pol 0, chan 31
    }

    /// Flag 120kHz on the edges of a coarse channel with `flag-edge-width`
    #[test]
    fn test_flag_edge_width() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--flag-edge-width", "120",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.fine_chan_flags[0]);
        assert!(flag_ctx.fine_chan_flags[1]);
        assert!(flag_ctx.fine_chan_flags[2]);
        assert!(!flag_ctx.fine_chan_flags[3]);
        assert!(!flag_ctx.fine_chan_flags[28]);
        assert!(flag_ctx.fine_chan_flags[29]);
        assert!(flag_ctx.fine_chan_flags[30]);
        assert!(flag_ctx.fine_chan_flags[31]);
    }

    /// Refuse to flag edge bandwidth that is not a multiple of fine channel resolution
    #[test]
    #[should_panic]
    #[allow(unused_variables)]
    fn test_error_flagging_width_not_multiple_of_fine_chan_width() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--emulate-cotter",
            "--flag-edge-width", "110",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();
    }

    /// Test that corrections work correctly with `--sel-ants`
    #[test]
    fn test_sel_ants_baselines() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--sel-ants", "1", "2", "3", "4",
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        // check baseline_idxs is the correct size
        assert_eq!(&birli_ctx.vis_sel.baseline_idxs.len(), &10);
        assert_eq!(
            &birli_ctx.vis_sel.baseline_idxs,
            &[128, 129, 130, 131, 255, 256, 257, 381, 382, 506]
        );

        birli_ctx.run().unwrap();
    }

    /// Test `--sel-ants` handles invalid antenna idxs
    #[test]
    fn test_sel_ants_invalid() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--sel-ants", "0", "11", "3", "999",
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    /// Test `--no-sel-flagged-ants` does not include flagged antennas in selection
    #[test]
    fn test_no_sel_flagged_ants() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-antennas", "1", "2",
            "--no-sel-flagged-ants",
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext {
            flag_ctx, vis_sel, ..
        } = BirliContext::from_args(&args).unwrap();

        let flagged_ants = flag_ctx.get_flagged_antenna_idxs();
        assert_eq!(flagged_ants.len(), 2);
        assert!(flagged_ants.contains(&1));
        assert!(flagged_ants.contains(&2));

        assert_eq!(vis_sel.baseline_idxs.len(), 8001);
        assert!(vis_sel.baseline_idxs.contains(&0));
        assert!(!vis_sel.baseline_idxs.contains(&1));
        assert!(!vis_sel.baseline_idxs.contains(&2));
        assert!(vis_sel.baseline_idxs.contains(&3));
        assert!(!vis_sel.baseline_idxs.contains(&128));
        assert!(!vis_sel.baseline_idxs.contains(&129));
        assert!(!vis_sel.baseline_idxs.contains(&130));
        assert!(!vis_sel.baseline_idxs.contains(&254));
        assert!(!vis_sel.baseline_idxs.contains(&255));
        assert!(!vis_sel.baseline_idxs.contains(&256));
        assert!(vis_sel.baseline_idxs.contains(&381));
        assert!(vis_sel.baseline_idxs.contains(&8255));
    }

    /// Test `--no-sel-autos` does not include autos
    #[test]
    fn test_no_sel_autos() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "--no-sel-autos",
            "--no-draw-progress",
            gpufits_paths[0],
            gpufits_paths[1],
        ];

        let BirliContext { vis_sel, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(vis_sel.baseline_idxs.len(), 8128);
        assert!(!vis_sel.baseline_idxs.contains(&0));
        assert!(vis_sel.baseline_idxs.contains(&1));
        assert!(vis_sel.baseline_idxs.contains(&2));
        assert!(vis_sel.baseline_idxs.contains(&3));
        assert!(!vis_sel.baseline_idxs.contains(&128));
        assert!(vis_sel.baseline_idxs.contains(&129));
        assert!(vis_sel.baseline_idxs.contains(&130));
        assert!(vis_sel.baseline_idxs.contains(&254));
        assert!(!vis_sel.baseline_idxs.contains(&255));
        assert!(vis_sel.baseline_idxs.contains(&256));
        assert!(!vis_sel.baseline_idxs.contains(&381));
        assert!(vis_sel.baseline_idxs.contains(&8127));
        assert!(!vis_sel.baseline_idxs.contains(&8255));
    }

    #[test]
    fn test_parse_missing_input() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // no gpufits
        let args = vec!["birli", "-m", metafits_path];

        match BirliContext::from_args(&args) {
            Err(BirliError::ClapError(inner)) => assert!(matches!(
                inner.kind(),
                clap::error::ErrorKind::MissingRequiredArgument
            )),
            Err(e) => panic!("expected missing required argument error, not {e}"),
            Ok(_) => panic!("expected error, but got Ok(_)"),
        }

        // no metafits
        let mut args = vec!["birli"];
        args.extend_from_slice(&gpufits_paths);

        match BirliContext::from_args(&args) {
            Err(BirliError::ClapError(inner)) => assert!(matches!(
                inner.kind(),
                clap::error::ErrorKind::MissingRequiredArgument
            )),
            Err(e) => panic!("expected missing required argument error, not {e}"),
            Ok(_) => panic!("expected error, but got Ok(_)"),
        }
    }

    #[test]
    fn test_parse_invalid_input() {
        let (metafits_path, mut gpufits_paths) = get_1254670392_avg_paths();

        // test nonexistent metafits

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", "nonexistent.metafits"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::MwalibError(_))
        ));

        gpufits_paths[0] = "nonexistent.fits";

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::MwalibError(_))
        ));
    }

    #[test]
    fn test_parse_invalid_time_selection() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "0", "999999"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "2", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_selection() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "1", "2"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { vis_sel, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(vis_sel.timestep_range, 1..3);
    }

    #[test]
    fn test_parse_invalid_time_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-times", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-times", "2", "--flag-init", "0", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.timestep_flags[2]);
        assert!(!flag_ctx.timestep_flags[1]);
    }

    #[test]
    fn test_parse_flag_init_time() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-init", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.timestep_flags[0]);
        assert!(!flag_ctx.timestep_flags[1]);
    }

    #[test]
    fn test_parse_flag_init_steps() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-init-steps", "1", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.timestep_flags[0]);
        assert!(!flag_ctx.timestep_flags[1]);
    }

    #[test]
    fn test_parse_flag_init_invalid_indivisible() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-init", "1",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_flag_end_time() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-init", "0", "--flag-end", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            corr_ctx, flag_ctx, ..
        } = BirliContext::from_args(&args).unwrap();
        let n = corr_ctx.num_common_timesteps;

        assert!(flag_ctx.timestep_flags[n - 1]);
        assert!(!flag_ctx.timestep_flags[n - 2]);
    }

    #[test]
    fn test_parse_flag_end_steps() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-init", "0", "--flag-end-steps", "1", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            corr_ctx, flag_ctx, ..
        } = BirliContext::from_args(&args).unwrap();
        let n = corr_ctx.num_common_timesteps;

        assert!(flag_ctx.timestep_flags[n - 1]);
        assert!(!flag_ctx.timestep_flags[n - 2]);
    }

    #[test]
    fn test_parse_flag_end_invalid_indivisible() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-end", "1",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_invalid_coarse_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-coarse-chans", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_coarse_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-coarse-chans", "2",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.coarse_chan_flags[2]);
        assert!(!flag_ctx.coarse_chan_flags[1]);
    }

    #[test]
    fn test_parse_invalid_fine_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-fine-chans", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_fine_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-fine-chans", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.fine_chan_flags[2]);
        assert!(!flag_ctx.fine_chan_flags[1]);
    }

    #[test]
    fn test_parse_invalid_antenna_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-antennas", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_antenna_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-antennas", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.antenna_flags[2]);
        assert!(!flag_ctx.antenna_flags[1]);
    }

    #[test]
    fn test_parse_no_flag_metafits_and_valid_antenna_flag() {
        let (_, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            // in this metafits file, all ants are flagged.
            "-m", "tests/data/1254670392_avg/1254670392.metafits",
            "--no-flag-metafits",
            "--flag-antennas", "2",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.antenna_flags[2]);
        assert!(!flag_ctx.antenna_flags[1]);
    }

    #[test]
    fn test_parse_flag_autos() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-autos"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.autos);
    }

    #[test]
    fn test_parse_invalid_avg_time() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-factor", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-res", "0.5"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-res", "4s"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::ClapError(_))
        ));
    }

    #[test]
    fn test_parse_invalid_avg_freq() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-factor", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "0.5"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "4s"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::ClapError(_))
        ));
    }

    #[test]
    fn test_parse_valid_avg_freq() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "120"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { avg_freq, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(avg_freq, 3);
    }

    #[test]
    fn test_parse_invalid_time_chunk() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // test when time_chunk is not a multiple of avg_time
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--avg-time-factor", "2",
            "--time-chunk", "1",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_chunk() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--time-chunk", "2"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            num_timesteps_per_chunk,
            ..
        } = BirliContext::from_args(&args).unwrap();

        assert_eq!(num_timesteps_per_chunk, Some(2));
    }

    #[test]
    fn test_parse_invalid_max_memory() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // test when max mem is less than 1 byte
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory",
            "0.0000000000000001",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        // test when max mem can't fit a single averaged timestep
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory", "0.32",
            "--avg-time-factor", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_max_memory() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory", "0.64",
            "--sel-time", "0", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            num_timesteps_per_chunk,
            ..
        } = BirliContext::from_args(&args).unwrap();

        assert_eq!(num_timesteps_per_chunk, Some(2));
    }

    #[test]
    fn test_parse_custom_phase_negative() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--phase-centre", "-1.0", "-2.0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert_eq!(
            birli_ctx.prep_ctx.phase_centre,
            RADec::from_degrees(-1., -2.)
        );
    }

    #[test]
    fn test_parse_sel_range_single() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "2-10", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            channel_range_sel, ..
        } = BirliContext::from_args(&args).unwrap();

        assert!(channel_range_sel.ranges.len() == 1);
        assert!(channel_range_sel.ranges[0].0 == 2);
        assert!(channel_range_sel.ranges[0].1 == 10);
    }

    #[test]
    fn test_parse_sel_two_ranges() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "2-4,6-8", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            channel_range_sel, ..
        } = BirliContext::from_args(&args).unwrap();

        assert!(channel_range_sel.ranges.len() == 2);
        assert!(channel_range_sel.ranges[0].0 == 2);
        assert!(channel_range_sel.ranges[0].1 == 4);
        assert!(channel_range_sel.ranges[1].0 == 6);
        assert!(channel_range_sel.ranges[1].1 == 8);
    }

    #[test]
    fn test_parse_sel_two_ranges_with_single() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "4,6-8", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            channel_range_sel, ..
        } = BirliContext::from_args(&args).unwrap();

        assert!(channel_range_sel.ranges.len() == 2);
        assert!(channel_range_sel.ranges[0].0 == 4);
        assert!(channel_range_sel.ranges[0].1 == 4);
        assert!(channel_range_sel.ranges[1].0 == 6);
        assert!(channel_range_sel.ranges[1].1 == 8);
    }

    #[test]
    fn test_handle_bad_chan_range_single() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "4,6-8,a", "--"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args).err(),
            Some(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_handle_bad_ch_range_double() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "4,6-a", "--"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args).err(),
            Some(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_handle_bad_ch_range_triple() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-chan-ranges", "4,6-7-8", "--"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args).err(),
            Some(BirliError::CLIError(_))
        ));
    }
}

#[cfg(test)]
mod channel_range_tests {
    use marlu::mwalib::CorrelatorContext;

    use crate::{cli::ChannelRanges, test_common::get_1119683928_picket_paths};

    #[test]
    fn test_picket_range() {
        let (metafits_path, gpufits_paths) = get_1119683928_picket_paths();

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        let channel_range_sel = ChannelRanges::all(&corr_ctx);

        assert!(channel_range_sel.ranges.len() == 12);
        assert!(channel_range_sel.ranges[0].0 == 0);
        assert!(channel_range_sel.ranges[0].1 == 1);
        assert!(channel_range_sel.ranges[1].0 == 2);
        assert!(channel_range_sel.ranges[1].1 == 3);
        // ...
        assert!(channel_range_sel.ranges[11].0 == 22);
        assert!(channel_range_sel.ranges[11].1 == 23);
    }
}

/// These are basically integration tests, but if I put them in a separate
/// module, then I don't get the coverage. :(
/// All these tests require aoflagger because they either require flagging or
/// use the --no-rfi option.
/// TODO: get unit test coverage to the point where this can be moved to a unit
/// test module.
#[cfg(test)]
#[cfg(feature = "aoflagger")]
mod tests_aoflagger {
    use std::path::PathBuf;

    use float_cmp::F32Margin;
    use marlu::{
        rubbl_casatables::{Table, TableOpenMode},
        RADec,
    };
    use tempfile::tempdir;

    use crate::{
        test_common::{compare_ms_with_csv, compare_uvfits_with_csv, get_1254670392_avg_paths},
        BirliContext,
    };

    #[test]
    fn test_birli_context_display_doesnt_crash() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("high_2019B_2458765_EOR0"));
        assert!(display.contains("Will not correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will flag with aoflagger"));
        assert!(display.contains("Will not correct geometry"));

        let comment = birli_ctx.prep_ctx.as_comment();

        assert!(!comment.contains("cable length corrections"));
        assert!(!comment.contains("digital gains"));
        assert!(!comment.contains("pfb gains"));
        assert!(comment.contains("aoflagging with"));
        assert!(!comment.contains("geometric corrections"));
    }

    #[test]
    fn compare_cotter_uvfits_nocorrect_rfi() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
            "--no-flag-dc",
            "--flag-init", "0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default(),
            false,
            false,
            true,
            false,
        );
    }

    #[test]
    fn compare_cotter_uvfits_nocorrect_norfi_timechunk1() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.chunked.uvfits");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--time-chunk", "1",
            "--sel-time", "0", "2"
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert_eq!(birli_ctx.prep_ctx.aoflagger_strategy, None);
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(1));

        birli_ctx.run().unwrap();

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
    fn compare_cotter_uvfits_geom_cable_rfi() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-flag-dc",
            "--flag-init", "0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
            false,
            true,
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
    fn compare_cotter_uvfits_geom_cable_rfi_phase_custom() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits.csv",
        );

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--phase-centre", "0.0", "0.0",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-flag-dc",
            "--flag-init", "0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(birli_ctx.prep_ctx.phase_centre, RADec { ra: 0., dec: 0. });
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
            false,
            true,
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
    #[ignore] // https://github.com/MWATelescope/Birli/issues/164
              // the --pointing-centre functionality only provides the correct RA/Dec at
              // the start of the observation, but the actual pointing center is fixed in
              // Azimuth/Elevation coordinates, not RA/Dec.
    fn compare_cotter_ms_geom_cable_rfi_phase_pointing() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms.csv",
        );

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--pointing-centre",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--sel-time", "0", "1",
            "--no-flag-dc",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        let pointing_centre =
            RADec::from_mwalib_tile_pointing(&birli_ctx.corr_ctx.metafits_context);
        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(birli_ctx.prep_ctx.phase_centre, pointing_centre);
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
            true,
            true,
        );
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
    fn compare_cotter_ms_nocorrect_norfi_cal() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.cal.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal", "tests/data/1254670392_avg/1254690096.bin",
            "--no-flag-dc",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    #[test]
    /// Handle when calibration solution is provided with 24 channels, but a subset of channels are provided
    fn compare_cotter_ms_none_norfi_cal_partial() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.cal.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.partial.ms.csv",
        );

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal", "tests/data/1254670392_avg/1254690096.bin",
        ];
        args.extend_from_slice(&gpufits_paths[21..]);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms \
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
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_norfi_nopfb() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.nopfb.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms.csv");

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--no-flag-dc",
            "--flag-init", "0",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct digital gains"));
        assert!(display.contains("Will not flag with aoflagger"));

        let comment = birli_ctx.prep_ctx.as_comment();
        assert!(comment.contains("digital gain"));
        assert!(!comment.contains("aoflag"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(7e-5),
            false,
            true,
        );
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms \
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
    ///  -sbpassband tests/data/subband-passband-32ch-cotter.txt \
    ///  -nostats \
    ///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///  tests/data/1254670392_avg/*gpubox*.fits
    /// ```
    ///
    /// then casa
    ///
    /// ```bash
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_norfi_nodigital_pfb_cotter_40() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.nodigital.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms.csv",
        );

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--no-digital-gains",
            "--pfb-gains", "cotter",
            "--emulate-cotter",
            "--flag-init", "0",
            "--no-flag-dc",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(birli_ctx.prep_ctx.passband_gains.is_some());
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct coarse pfb passband gains"));

        let comment = birli_ctx.prep_ctx.as_comment();
        assert!(comment.contains("pfb gains"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-2),
            false,
            true,
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
    #[ignore]
    fn compare_cotter_ms_corrected() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.dequacked.ms.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-flag-dc",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will flag with aoflagger"));
        assert!(display.contains("Will correct geometry"));

        let comment = birli_ctx.prep_ctx.as_comment();

        assert!(comment.contains("cable length corrections"));
        assert!(!comment.contains("digital gains"));
        assert!(!comment.contains("pfb gains"));
        assert!(comment.contains("aoflagging with"));
        assert!(comment.contains("geometric corrections"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
            false,
        );
    }

    /// Test generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.dequacked.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -noantennapruning \
    ///   -nocablelength \
    ///   -nogeom \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   -timeres 4 \
    ///   -freqres 160 \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// then the following casa commands:
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.dequacked.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.dequacked.ms.csv",
        );

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-res", "4",
            "--avg-freq-res", "160",
            "--no-flag-dc",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }

    /// Same as above but with factors instead of resolution
    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz_factors() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.dequacked.ms.csv",
        );

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--no-flag-dc",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }

    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz_max_mem() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--sel-time", "0", "2",
            "--time-chunk", "2",
            "--no-flag-dc",
            "--flag-init", "0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(2));
        assert_eq!(birli_ctx.vis_sel.timestep_range, 0..3);

        birli_ctx.run().unwrap();

        let main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        assert_eq!(main_table.n_rows(), 2 * 8256);

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }

    /// compare birli uvfits agains pyuvdata, no corrections
    /// Test generated with pyuvdat, see: tests/data/README.md#van-vleck-test-data
    #[test]
    fn compare_pyuvdata_1196175296_mwa_ord_none() {
        let tmp_dir = tempdir().unwrap();
        let uvf_path = tmp_dir.path().join("birli_1196175296.uvfits");
        let metafits_path = PathBuf::from("tests/data/1196175296_mwa_ord/1196175296.metafits");

        let expected_csv_path =
            PathBuf::from("tests/data/1196175296_mwa_ord/pyuvdata_1196175296.none.csv");

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path.to_str().unwrap(),
            "-u", uvf_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-rfi",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-flag-dc",
            "--flag-edge-width", "0",
            "--flag-init", "0",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        println!("{:?}", &args);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path.display().to_string()
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will not correct Van Vleck"));
        assert!(display.contains("Will not correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will not flag with aoflagger"));
        assert!(display.contains("Will not correct geometry"));

        let comment = birli_ctx.prep_ctx.as_comment();

        assert!(!comment.contains("Van Vleck corrections"));
        assert!(!comment.contains("cable length corrections"));
        assert!(!comment.contains("digital gains"));
        assert!(!comment.contains("pfb gains"));
        assert!(!comment.contains("aoflagging with"));
        assert!(!comment.contains("geometric corrections"));

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvf_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
            false,
            true,
            true, // todo: match uvws
        );
    }

    /// compare birli uvfits agains pyuvdata, no corrections
    /// Test generated with pyuvdat, see: tests/data/README.md#van-vleck-test-data
    #[test]
    fn compare_pyuvdata_1254670392_avg_none() {
        let tmp_dir = tempdir().unwrap();
        let uvf_path = tmp_dir.path().join("birli_1254670392.uvfits");
        // let uvf_path = PathBuf::from("birli_1254670392.none.uvfits");
        let metafits_path = PathBuf::from("tests/data/1254670392_avg/1254670392.metafits");

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/pyuvdata_1254670392.none.csv");

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path.to_str().unwrap(),
            "-u", uvf_path.to_str().unwrap(),
            "--sel-time", "0", "0",
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-rfi",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-flag-dc",
            "--flag-edge-width", "0",
            "--flag-init", "0",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits",
        ];
        println!("{:?}", &args);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path.display().to_string()
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will not correct Van Vleck"));
        assert!(display.contains("Will not correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will not flag with aoflagger"));
        assert!(display.contains("Will not correct geometry"));

        let comment = birli_ctx.prep_ctx.as_comment();

        assert!(!comment.contains("Van Vleck corrections"));
        assert!(!comment.contains("cable length corrections"));
        assert!(!comment.contains("digital gains"));
        assert!(!comment.contains("pfb gains"));
        assert!(!comment.contains("aoflagging with"));
        assert!(!comment.contains("geometric corrections"));

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvf_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
            false,
            true,
            true, // todo: match uvws
        );
    }

    /// compare birli uvfits agains VV with no chebychev approximation
    /// Test generated with pyuvdat, see: tests/data/README.md#van-vleck-test-data
    #[test]
    fn compare_pyuvdata_vvnoc() {
        let tmp_dir = tempdir().unwrap();
        let uvf_path = tmp_dir.path().join("1254670392.uvfits");

        let metafits_path = "tests/data/1254670392_avg/1254670392.fixed.metafits";

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/pyuvdata_1254670392.vvnoc.csv");

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvf_path.to_str().unwrap(),
            "--sel-time", "0", "0",
            "--van-vleck",
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-rfi",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-flag-dc",
            "--flag-edge-width", "0",
            "--flag-init", "0",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits",
        ];
        println!("{:?}", &args);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_none());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct Van Vleck"));
        assert!(display.contains("Will not correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will not flag with aoflagger"));
        assert!(display.contains("Will not correct geometry"));

        let comment = birli_ctx.prep_ctx.as_comment();

        assert!(comment.contains("Van Vleck corrections"));
        assert!(!comment.contains("cable length corrections"));
        assert!(!comment.contains("digital gains"));
        assert!(!comment.contains("pfb gains"));
        assert!(!comment.contains("aoflagging with"));
        assert!(!comment.contains("geometric corrections"));

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvf_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
            false,
            true,
            true,
        );
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger_flagset {
    use crate::{io::mwaf::FlagFileSet, BirliContext};
    use itertools::izip;
    use marlu::mwalib::CorrelatorContext;
    use ndarray::prelude::*;
    use tempfile::tempdir;

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

            let left_header = &$left_flagset.header;
            let right_header = &$right_flagset.header;
            assert_eq!(left_header.obs_id, right_header.obs_id);
            assert_eq!(left_header.num_channels, right_header.num_channels);
            assert_eq!(left_header.num_ants, right_header.num_ants);
            assert_eq!(left_header.num_timesteps as usize, num_common_timesteps);
            assert_eq!(left_header.num_pols, right_header.num_pols);
            assert_eq!(left_header.num_rows, right_header.num_rows);
            assert_eq!(left_header.num_rows as usize, num_rows);

            let (left_flags, left_flags_offset) = $left_flagset.read_flags().unwrap().into_raw_vec_and_offset();
            assert_eq!(left_flags_offset.unwrap_or(0), 0);
            let (right_flags, right_flags_offset) = $right_flagset.read_flags().unwrap().into_raw_vec_and_offset();
            assert_eq!(right_flags_offset.unwrap_or(0), 0);
            assert_eq!(left_flags.len(), right_flags.len());
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
        };
    }

    #[test]
    fn aoflagger_outputs_flags() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "none",
            "-f", mwaf_path_template.to_str().unwrap(),
            "--no-flag-dc",
            "--flag-init", "0",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(birli_ctx.prep_ctx.passband_gains.is_none());
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.flag_template,
            Some(mwaf_path_template.to_str().unwrap().into())
        );

        let corr_ctx =
            CorrelatorContext::new(&birli_ctx.io_ctx.metafits_in, &birli_ctx.io_ctx.gpufits_in)
                .unwrap();

        birli_ctx.run().unwrap();

        let gpubox_ids: Vec<usize> = corr_ctx
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| corr_ctx.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let birli_flag_file_set = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();

        let (cotter_flag_file_set, date) = FlagFileSet::open_cotter(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();
        assert_eq!(date, "2020-08-10");

        assert_flagsets_eq!(
            corr_ctx,
            birli_flag_file_set,
            cotter_flag_file_set,
            gpubox_ids
        );
    }

    #[test]
    /// Test data generated on a machine with aoflagger 3.2:
    /// ```bash
    /// cd tests/data/1247842824_flags
    /// birli \
    ///   --pfb-gains none \
    ///   --sel-time 1 1 \
    ///   --flag-init-steps 0 \
    ///   --metafits 1247842824.metafits \
    ///   -f 'FlagfileBirli%%_ts1_ao3.2.mwaf' \
    ///   1247842824_20190722150008_gpubox01_00.fits
    /// birli \
    ///   --pfb-gains none \
    ///   --sel-time 2 2 \
    ///   --flag-init-steps 0 \
    ///   --metafits 1247842824.metafits \
    ///   -f 'FlagfileBirli%%_ts2_ao3.2.mwaf' \
    ///   1247842824_20190722150008_gpubox01_00.fits
    /// ```
    fn aoflagger_outputs_flags_chunked() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        // There are two "chunked" mwaf files; one for timestep 1, another for
        // timestep 2. When running Birli with chunking, the flags for the first
        // chunk should match the first file, and second chunk second file.

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--time-chunk", "1",
            "--sel-time", "1", "2",
            "--flag-init-steps", "0",
            "-f", mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        BirliContext::from_args(&args).unwrap();

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(birli_ctx.prep_ctx.passband_gains.is_none());
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(birli_ctx.prep_ctx.aoflagger_strategy.is_some());
        assert_eq!(
            birli_ctx.io_ctx.metafits_in.display().to_string(),
            metafits_path
        );
        assert_eq!(
            birli_ctx.io_ctx.flag_template,
            Some(mwaf_path_template.to_str().unwrap().into())
        );

        let corr_ctx =
            CorrelatorContext::new(&birli_ctx.io_ctx.metafits_in, &birli_ctx.io_ctx.gpufits_in)
                .unwrap();

        birli_ctx.run().unwrap();

        let gpubox_ids: Vec<usize> = corr_ctx
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| corr_ctx.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let our_flags = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap()
        .read_flags()
        .unwrap();

        let disk_flags_ts1 = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileBirli%%_ts1.mwaf",
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap()
        .read_flags()
        .unwrap();
        let disk_flags_ts2 = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileBirli%%_ts2.mwaf",
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap()
        .read_flags()
        .unwrap();

        // there is a slight difference between aoflagger 3.1 and 3.2 :(
        assert_eq!(
            our_flags.slice(s![0..1, ..20, ..]),
            disk_flags_ts1.slice(s![.., ..20, ..])
        );
        assert_eq!(
            our_flags.slice(s![1..2, ..20, ..]),
            disk_flags_ts2.slice(s![.., ..20, ..])
        );
    }
}
