use birli::cli::main_with_args;
use log::trace;
use std::env;

fn main() {
    env_logger::init_from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    );
    trace!("start main");
    main_with_args(env::args());
    trace!("end main");
}
