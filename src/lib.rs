//! # pid-control
//!
//! A research-grade PID controller library with anti-windup, bumpless transfer,
//! cascade control, and Ziegler-Nichols auto-tuning. Pure Rust, no external dependencies.

pub mod pid;
pub mod anti_windup;
pub mod tuning;
pub mod cascade;
pub mod auto_tune;

pub use pid::PidController;
pub use anti_windup::AntiWindupConfig;
pub use cascade::CascadeController;
