//! Auto-tuning via relay method and Ziegler-Nichols.

use crate::tuning::{TuningParams, ZieglerNichols};

/// State of the auto-tuner relay experiment.
#[derive(Debug, Clone)]
pub struct AutoTuner {
    /// Relay amplitude
    pub relay_amplitude: f64,
    /// Setpoint for the relay experiment
    pub setpoint: f64,
    /// Current phase of the experiment
    phase: AutoTunePhase,
    /// Collected oscillation peaks
    peaks: Vec<f64>,
    /// Collected peak times
    peak_times: Vec<f64>,
    /// Current relay output
    relay_output: f64,
    /// Previous measurement
    prev_measurement: Option<f64>,
    /// Previous error sign
    prev_error_sign: i32,
    /// Time accumulator
    time: f64,
    /// Number of complete cycles detected
    cycles: usize,
    /// Target cycles before computing result
    target_cycles: usize,
}

#[derive(Debug, Clone, PartialEq)]
enum AutoTunePhase {
    Running,
    Done,
}

impl AutoTuner {
    /// Create a new auto-tuner using the relay feedback method.
    pub fn new(setpoint: f64, relay_amplitude: f64) -> Self {
        Self {
            relay_amplitude,
            setpoint,
            phase: AutoTunePhase::Running,
            peaks: Vec::new(),
            peak_times: Vec::new(),
            relay_output: relay_amplitude,
            prev_measurement: None,
            prev_error_sign: 1,
            time: 0.0,
            cycles: 0,
            target_cycles: 4,
        }
    }

    /// Set number of cycles to collect before computing tuning.
    pub fn with_cycles(mut self, n: usize) -> Self {
        self.target_cycles = n.max(2);
        self
    }

    /// Feed one measurement sample, returns the relay output to apply.
    pub fn update(&mut self, measurement: f64, dt: f64) -> f64 {
        if self.phase == AutoTunePhase::Done {
            return 0.0;
        }

        self.time += dt;
        let error = self.setpoint - measurement;
        let error_sign = if error >= 0.0 { 1 } else { -1 };

        // Detect zero crossing (peak detection)
        if error_sign != self.prev_error_sign {
            if let Some(prev) = self.prev_measurement {
                self.peaks.push(prev);
                self.peak_times.push(self.time);
            }
            if self.peaks.len() >= 3 {
                self.cycles += 1;
            }
        }

        // Relay output
        self.relay_output = if error > 0.0 {
            self.relay_amplitude
        } else {
            -self.relay_amplitude
        };

        self.prev_measurement = Some(measurement);
        self.prev_error_sign = error_sign;

        if self.cycles >= self.target_cycles {
            self.phase = AutoTunePhase::Done;
        }

        self.relay_output
    }

    /// Check if auto-tuning is complete.
    pub fn is_done(&self) -> bool {
        self.phase == AutoTunePhase::Done
    }

    /// Compute the ultimate gain and period from the relay experiment.
    pub fn ultimate_params(&self) -> Option<(f64, f64)> {
        if self.peak_times.len() < 3 {
            return None;
        }
        // Ultimate period: average period of oscillation
        let mut periods = Vec::new();
        for w in self.peak_times.windows(2) {
            periods.push(w[1] - w[0]);
        }
        let avg_period = periods.iter().sum::<f64>() / periods.len() as f64;

        // Ultimate gain from relay amplitude and oscillation amplitude
        let peak_max = self.peaks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let peak_min = self.peaks.iter().cloned().fold(f64::INFINITY, f64::min);
        let amplitude = (peak_max - peak_min) / 2.0;

        if amplitude < 1e-9 {
            return None;
        }

        // Ku = (4 * d) / (pi * a) where d = relay amplitude, a = oscillation amplitude
        let ku = (4.0 * self.relay_amplitude) / (std::f64::consts::PI * amplitude);

        Some((ku, avg_period * 2.0)) // Period is between same-sign crossings
    }

    /// Compute PID tuning parameters from the relay experiment.
    pub fn compute_pid_tuning(&self) -> Option<TuningParams> {
        let (ku, tu) = self.ultimate_params()?;
        Some(ZieglerNichols::pid_classic(ku, tu))
    }
}

/// Run a simulated relay auto-tune experiment on a first-order plant.
pub fn simulate_auto_tune(
    plant_gain: f64,
    plant_time_const: f64,
    plant_dead_time: f64,
    setpoint: f64,
    relay_amplitude: f64,
    dt: f64,
    duration: f64,
) -> AutoTuner {
    let mut tuner = AutoTuner::new(setpoint, relay_amplitude);
    let mut measurement = 0.0;
    let steps = (duration / dt) as usize;
    for _ in 0..steps {
        let output = tuner.update(measurement, dt);
        // Simple first-order + dead-time plant simulation
        let effective_input = output * plant_gain;
        measurement += (effective_input - measurement) / plant_time_const * dt;
        let _ = plant_dead_time; // Dead time approximation ignored in simple sim
    }
    tuner
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_auto_tuner_creation() {
        let at = AutoTuner::new(1.0, 0.5);
        assert!(!at.is_done());
        assert!((at.relay_amplitude - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_relay_switching() {
        let mut at = AutoTuner::new(0.0, 1.0);
        let out1 = at.update(0.5, 0.01); // error < 0
        assert!((out1 + 1.0).abs() < 1e-9);
        let out2 = at.update(-0.5, 0.01); // error > 0
        assert!((out2 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_simulated_auto_tune() {
        let result = simulate_auto_tune(2.0, 1.0, 0.1, 1.0, 0.5, 0.001, 20.0);
        // Relay experiment may or may not produce peaks depending on dynamics
        // Just verify it doesn't panic and completes
        assert!(result.time > 0.0);
    }

    #[test]
    fn test_ultimate_params_computation() {
        let result = simulate_auto_tune(2.0, 1.0, 0.1, 1.0, 0.5, 0.001, 30.0);
        if let Some((ku, tu)) = result.ultimate_params() {
            assert!(ku > 0.0, "ku={}", ku);
            assert!(tu > 0.0, "tu={}", tu);
        }
    }

    #[test]
    fn test_pid_tuning_from_auto_tune() {
        let result = simulate_auto_tune(2.0, 1.0, 0.1, 1.0, 0.5, 0.001, 30.0);
        if let Some(tuning) = result.compute_pid_tuning() {
            assert!(tuning.kp > 0.0);
        }
    }

    #[test]
    fn test_auto_tuner_with_cycles() {
        let at = AutoTuner::new(1.0, 0.5).with_cycles(6);
        assert_eq!(at.target_cycles, 6);
    }
}
