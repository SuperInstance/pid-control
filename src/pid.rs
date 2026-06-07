//! Core PID controller implementation.

/// PID controller with configurable gains and output limits.
#[derive(Debug, Clone)]
pub struct PidController {
    /// Proportional gain
    pub kp: f64,
    /// Integral gain
    pub ki: f64,
    /// Derivative gain
    pub kd: f64,
    /// Output minimum
    pub out_min: f64,
    /// Output maximum
    pub out_max: f64,
    /// Setpoint weight for proportional term (0..1)
    pub setpoint_weight: f64,

    // Internal state
    integral: f64,
    prev_error: f64,
    prev_measurement: f64,
    first_update: bool,
}

impl PidController {
    /// Create a new PID controller with the given gains.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            out_min: f64::NEG_INFINITY,
            out_max: f64::INFINITY,
            setpoint_weight: 1.0,
            integral: 0.0,
            prev_error: 0.0,
            prev_measurement: 0.0,
            first_update: true,
        }
    }

    /// Set output limits.
    pub fn with_output_limits(mut self, min: f64, max: f64) -> Self {
        self.out_min = min;
        self.out_max = max;
        self
    }

    /// Set the setpoint weight for the proportional term (derivative-on-measurement).
    pub fn with_setpoint_weight(mut self, w: f64) -> Self {
        self.setpoint_weight = w.clamp(0.0, 1.0);
        self
    }

    /// Compute one PID step, returning the control output.
    pub fn update(&mut self, setpoint: f64, measurement: f64, dt: f64) -> f64 {
        let error = setpoint - measurement;

        // Proportional term (with setpoint weighting)
        let p_term = self.kp * (self.setpoint_weight * setpoint - measurement);

        // Integral term
        self.integral += error * dt;

        // Derivative term (on measurement to avoid derivative kick)
        let d_term = if self.first_update {
            self.first_update = false;
            0.0
        } else {
            self.kd * -(measurement - self.prev_measurement) / dt
        };

        // Raw output before clamping
        let output = p_term + self.ki * self.integral + d_term;

        // Clamp output
        let clamped = output.clamp(self.out_min, self.out_max);

        // Anti-windup: back-calculate integral if output was clamped
        if output != clamped && self.ki != 0.0 {
            self.integral -= (output - clamped) / self.ki;
        }

        self.prev_error = error;
        self.prev_measurement = measurement;

        clamped
    }

    /// Reset the controller state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.prev_measurement = 0.0;
        self.first_update = true;
    }

    /// Get the current integral accumulator.
    pub fn integral_term(&self) -> f64 {
        self.integral
    }

    /// Get the last error value.
    pub fn last_error(&self) -> f64 {
        self.prev_error
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pid_basic_creation() {
        let pid = PidController::new(1.0, 0.1, 0.01);
        assert!((pid.kp - 1.0).abs() < 1e-9);
        assert!((pid.ki - 0.1).abs() < 1e-9);
        assert!((pid.kd - 0.01).abs() < 1e-9);
    }

    #[test]
    fn test_setpoint_tracking() {
        let mut pid = PidController::new(2.0, 1.0, 0.5);
        let mut measurement = 0.0;
        let dt = 0.01;
        for _ in 0..1000 {
            let output = pid.update(1.0, measurement, dt);
            // Simulate simple plant: first-order system
            measurement += output * dt;
        }
        // Should converge close to setpoint
        assert!((measurement - 1.0).abs() < 0.05, "measurement={}", measurement);
    }

    #[test]
    fn test_output_clamping() {
        let mut pid = PidController::new(100.0, 0.0, 0.0).with_output_limits(-10.0, 10.0);
        let output = pid.update(1.0, 0.0, 0.01);
        // P = 100 * 1.0 = 100, but clamped to 10
        assert!((output - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_disturbance_rejection() {
        let mut pid = PidController::new(2.0, 0.5, 0.1);
        let mut measurement = 0.0;
        let dt = 0.01;
        // Run to steady state
        for _ in 0..500 {
            let output = pid.update(1.0, measurement, dt);
            measurement += output * dt;
        }
        // Apply disturbance
        measurement -= 0.5;
        // Let controller recover
        for _ in 0..500 {
            let output = pid.update(1.0, measurement, dt);
            measurement += output * dt;
        }
        assert!((measurement - 1.0).abs() < 0.1, "after disturbance measurement={}", measurement);
    }

    #[test]
    fn test_reset_clears_state() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        pid.update(1.0, 0.0, 0.01);
        pid.update(1.0, 0.5, 0.01);
        assert!(pid.integral_term().abs() > 0.0);
        pid.reset();
        assert!(pid.integral_term().abs() < 1e-9);
        assert!(pid.first_update);
    }

    #[test]
    fn test_zero_gains_produce_zero_or_only_integral() {
        let mut pid = PidController::new(0.0, 0.0, 0.0);
        let out = pid.update(1.0, 0.0, 0.01);
        assert!(out.abs() < 1e-9);
    }

    #[test]
    fn test_pure_p_controller_steady_state_error() {
        let mut pid = PidController::new(1.0, 0.0, 0.0);
        let mut measurement = 0.0;
        let dt = 0.01;
        for _ in 0..2000 {
            let output = pid.update(1.0, measurement, dt);
            // Plant with unit gain
            measurement += (output - measurement) * dt;
        }
        // Pure P controller has steady-state error
        let error = (1.0 - measurement).abs();
        assert!(error > 0.001, "pure P should have steady-state error, error={}", error);
    }

    #[test]
    fn test_pi_eliminated_steady_state_error() {
        let mut pid = PidController::new(1.0, 0.5, 0.0);
        let mut measurement = 0.0;
        let dt = 0.01;
        for _ in 0..5000 {
            let output = pid.update(1.0, measurement, dt);
            measurement += (output - measurement) * dt;
        }
        assert!((measurement - 1.0).abs() < 0.01, "PI should eliminate steady-state error, got {}", measurement);
    }
}
