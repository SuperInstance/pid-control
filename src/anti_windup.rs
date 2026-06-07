//! Anti-windup strategies for PID controllers.

/// Anti-windup configuration.
#[derive(Debug, Clone)]
pub struct AntiWindupConfig {
    /// Integral term minimum
    pub integral_min: f64,
    /// Integral term maximum
    pub integral_max: f64,
    /// Back-calculation coefficient
    pub back_calc_coeff: f64,
    /// Enable conditional integration
    pub conditional: bool,
}

impl AntiWindupConfig {
    /// Create clamping anti-windup with given integral limits.
    pub fn clamping(min: f64, max: f64) -> Self {
        Self {
            integral_min: min,
            integral_max: max,
            back_calc_coeff: 1.0,
            conditional: false,
        }
    }

    /// Create back-calculation anti-windup.
    pub fn back_calc(coeff: f64) -> Self {
        Self {
            integral_min: f64::NEG_INFINITY,
            integral_max: f64::INFINITY,
            back_calc_coeff: coeff,
            conditional: false,
        }
    }

    /// Create conditional integration anti-windup.
    pub fn conditional() -> Self {
        Self {
            integral_min: f64::NEG_INFINITY,
            integral_max: f64::INFINITY,
            back_calc_coeff: 1.0,
            conditional: true,
        }
    }

    /// Apply clamping to an integral value.
    pub fn clamp_integral(&self, integral: f64) -> f64 {
        integral.clamp(self.integral_min, self.integral_max)
    }

    /// Compute back-calculation correction.
    pub fn back_calc_correction(&self, saturated: bool, output_diff: f64) -> f64 {
        if saturated {
            output_diff * self.back_calc_coeff
        } else {
            0.0
        }
    }

    /// Check if integration should be conditional.
    pub fn should_integrate(&self, error: f64, output: f64, out_min: f64, out_max: f64) -> bool {
        if !self.conditional {
            return true;
        }
        // Don't integrate if output is saturated and error would make it worse
        let saturated_high = output >= out_max;
        let saturated_low = output <= out_min;
        !(saturated_high && error > 0.0 || saturated_low && error < 0.0)
    }
}

/// PID controller with explicit anti-windup support.
#[derive(Debug, Clone)]
pub struct AntiWindupPid {
    pub kp: f64,
    pub ki: f64,
    pub kd: f64,
    pub out_min: f64,
    pub out_max: f64,
    pub config: AntiWindupConfig,

    integral: f64,
    prev_measurement: f64,
    first: bool,
}

impl AntiWindupPid {
    pub fn new(kp: f64, ki: f64, kd: f64, config: AntiWindupConfig) -> Self {
        Self {
            kp, ki, kd,
            out_min: f64::NEG_INFINITY,
            out_max: f64::INFINITY,
            config,
            integral: 0.0,
            prev_measurement: 0.0,
            first: true,
        }
    }

    pub fn with_output_limits(mut self, min: f64, max: f64) -> Self {
        self.out_min = min;
        self.out_max = max;
        self
    }

    /// Bumpless transfer: set the integral to match a given output.
    pub fn bumpless_transfer(&mut self, target_output: f64, measurement: f64) {
        // output ≈ kp*(setpoint - measurement) + ki*integral
        // We don't know setpoint here, so set integral for zero P contribution
        if self.ki != 0.0 {
            self.integral = target_output / self.ki;
        }
        self.prev_measurement = measurement;
    }

    pub fn update(&mut self, setpoint: f64, measurement: f64, dt: f64) -> f64 {
        let error = setpoint - measurement;

        let p_term = self.kp * error;

        // Conditional integration check
        let raw_i = self.integral + error * dt;
        let should = self.config.should_integrate(error, p_term, self.out_min, self.out_max);
        if should {
            self.integral = self.config.clamp_integral(raw_i);
        }

        let d_term = if self.first {
            self.first = false;
            0.0
        } else {
            self.kd * -(measurement - self.prev_measurement) / dt
        };

        let output = p_term + self.ki * self.integral + d_term;
        let clamped = output.clamp(self.out_min, self.out_max);

        // Back-calculation correction
        if output != clamped {
            let correction = self.config.back_calc_correction(true, output - clamped);
            if self.ki != 0.0 {
                self.integral -= correction / self.ki;
                self.integral = self.config.clamp_integral(self.integral);
            }
        }

        self.prev_measurement = measurement;
        clamped
    }

    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_measurement = 0.0;
        self.first = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clamping_anti_windup() {
        let config = AntiWindupConfig::clamping(-5.0, 5.0);
        let mut pid = AntiWindupPid::new(1.0, 10.0, 0.0, config)
            .with_output_limits(-10.0, 10.0);
        // Large error should saturate but integral should be clamped
        for _ in 0..100 {
            pid.update(10.0, 0.0, 0.01);
        }
        assert!(pid.integral.abs() <= 5.1, "integral clamped: {}", pid.integral);
    }

    #[test]
    fn test_back_calc_anti_windup() {
        let config = AntiWindupConfig::back_calc(1.0);
        let mut pid = AntiWindupPid::new(1.0, 5.0, 0.0, config)
            .with_output_limits(-10.0, 10.0);
        for _ in 0..50 {
            pid.update(5.0, 0.0, 0.01);
        }
        // After prolonged saturation, integral should not grow unboundedly
        assert!(pid.integral.abs() < 100.0);
    }

    #[test]
    fn test_conditional_integration() {
        let config = AntiWindupConfig::conditional();
        let mut pid = AntiWindupPid::new(1.0, 5.0, 0.0, config)
            .with_output_limits(-10.0, 10.0);
        // Drive to saturation
        for _ in 0..200 {
            pid.update(100.0, 0.0, 0.01);
        }
        // Integral should stop growing once saturated
        let int_after_sat = pid.integral;
        for _ in 0..100 {
            pid.update(100.0, 0.0, 0.01);
        }
        // Should not grow much more (conditional stops integration)
        assert!((pid.integral - int_after_sat).abs() < 1.0);
    }

    #[test]
    fn test_bumpless_transfer() {
        let config = AntiWindupConfig::clamping(-100.0, 100.0);
        let mut pid = AntiWindupPid::new(1.0, 1.0, 0.0, config);
        pid.bumpless_transfer(5.0, 0.0);
        // Integral should be set so that ki*integral ≈ 5.0
        assert!((pid.ki * pid.integral - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_no_windup_when_not_saturated() {
        let config = AntiWindupConfig::clamping(-100.0, 100.0);
        let mut pid = AntiWindupPid::new(1.0, 1.0, 0.0, config);
        let mut integral_grew = false;
        let prev_int = pid.integral;
        pid.update(0.1, 0.0, 0.01);
        if pid.integral > prev_int {
            integral_grew = true;
        }
        assert!(integral_grew);
    }
}
