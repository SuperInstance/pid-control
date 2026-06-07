//! Cascade PID controller (outer + inner loop).

use crate::pid::PidController;

#[allow(clippy::new_without_default)]
/// Cascade controller with two PID loops (outer drives inner setpoint).
#[derive(Debug, Clone)]
pub struct CascadeController {
    /// Outer (primary) loop
    pub outer: PidController,
    /// Inner (secondary) loop
    pub inner: PidController,
}

impl CascadeController {
    /// Create a cascade controller with outer and inner PIDs.
    pub fn new(outer: PidController, inner: PidController) -> Self {
        Self { outer, inner }
    }

    /// Update cascade: outer PID produces setpoint for inner PID.
    pub fn update(&mut self, outer_setpoint: f64, outer_measurement: f64, inner_measurement: f64, dt: f64) -> f64 {
        let inner_setpoint = self.outer.update(outer_setpoint, outer_measurement, dt);
        self.inner.update(inner_setpoint, inner_measurement, dt)
    }

    /// Reset both controllers.
    pub fn reset(&mut self) {
        self.outer.reset();
        self.inner.reset();
    }
}

/// Multi-loop cascade controller (arbitrary depth).
#[derive(Debug, Clone)]
pub struct MultiCascade {
    /// PIDs ordered outer-to-inner
    pub controllers: Vec<PidController>,
}

impl MultiCascade {
    pub fn new(controllers: Vec<PidController>) -> Self {
        Self { controllers }
    }

    /// Update multi-loop cascade.
    /// `measurements` are ordered outer-to-inner (same order as controllers).
    pub fn update(&mut self, setpoint: f64, measurements: &[f64], dt: f64) -> f64 {
        assert_eq!(measurements.len(), self.controllers.len(), "measurements must match controller count");
        let mut sp = setpoint;
        for (ctrl, meas) in self.controllers.iter_mut().zip(measurements.iter()) {
            sp = ctrl.update(sp, *meas, dt);
        }
        sp
    }

    /// Reset all controllers.
    pub fn reset(&mut self) {
        for c in &mut self.controllers {
            c.reset();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cascade_basic() {
        let outer = PidController::new(1.0, 0.1, 0.01);
        let inner = PidController::new(2.0, 0.5, 0.1);
        let mut cascade = CascadeController::new(outer, inner);
        let output = cascade.update(1.0, 0.0, 0.0, 0.01);
        assert!(output != 0.0);
    }

    #[test]
    fn test_cascade_setpoint_tracking() {
        let outer = PidController::new(0.5, 0.2, 0.01);
        let inner = PidController::new(2.0, 1.0, 0.1);
        let mut cascade = CascadeController::new(outer, inner);
        let mut outer_meas = 0.0;
        let mut inner_meas = 0.0;
        let dt = 0.01;
        for _ in 0..3000 {
            let output = cascade.update(1.0, outer_meas, inner_meas, dt);
            inner_meas += output * dt;
            outer_meas = inner_meas; // Outer measurement is inner output
        }
        assert!((outer_meas - 1.0).abs() < 0.1, "outer_meas={}", outer_meas);
    }

    #[test]
    fn test_cascade_reset() {
        let outer = PidController::new(1.0, 1.0, 0.0);
        let inner = PidController::new(1.0, 1.0, 0.0);
        let mut cascade = CascadeController::new(outer, inner);
        cascade.update(1.0, 0.0, 0.0, 0.01);
        cascade.reset();
        assert!(cascade.outer.integral_term().abs() < 1e-9);
        assert!(cascade.inner.integral_term().abs() < 1e-9);
    }

    #[test]
    fn test_multi_cascade() {
        let pids = vec![
            PidController::new(0.5, 0.1, 0.01),
            PidController::new(1.0, 0.2, 0.05),
            PidController::new(2.0, 0.5, 0.1),
        ];
        let mut mc = MultiCascade::new(pids);
        let measurements = [0.0, 0.0, 0.0];
        let output = mc.update(1.0, &measurements, 0.01);
        assert!(output != 0.0);
    }

    #[test]
    fn test_cascade_disturbance_rejection_inner() {
        let outer = PidController::new(0.5, 0.2, 0.01);
        let inner = PidController::new(3.0, 1.0, 0.2);
        let mut cascade = CascadeController::new(outer, inner);
        let mut outer_meas = 0.0;
        let mut inner_meas = 0.0;
        let dt = 0.01;
        // Reach steady state
        for _ in 0..2000 {
            let output = cascade.update(1.0, outer_meas, inner_meas, dt);
            inner_meas += output * dt;
            outer_meas = inner_meas;
        }
        // Inner disturbance
        inner_meas -= 0.3;
        for _ in 0..1000 {
            let output = cascade.update(1.0, outer_meas, inner_meas, dt);
            inner_meas += output * dt;
            outer_meas = inner_meas;
        }
        assert!((outer_meas - 1.0).abs() < 0.15, "after inner disturbance: {}", outer_meas);
    }
}
