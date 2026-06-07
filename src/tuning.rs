//! Ziegler-Nichols and other tuning methods.

/// Tuning parameters (Kp, Ki, Kd).
#[derive(Debug, Clone, Copy)]
pub struct TuningParams {
    pub kp: f64,
    pub ki: f64,
    pub kd: f64,
}

/// Result of a step response analysis.
#[derive(Debug, Clone)]
pub struct StepResponse {
    /// Process gain
    pub gain: f64,
    /// Time constant (seconds)
    pub time_constant: f64,
    /// Dead time (seconds)
    pub dead_time: f64,
}

impl StepResponse {
    /// Estimate from step test data. `data` is (time, value) pairs after a unit step input.
    pub fn from_step_data(data: &[(f64, f64)], step_size: f64) -> Option<Self> {
        if data.len() < 3 {
            return None;
        }
        let final_val = data.last()?.1;
        let gain = final_val / step_size;

        // Find 63.2% point for time constant
        let target = gain * step_size * 0.632;
        let mut tc_idx = 0;
        for (i, &(_, v)) in data.iter().enumerate() {
            if v >= target {
                tc_idx = i;
                break;
            }
        }
        let time_constant = data.get(tc_idx)?.0;

        // Estimate dead time as time to first significant response
        let threshold = gain * step_size * 0.01;
        let mut dead_time = 0.0;
        for &(t, v) in data {
            if v > threshold {
                dead_time = t;
                break;
            }
        }

        Some(StepResponse { gain, time_constant, dead_time })
    }
}

/// Ziegler-Nichols tuning methods.
pub struct ZieglerNichols;

impl ZieglerNichols {
    /// Classic Ziegler-Nichols PID tuning from ultimate gain and period.
    /// `ku` = ultimate gain, `tu` = ultimate period.
    pub fn pid_classic(ku: f64, tu: f64) -> TuningParams {
        TuningParams {
            kp: 0.6 * ku,
            ki: 1.2 * ku / tu,
            kd: 0.075 * ku * tu,
        }
    }

    /// Ziegler-Nichols PI tuning.
    pub fn pi_classic(ku: f64, tu: f64) -> TuningParams {
        TuningParams {
            kp: 0.45 * ku,
            ki: 0.54 * ku / tu,
            kd: 0.0,
        }
    }

    /// Ziegler-Nichols P-only tuning.
    pub fn p_classic(ku: f64) -> TuningParams {
        TuningParams {
            kp: 0.5 * ku,
            ki: 0.0,
            kd: 0.0,
        }
    }

    /// Some-overshoot PID (modified Z-N).
    pub fn pid_some_overshoot(ku: f64, tu: f64) -> TuningParams {
        TuningParams {
            kp: 0.33 * ku,
            ki: 0.66 * ku / tu,
            kd: 0.105 * ku * tu,
        }
    }

    /// No-overshoot PID (modified Z-N).
    pub fn pid_no_overshoot(ku: f64, tu: f64) -> TuningParams {
        TuningParams {
            kp: 0.2 * ku,
            ki: 0.4 * ku / tu,
            kd: 0.066 * ku * tu,
        }
    }

    /// Tyreus-Luyben PID tuning (more conservative than Z-N).
    pub fn tyreus_luyben_pid(ku: f64, tu: f64) -> TuningParams {
        TuningParams {
            kp: ku / 3.2,
            ki: ku / (2.2 * tu),
            kd: ku * tu / 12.0,
        }
    }

    /// Cohen-Coon tuning from process reaction curve.
    pub fn cohen_coon(gain: f64, time_const: f64, dead_time: f64) -> TuningParams {
        let r = dead_time / time_const;
        let kc = (1.0 / gain) * (time_const / dead_time) * (0.9 + r / 12.0);
        let ti = dead_time * (30.0 + 3.0 * r) / (9.0 + 20.0 * r);
        let td = dead_time * (4.0 + 11.0 * r) / (30.0 * (1.0 + 4.0 * r));
        TuningParams {
            kp: kc,
            ki: kc / ti,
            kd: kc * td,
        }
    }

    /// Internal Model Control (IMC) PID tuning.
    pub fn imc_pid(gain: f64, time_const: f64, dead_time: f64, lambda: f64) -> TuningParams {
        let tc = f64::max(time_const, dead_time / 2.0);
        let kc = tc / (gain * (lambda + dead_time / 2.0));
        let ti = tc;
        let td = dead_time / 2.0;
        TuningParams {
            kp: kc,
            ki: kc / ti,
            kd: kc * td,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zn_pid_classic() {
        let t = ZieglerNichols::pid_classic(10.0, 2.0);
        assert!((t.kp - 6.0).abs() < 1e-9);
        assert!((t.ki - 6.0).abs() < 1e-9);
        assert!((t.kd - 1.5).abs() < 1e-9);
    }

    #[test]
    fn test_zn_pi_classic() {
        let t = ZieglerNichols::pi_classic(10.0, 2.0);
        assert!((t.kp - 4.5).abs() < 1e-9);
        assert!((t.ki - 2.7).abs() < 1e-9);
        assert!((t.kd).abs() < 1e-9);
    }

    #[test]
    fn test_zn_p_classic() {
        let t = ZieglerNichols::p_classic(10.0);
        assert!((t.kp - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_cohen_coon_tuning() {
        let t = ZieglerNichols::cohen_coon(1.0, 10.0, 2.0);
        assert!(t.kp > 0.0);
        assert!(t.ki > 0.0);
        assert!(t.kd > 0.0);
    }

    #[test]
    fn test_imc_tuning() {
        let t = ZieglerNichols::imc_pid(1.0, 10.0, 2.0, 1.0);
        assert!(t.kp > 0.0);
        assert!(t.ki > 0.0);
        assert!(t.kd > 0.0);
    }

    #[test]
    fn test_step_response_analysis() {
        let data: Vec<(f64, f64)> = vec![
            (0.0, 0.0), (0.5, 0.0), (1.0, 0.1), (2.0, 0.4),
            (3.0, 0.7), (5.0, 0.9), (8.0, 0.98), (10.0, 1.0),
        ];
        let resp = StepResponse::from_step_data(&data, 1.0).unwrap();
        assert!(resp.gain > 0.0);
        assert!(resp.time_constant > 0.0);
    }

    #[test]
    fn test_tyreus_luyben() {
        let t = ZieglerNichols::tyreus_luyben_pid(10.0, 2.0);
        assert!(t.kp > 0.0);
        assert!(t.ki > 0.0);
        assert!(t.kd > 0.0);
        // More conservative: lower Kp than classic Z-N
        let classic = ZieglerNichols::pid_classic(10.0, 2.0);
        assert!(t.kp < classic.kp);
    }

    #[test]
    fn test_no_overshoot_more_conservative() {
        let no_os = ZieglerNichols::pid_no_overshoot(10.0, 2.0);
        let some_os = ZieglerNichols::pid_some_overshoot(10.0, 2.0);
        assert!(no_os.kp < some_os.kp);
    }
}
