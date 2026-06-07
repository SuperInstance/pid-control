# pid-control

Research-grade PID controller library for Rust with anti-windup, bumpless transfer, cascade control, and Ziegler-Nichols auto-tuning.

## Features

- **Core PID**: Proportional, integral, derivative terms with configurable gains and output limits
- **Anti-Windup**: Clamping, back-calculation, and conditional integration strategies
- **Bumpless Transfer**: Seamless mode switching without output spikes
- **Tuning Methods**: Ziegler-Nichols (classic, some-overshoot, no-overshoot), Tyreus-Luyben, Cohen-Coon, IMC
- **Cascade Control**: Dual-loop and multi-loop cascade PID architectures
- **Auto-Tuning**: Relay feedback method for automatic PID parameter estimation
- **Zero Dependencies**: Pure `std` Rust, no external crates

## Usage

```rust
use pid_control::PidController;

let mut pid = PidController::new(2.0, 0.5, 0.1)
    .with_output_limits(-10.0, 10.0);

let dt = 0.01;
let mut measurement = 0.0;
for _ in 0..1000 {
    let output = pid.update(1.0, measurement, dt);
    measurement += output * dt;
}
```

## Modules

- `pid` — Core PID controller
- `anti_windup` — Anti-windup strategies
- `tuning` — Ziegler-Nichols and other tuning methods
- `cascade` — Cascade and multi-loop controllers
- `auto_tune` — Relay feedback auto-tuning

## License

MIT
