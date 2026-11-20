use std::f64::consts::PI;

use crate::circuit::{Circuit, SimulatorError, StateVectorSimulator, apply_gate};
use crate::state::StateVector;

#[derive(Clone, Debug)]
pub struct QftResult {
    pub circuit: Circuit,
    pub state: StateVector,
}

/// Build and simulate QFT on `targets`. `targets[0]` is the most significant bit.
pub fn qft(targets: &[usize], qubits: usize) -> Result<QftResult, QftError> {
    build_qft(targets, qubits, false)
}

/// Build the inverse QFT, which is the adjoint of `qft`.
pub fn inverse_qft(targets: &[usize], qubits: usize) -> Result<QftResult, QftError> {
    build_qft(targets, qubits, true)
}

fn build_qft(targets: &[usize], qubits: usize, inverse: bool) -> Result<QftResult, QftError> {
    let mut circuit = Circuit::new(qubits);
    let mut seen = Vec::new();
    for &q in targets {
        if q >= qubits {
            return Err(QftError::QubitOutOfRange { index: q, qubits });
        }
        if seen.contains(&q) {
            return Err(QftError::DuplicateQubit { index: q });
        }
        seen.push(q);
    }

    for stage in 0..targets.len() {
        circuit.h(targets[stage]);
        for offset in 1..targets.len() - stage {
            let angle = PI / (1usize << offset) as f64;
            let controlled_angle = if inverse { -angle } else { angle };
            circuit.cp(targets[stage + offset], targets[stage], controlled_angle);
        }
    }
    for i in 0..targets.len() / 2 {
        circuit.swap(targets[i], targets[targets.len() - 1 - i]);
    }

    let state = StateVectorSimulator::new()
        .run_state(&circuit)
        .map_err(QftError::Simulator)?;
    Ok(QftResult { circuit, state })
}

/// Effective number of Fourier modes above a probability threshold.
pub fn qft_sparsity(result: &QftResult, threshold: f64) -> usize {
    result
        .state
        .amplitudes()
        .iter()
        .filter(|amplitude| amplitude.norm_sqr() > threshold)
        .count()
}

/// Gate density: average number of gates per QFT target.
pub fn qft_gate_density(result: &QftResult) -> f64 {
    let targets = result.circuit.qubits().max(1);
    result.circuit.gates().len() as f64 / targets as f64
}

/// Average basis-state probability of a QFT output.
pub fn qft_average_probability(result: &QftResult) -> f64 {
    let amplitudes = result.state.amplitudes();
    if amplitudes.is_empty() {
        return 0.0;
    }
    amplitudes.iter().map(|a| a.norm_sqr()).sum::<f64>() / amplitudes.len() as f64
}

/// Ratio of QFT gates to inverse-QFT gates; returns 1.0 when inverse has no gates.
pub fn qft_gate_symmetry(result: &QftResult, inverse: &QftResult) -> f64 {
    if inverse.circuit.gates().is_empty() {
        return 1.0;
    }
    result.circuit.gates().len() as f64 / inverse.circuit.gates().len() as f64
}

/// Gate-count balance between forward and inverse QFT: `1 - |a-b| / max(a,b)`.
pub fn qft_gate_balance(result: &QftResult, inverse: &QftResult) -> f64 {
    let forward = result.circuit.gates().len();
    let backward = inverse.circuit.gates().len();
    let largest = forward.max(backward);
    if largest == 0 {
        return 1.0;
    }
    1.0 - (forward.abs_diff(backward) as f64 / largest as f64)
}

/// Average gate density across the QFT and inverse-QFT circuits.
pub fn qft_density_balance(result: &QftResult, inverse: &QftResult) -> f64 {
    (qft_gate_density(result) + qft_gate_density(inverse)) / 2.0
}

/// Phase symmetry score: `1 - |mean phase of odd modes - mean phase of even modes| / pi`.
pub fn qft_phase_balance(result: &QftResult) -> f64 {
    let mut even_sum = 0.0;
    let mut even_count = 0usize;
    let mut odd_sum = 0.0;
    let mut odd_count = 0usize;
    for (index, amplitude) in result.state.amplitudes().iter().enumerate() {
        if amplitude.norm_sqr() <= f64::EPSILON {
            continue;
        }
        if index % 2 == 0 {
            even_sum += amplitude.arg();
            even_count += 1;
        } else {
            odd_sum += amplitude.arg();
            odd_count += 1;
        }
    }
    if even_count == 0 || odd_count == 0 {
        return 1.0;
    }
    let gap = (even_sum / even_count as f64 - odd_sum / odd_count as f64).abs();
    (1.0 - gap / std::f64::consts::PI).clamp(0.0, 1.0)
}

/// Largest phase gap between adjacent Fourier modes, normalized by pi.
pub fn qft_adjacent_phase_gap(result: &QftResult) -> f64 {
    let phases: Vec<f64> = result
        .state
        .amplitudes()
        .iter()
        .filter_map(|amplitude| (amplitude.norm_sqr() > f64::EPSILON).then(|| amplitude.arg()))
        .collect();
    phases
        .windows(2)
        .map(|pair| (pair[0] - pair[1]).abs())
        .fold(0.0_f64, f64::max)
        / std::f64::consts::PI
}

/// Average absolute phase difference between adjacent Fourier modes, normalized by pi.
pub fn qft_adjacent_phase_step(result: &QftResult) -> f64 {
    let phases: Vec<f64> = result
        .state
        .amplitudes()
        .iter()
        .filter_map(|amplitude| (amplitude.norm_sqr() > f64::EPSILON).then(|| amplitude.arg()))
        .collect();
    if phases.len() < 2 {
        return 0.0;
    }
    phases
        .windows(2)
        .map(|pair| (pair[1] - pair[0]).abs())
        .sum::<f64>()
        / (phases.len() - 1) as f64
        / std::f64::consts::PI
}

/// Circular variance of nonzero Fourier-mode phases, in `[0, 1]`.
pub fn qft_phase_circular_variance(result: &QftResult) -> f64 {
    let mut cos_sum = 0.0;
    let mut sin_sum = 0.0;
    let mut count = 0usize;
    for amplitude in result.state.amplitudes() {
        if amplitude.norm_sqr() <= f64::EPSILON {
            continue;
        }
        let phase = amplitude.arg();
        cos_sum += phase.cos();
        sin_sum += phase.sin();
        count += 1;
    }
    if count == 0 {
        return 0.0;
    }
    let mean = (cos_sum / count as f64).hypot(sin_sum / count as f64);
    (1.0 - mean).clamp(0.0, 1.0)
}

/// Fidelity between the QFT output and the inverse-QFT output of a state.
pub fn qft_roundtrip_fidelity(
    targets: &[usize],
    qubits: usize,
    input: &StateVector,
) -> Result<f64, QftRoundTripError> {
    let forward = qft(targets, qubits)?;
    let mut state = input.clone();
    for gate in forward.circuit.gates() {
        apply_gate(gate, &mut state)?;
    }
    let inverse = inverse_qft(targets, qubits)?;
    for gate in inverse.circuit.gates() {
        apply_gate(gate, &mut state)?;
    }
    Ok(state.fidelity(input)?)
}

/// Difference between QFT gate density and inverse-QFT gate density.
pub fn qft_density_gap(result: &QftResult, inverse: &QftResult) -> f64 {
    (qft_gate_density(result) - qft_gate_density(inverse)).abs()
}

#[derive(Debug, thiserror::Error)]
pub enum QftRoundTripError {
    #[error(transparent)]
    Qft(#[from] QftError),
    #[error(transparent)]
    Simulator(#[from] SimulatorError),
    #[error(transparent)]
    State(#[from] crate::state::StateVectorError),
}

#[derive(Debug, thiserror::Error)]
pub enum QftError {
    #[error("qubit index {index} is out of range for {qubits} qubits")]
    QubitOutOfRange { index: usize, qubits: usize },
    #[error("duplicate QFT target qubit {index}")]
    DuplicateQubit { index: usize },
    #[error(transparent)]
    Simulator(#[from] SimulatorError),
}

/// Mean state fidelity of QFT outputs against the ideal basis and uniform references.
pub fn qft_gate_fidelity(result: &QftResult, qubits: usize) -> Result<f64, QftRoundTripError> {
    let ideal = StateVector::zero(qubits)?;
    let fidelity = result.state.fidelity(&ideal)?;
    Ok(fidelity.clamp(0.0, 1.0))
}

/// Resultant length of nonzero Fourier-mode phases in `[0, 1]`.
pub fn qft_phase_resultant_length(result: &QftResult) -> f64 {
    let mut cos_sum = 0.0;
    let mut sin_sum = 0.0;
    let mut count = 0usize;
    for amplitude in result.state.amplitudes() {
        if amplitude.norm_sqr() <= f64::EPSILON {
            continue;
        }
        let phase = amplitude.arg();
        cos_sum += phase.cos();
        sin_sum += phase.sin();
        count += 1;
    }
    if count == 0 {
        return 0.0;
    }
    (cos_sum / count as f64)
        .hypot(sin_sum / count as f64)
        .clamp(0.0, 1.0)
}
