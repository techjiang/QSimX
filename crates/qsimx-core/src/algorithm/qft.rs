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
