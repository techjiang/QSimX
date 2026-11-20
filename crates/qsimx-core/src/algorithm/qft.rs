use std::f64::consts::PI;

use crate::circuit::{Circuit, SimulatorError, StateVectorSimulator};
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

#[derive(Debug, thiserror::Error)]
pub enum QftError {
    #[error("qubit index {index} is out of range for {qubits} qubits")]
    QubitOutOfRange { index: usize, qubits: usize },
    #[error("duplicate QFT target qubit {index}")]
    DuplicateQubit { index: usize },
    #[error(transparent)]
    Simulator(#[from] SimulatorError),
}
