use crate::algorithm::qft::{inverse_qft, qft};
use crate::circuit::{Circuit, SimulatorError, StateVectorSimulator};
use crate::state::StateVector;

#[derive(Clone, Debug)]
pub enum GroverOracle {
    Marked(usize),
    Circuit(Circuit),
}

impl GroverOracle {
    pub fn marked(state: usize) -> Self {
        Self::Marked(state)
    }

    pub fn circuit(circuit: Circuit) -> Self {
        Self::Circuit(circuit)
    }

}

#[derive(Clone, Debug, PartialEq)]
pub struct GroverResult {
    pub circuit: Circuit,
    pub state: StateVector,
    pub marked_probability: f64,
    pub iterations: usize,
}

/// Run Grover search with an exact amplitude amplification step.
pub fn grover_search(
    qubits: usize,
    oracle: &GroverOracle,
    iterations: usize,
) -> Result<GroverResult, GroverError> {
    if qubits == 0 {
        return Err(GroverError::EmptySearch);
    }
    let search_space = 1usize << qubits;
    let marked = match oracle {
        GroverOracle::Marked(state) => {
            if *state >= search_space {
                return Err(GroverError::MarkedOutOfRange {
                    index: *state,
                    qubits,
                });
            }
            *state
        }
        GroverOracle::Circuit(circuit) => {
            if circuit.qubits() != qubits {
                return Err(GroverError::OracleQubitMismatch {
                    expected: qubits,
                    actual: circuit.qubits(),
                });
            }
            0
        }
    };

    let targets: Vec<usize> = (0..qubits).collect();
    let mut circuit = Circuit::new(qubits);
    for q in 0..qubits {
        circuit.h(q);
    }

    for iteration in 0..iterations {
        match oracle {
            GroverOracle::Marked(_) => apply_phase_oracle(&mut circuit, qubits, marked),
            GroverOracle::Circuit(oracle_circuit) => {
                for gate in oracle_circuit.gates() {
                    circuit.push(gate.clone());
                }
            }
        }
        if iteration + 1 < iterations {
            apply_diffusion(&mut circuit, &targets)?;
        }
    }

    let state = StateVectorSimulator::new()
        .run_state(&circuit)
        .map_err(GroverError::Simulator)?;
    let marked_probability = match oracle {
        GroverOracle::Marked(_) => state.probability(marked),
        GroverOracle::Circuit(_) => state
            .probabilities()
            .into_iter()
            .fold(0.0_f64, f64::max),
    };

    Ok(GroverResult {
        circuit,
        state,
        marked_probability,
        iterations,
    })
}

fn apply_phase_oracle(circuit: &mut Circuit, qubits: usize, marked: usize) {
    for q in 0..qubits {
        if marked & (1 << q) == 0 {
            circuit.x(q);
        }
    }
    if qubits == 1 {
        circuit.z(0);
    } else {
        multi_controlled_z(circuit, &((0..qubits).collect::<Vec<_>>()));
    }
    for q in 0..qubits {
        if marked & (1 << q) == 0 {
            circuit.x(q);
        }
    }
}

fn apply_diffusion(circuit: &mut Circuit, targets: &[usize]) -> Result<(), GroverError> {
    for &q in targets {
        circuit.x(q);
    }
    multi_controlled_z(circuit, targets);
    for &q in targets {
        circuit.x(q);
    }
    let qubits = targets.iter().copied().max().unwrap_or(0) + 1;
    let forward = qft(targets, qubits).map_err(|_| GroverError::InvalidDiffuser)?;
    let inverse = inverse_qft(targets, qubits).map_err(|_| GroverError::InvalidDiffuser)?;
    append_gates(circuit, &forward.circuit);
    for gate in inverse.circuit.gates() {
        circuit.push(gate.clone());
    }
    Ok(())
}

fn append_gates(circuit: &mut Circuit, source: &Circuit) {
    for gate in source.gates() {
        circuit.push(gate.clone());
    }
}

fn multi_controlled_z(circuit: &mut Circuit, qubits: &[usize]) {
    match qubits.len() {
        0 => {}
        1 => {
            circuit.z(qubits[0]);
        }
        2 => {
            circuit.cz(qubits[0], qubits[1]);
        }
        _ => {
            let target = qubits[qubits.len() - 1];
            circuit.h(target);
            for pair in qubits.windows(2) {
                circuit.toffoli(pair[0], pair[1], target);
            }
            circuit.h(target);
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum GroverError {
    #[error("Grover search needs at least one qubit")]
    EmptySearch,
    #[error("marked state {index} is out of range for {qubits} qubits")]
    MarkedOutOfRange { index: usize, qubits: usize },
    #[error("oracle circuit has {actual} qubits, expected {expected}")]
    OracleQubitMismatch { expected: usize, actual: usize },
    #[error("failed to build Grover diffuser")]
    InvalidDiffuser,
    #[error(transparent)]
    Simulator(#[from] SimulatorError),
}
