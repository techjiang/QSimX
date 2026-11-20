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

#[derive(Clone, Debug)]
pub struct GroverResult {
    pub circuit: Circuit,
    pub state: StateVector,
    pub marked_probability: f64,
    pub iterations: usize,
}

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

    let mut circuit = Circuit::new(qubits);
    for q in 0..qubits {
        circuit.h(q);
    }

    for iteration in 0..iterations {
        match oracle {
            GroverOracle::Marked(_) => apply_phase_oracle(&mut circuit, qubits, marked)?,
            GroverOracle::Circuit(oracle_circuit) => append_gates(&mut circuit, oracle_circuit),
        }
        if iteration + 1 < iterations {
            apply_diffusion(&mut circuit, qubits)?;
        }
    }

    let state = StateVectorSimulator::new()
        .run_state(&circuit)
        .map_err(GroverError::Simulator)?;
    let marked_probability = match oracle {
        GroverOracle::Marked(_) => state.probability(marked),
        GroverOracle::Circuit(_) => state.probabilities().into_iter().fold(0.0_f64, f64::max),
    };

    Ok(GroverResult {
        circuit,
        state,
        marked_probability,
        iterations,
    })
}

fn apply_phase_oracle(
    circuit: &mut Circuit,
    qubits: usize,
    marked: usize,
) -> Result<(), GroverError> {
    for q in 0..qubits {
        if marked & (1 << q) == 0 {
            circuit.x(q);
        }
    }
    match qubits {
        1 => {
            circuit.z(0);
        }
        2 => {
            circuit.cz(0, 1);
        }
        _ => {
            for q in 0..qubits {
                circuit.h(q);
            }
            for q in 0..qubits {
                circuit.x(q);
            }
            for q in 1..qubits {
                circuit.h(q);
            }
            for control in 1..qubits {
                for target in 0..control {
                    circuit.cp(control, target, std::f64::consts::PI);
                }
            }
            for q in 1..qubits {
                circuit.h(q);
            }
            for q in 0..qubits {
                circuit.x(q);
            }
            for q in 0..qubits {
                circuit.h(q);
            }
        }
    }
    for q in 0..qubits {
        if marked & (1 << q) == 0 {
            circuit.x(q);
        }
    }
    Ok(())
}

fn apply_diffusion(circuit: &mut Circuit, qubits: usize) -> Result<(), GroverError> {
    for q in 0..qubits {
        circuit.h(q);
    }
    for q in 0..qubits {
        circuit.x(q);
    }
    phase_flip_all_ones(circuit, qubits)?;
    for q in 0..qubits {
        circuit.x(q);
    }
    for q in 0..qubits {
        circuit.h(q);
    }
    Ok(())
}

fn append_gates(circuit: &mut Circuit, source: &Circuit) {
    for gate in source.gates() {
        circuit.push(gate.clone());
    }
}

/// Flip the phase of the all-ones basis state.
fn phase_flip_all_ones(circuit: &mut Circuit, qubits: usize) -> Result<(), GroverError> {
    match qubits {
        0 => {}
        1 => {
            circuit.z(0);
        }
        2 => {
            circuit.cz(0, 1);
        }
        _ => {
            for q in 0..qubits {
                circuit.x(q);
            }
            for control in 0..qubits {
                for target in (control + 1)..qubits {
                    circuit.cp(control, target, std::f64::consts::PI);
                }
            }
            for q in 0..qubits {
                circuit.x(q);
            }
        }
    }
    Ok(())
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
