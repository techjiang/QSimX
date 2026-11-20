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

/// Recommended Grover iteration count for a single marked state.
///
/// Clamped to `sqrt(N)` so small search spaces do not overshoot.
pub type GroverAnalysisResult<T> = Result<T, GroverError>;

pub fn recommended_grover_iterations(qubits: usize) -> usize {
    if qubits == 0 {
        return 1;
    }
    grover_iterations(qubits, 1)
}

/// Iteration count for `solutions` marked states in a `qubits`-bit search space.
/// Uses `floor(pi/4 * sqrt(N/M))` and never returns zero.
pub fn grover_iterations(qubits: usize, solutions: usize) -> usize {
    if qubits == 0 || solutions == 0 {
        return 1;
    }
    let space = 1usize << qubits;
    let ratio = (space / solutions.max(1)).max(1) as f64;
    ((std::f64::consts::PI / 4.0) * ratio.sqrt())
        .floor()
        .max(1.0) as usize
}

/// Grover success probability after each iteration, starting from 0 iterations.
pub fn grover_success_curve(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> Result<Vec<f64>, GroverError> {
    (0..=max_iterations)
        .map(|iterations| {
            grover_search(qubits, oracle, iterations).map(|result| result.marked_probability)
        })
        .collect()
}

/// Total absolute probability movement across all Grover iterations.
pub fn grover_success_path_length(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> Result<f64, GroverError> {
    let curve = grover_success_curve(qubits, oracle, max_iterations)?;
    Ok(curve.windows(2).map(|pair| (pair[1] - pair[0]).abs()).sum())
}

/// Amplification gain over the uniform prior for a single marked state.
pub fn grover_amplification_gain(qubits: usize, probability: f64) -> f64 {
    if qubits == 0 {
        return 1.0;
    }
    if qubits >= usize::BITS as usize {
        return probability * (usize::MAX as f64);
    }
    probability * (1usize << qubits) as f64
}

/// Cumulative Grover probability after each iteration, including the initial state.
pub fn grover_success_cumulative(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> Result<Vec<f64>, GroverError> {
    let curve = grover_success_curve(qubits, oracle, max_iterations)?;
    let mut cumulative = 0.0;
    Ok(curve
        .into_iter()
        .map(|probability| {
            cumulative += probability;
            cumulative
        })
        .collect())
}

/// Probability gain from the first Grover iteration to the recommended count.
pub fn grover_recommended_gain(qubits: usize, oracle: &GroverOracle) -> Result<f64, GroverError> {
    let initial = grover_search(qubits, oracle, 0)?.marked_probability;
    let final_probability =
        grover_search(qubits, oracle, recommended_grover_iterations(qubits))?.marked_probability;
    Ok(final_probability - initial)
}

/// Skewness of the Grover success probabilities across the evaluated iterations.
pub fn grover_success_skewness(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> Result<f64, GroverError> {
    let curve = grover_success_curve(qubits, oracle, max_iterations)?;
    let mean = curve.iter().sum::<f64>() / curve.len() as f64;
    let moment2 = curve.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / curve.len() as f64;
    if moment2 <= f64::EPSILON {
        return Ok(0.0);
    }
    let moment3 = curve.iter().map(|p| (p - mean).powi(3)).sum::<f64>() / curve.len() as f64;
    Ok(moment3 / moment2.powf(1.5))
}

/// Excess kurtosis of the Grover success probabilities across the evaluated iterations.
pub fn grover_success_kurtosis(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> Result<f64, GroverError> {
    let curve = grover_success_curve(qubits, oracle, max_iterations)?;
    let mean = curve.iter().sum::<f64>() / curve.len() as f64;
    let moment2 = curve.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / curve.len() as f64;
    if moment2 <= f64::EPSILON {
        return Ok(0.0);
    }
    let moment4 = curve.iter().map(|p| (p - mean).powi(4)).sum::<f64>() / curve.len() as f64;
    Ok(moment4 / moment2.powi(2) - 3.0)
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
