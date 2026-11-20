use crate::algorithm::{
    GroverAnalysisResult, GroverOracle, QftResult, grover_search, recommended_grover_iterations,
};
use crate::state::StateVector;

/// Grover iterations at which the marked-state probability reaches its first local maximum.
pub fn grover_peak_iteration(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> GroverAnalysisResult<usize> {
    let mut best_iteration = 0;
    let mut best_probability = -1.0;
    for iteration in 0..=max_iterations {
        let probability = grover_search(qubits, oracle, iteration)?.marked_probability;
        if probability > best_probability {
            best_probability = probability;
            best_iteration = iteration;
        }
    }
    Ok(best_iteration)
}

/// Difference between the recommended Grover iteration count and the observed peak.
pub fn grover_peak_delta(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> GroverAnalysisResult<isize> {
    Ok(recommended_grover_iterations(qubits) as isize
        - grover_peak_iteration(qubits, oracle, max_iterations)? as isize)
}

/// Largest basis-state probability imbalance between a QFT output and its inverse.
pub fn qft_probability_imbalance(forward: &QftResult, inverse: &QftResult) -> f64 {
    forward
        .state
        .probabilities()
        .iter()
        .zip(inverse.state.probabilities().iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max)
}

/// Whether applying a QFT circuit and then its inverse returns the original state.
pub fn qft_pair_is_reversible(
    forward: &QftResult,
    inverse: &QftResult,
    input: &StateVector,
) -> bool {
    use crate::circuit::apply_gate;

    let mut state = input.clone();
    for gate in forward.circuit.gates() {
        if apply_gate(gate, &mut state).is_err() {
            return false;
        }
    }
    for gate in inverse.circuit.gates() {
        if apply_gate(gate, &mut state).is_err() {
            return false;
        }
    }
    state
        .fidelity(input)
        .map(|fidelity| fidelity > 1.0 - 1e-9)
        .unwrap_or(false)
}

/// Smallest Grover iteration whose marked-state probability is at least `target`.
pub fn grover_first_success_iteration(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
    target: f64,
) -> GroverAnalysisResult<Option<usize>> {
    for iteration in 0..=max_iterations {
        if grover_search(qubits, oracle, iteration)?.marked_probability >= target {
            return Ok(Some(iteration));
        }
    }
    Ok(None)
}

/// Largest absolute probability change between adjacent Grover iterations.
pub fn grover_peak_step(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> GroverAnalysisResult<f64> {
    let curve = crate::algorithm::grover_success_curve(qubits, oracle, max_iterations)?;
    Ok(curve
        .windows(2)
        .map(|pair| (pair[1] - pair[0]).abs())
        .fold(0.0_f64, f64::max))
}

/// Time-averaged marked-state probability over all Grover iterations up to `max_iterations`.
pub fn grover_success_area(
    qubits: usize,
    oracle: &GroverOracle,
    max_iterations: usize,
) -> GroverAnalysisResult<f64> {
    let curve = crate::algorithm::grover_success_curve(qubits, oracle, max_iterations)?;
    Ok(curve.iter().sum::<f64>() / curve.len() as f64)
}
