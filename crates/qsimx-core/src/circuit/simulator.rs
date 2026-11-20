use std::collections::HashMap;

use crate::circuit::Circuit;
use crate::gate::Gate;
use crate::state::{Complex, StateVector, StateVectorError};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SimulatorError {
    #[error(transparent)]
    State(#[from] StateVectorError),
}

#[derive(Clone, Debug)]
pub struct SimulationResult {
    state: StateVector,
    counts: HashMap<usize, usize>,
    shots: usize,
}

impl SimulationResult {
    pub fn probabilities(&self) -> Vec<f64> {
        self.state.probabilities()
    }
    pub fn state_vector(&self) -> &[Complex] {
        self.state.amplitudes()
    }
    pub fn state(&self) -> &StateVector {
        &self.state
    }
    pub fn counts(&self) -> &HashMap<usize, usize> {
        &self.counts
    }
    pub fn shots(&self) -> usize {
        self.shots
    }

    pub fn bloch_vector(&self, qubit: usize) -> Result<[f64; 3], StateVectorError> {
        self.state.bloch_vector(qubit)
    }

    pub fn sorted_counts(&self) -> Vec<(usize, usize)> {
        let mut counts: Vec<_> = self.counts.iter().map(|(&k, &v)| (k, v)).collect();
        counts.sort_by_key(|&(state, _)| state);
        counts
    }

    /// Shannon entropy of the measurement distribution, in bits.
    pub fn shannon_entropy(&self) -> f64 {
        if self.shots == 0 {
            return 0.0;
        }
        -self
            .counts
            .values()
            .map(|&count| {
                let p = count as f64 / self.shots as f64;
                p * p.log2()
            })
            .sum::<f64>()
    }

    /// Squared overlap between this result and the simulated state of another circuit.
    pub fn state_fidelity(&self, other_state: &crate::state::StateVector) -> Result<f64, SimulatorError> {
        Ok(self.state.fidelity(other_state)?)
    }

    /// Empirical probability of one measurement outcome.
    pub fn probability_by_state(&self, state: usize) -> f64 {
        if self.shots == 0 {
            return 0.0;
        }
        self.counts.get(&state).copied().unwrap_or(0) as f64 / self.shots as f64
    }

    /// Renyi entropy with alpha = 2, also called collision entropy.
    pub fn min_entropy(&self) -> f64 {
        if self.shots == 0 {
            return 0.0;
        }
        -self
            .counts
            .values()
            .map(|&count| {
                let p = count as f64 / self.shots as f64;
                p * p
            })
            .sum::<f64>()
            .log2()
    }

    /// Highest-frequency measurement outcomes.
    pub fn top_counts(&self, limit: usize) -> Vec<(usize, usize)> {
        let mut counts: Vec<_> = self.counts.iter().map(|(&k, &v)| (k, v)).collect();
        counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        counts.truncate(limit);
        counts
    }

    /// Most likely measurement outcome and its probability estimate.
    pub fn most_likely(&self) -> Option<(usize, f64)> {
        if self.shots == 0 {
            return None;
        }
        self.top_counts(1)
            .first()
            .map(|&(state, count)| (state, count as f64 / self.shots as f64))
    }

    /// Render the measurement histogram as compact JSON.
    pub fn histogram_json(&self, width: usize) -> String {
        let entries: Vec<String> = self
            .sorted_counts()
            .iter()
            .map(|(state, count)| {
                format!(
                    "{{\"state\":\"{state:0width$b}\",\"count\":{count}}}"
                )
            })
            .collect();
        format!("{{\"shots\":{},\"counts\":[{}]}}", self.shots, entries.join(","))
    }
}

#[derive(Debug, Default)]
pub struct StateVectorSimulator;

impl StateVectorSimulator {
    pub fn new() -> Self {
        Self
    }

    pub fn run(&self, circuit: &Circuit, shots: usize) -> Result<SimulationResult, SimulatorError> {
        let mut state = StateVector::zero(circuit.qubits())?;
        self.apply(circuit, &mut state)?;
        let mut counts = HashMap::new();
        for sample in state.sample(shots, 0x5158_494d) {
            *counts.entry(sample).or_insert(0) += 1;
        }
        Ok(SimulationResult {
            state,
            counts,
            shots,
        })
    }

    pub fn run_state(&self, circuit: &Circuit) -> Result<StateVector, SimulatorError> {
        let mut state = StateVector::zero(circuit.qubits())?;
        self.apply(circuit, &mut state)?;
        Ok(state)
    }

    fn apply(&self, circuit: &Circuit, state: &mut StateVector) -> Result<(), SimulatorError> {
        for gate in circuit.gates() {
            gate.apply(state)?;
        }
        Ok(())
    }
}

pub fn apply_gate(gate: &Gate, state: &mut StateVector) -> Result<(), SimulatorError> {
    gate.apply(state)?;
    Ok(())
}
