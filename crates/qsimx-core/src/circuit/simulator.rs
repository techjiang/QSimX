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
