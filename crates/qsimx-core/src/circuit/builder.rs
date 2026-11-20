use std::fmt::Write;

use crate::gate::Gate;
use crate::state::Complex;
use thiserror::Error;

#[derive(Debug, Error, Eq, PartialEq)]
pub enum QasmError {
    #[error("circuit has no qubits")]
    EmptyCircuit,
    #[error("qubit index {index} is out of range for {qubits} qubits")]
    QubitOutOfRange { index: usize, qubits: usize },
}

#[derive(Clone, Debug)]
pub struct Circuit {
    qubits: usize,
    gates: Vec<Gate>,
    layers: Vec<usize>,
}

impl Circuit {
    pub fn new(qubits: usize) -> Self {
        Self {
            qubits,
            gates: Vec::new(),
            layers: vec![0; qubits],
        }
    }

    pub fn qubits(&self) -> usize {
        self.qubits
    }

    pub fn gates(&self) -> &[Gate] {
        &self.gates
    }

    pub fn push(&mut self, gate: Gate) -> &mut Self {
        let gate_qubits = gate.qubits();
        let gate_depth = gate_qubits
            .iter()
            .filter_map(|&q| self.layers.get(q).copied())
            .max()
            .unwrap_or(0)
            + 1;
        for q in gate_qubits {
            if q < self.layers.len() {
                self.layers[q] = gate_depth;
            }
        }
        self.gates.push(gate);
        self
    }

    /// Qubits that have appeared in at least one gate so far, in index order.
    pub fn distinct_touching_qubits(&self) -> Vec<usize> {
        let mut seen = vec![false; self.qubits];
        for gate in &self.gates {
            for qubit in gate.qubits() {
                if qubit < seen.len() {
                    seen[qubit] = true;
                }
            }
        }
        seen.iter()
            .enumerate()
            .filter_map(|(q, &used)| used.then_some(q))
            .collect()
    }

    /// Append a gate after every currently active qubit has been touched.
    pub fn add_layered(&mut self, gate: Gate) -> &mut Self {
        let target_layer = self
            .distinct_touching_qubits()
            .into_iter()
            .max()
            .map_or(1, |q| self.layers.get(q).copied().unwrap_or(0) + 1);
        self.push(gate);
        if let Some(layer) = self.layers.get_mut(0) {
            *layer = (*layer).max(target_layer);
        }
        self
    }

    pub fn h(&mut self, q: usize) -> &mut Self {
        self.push(Gate::H(q))
    }
    pub fn x(&mut self, q: usize) -> &mut Self {
        self.push(Gate::X(q))
    }
    pub fn y(&mut self, q: usize) -> &mut Self {
        self.push(Gate::Y(q))
    }
    pub fn z(&mut self, q: usize) -> &mut Self {
        self.push(Gate::Z(q))
    }
    pub fn s(&mut self, q: usize) -> &mut Self {
        self.push(Gate::S(q))
    }
    pub fn sdg(&mut self, q: usize) -> &mut Self {
        self.push(Gate::SDG(q))
    }
    pub fn t(&mut self, q: usize) -> &mut Self {
        self.push(Gate::T(q))
    }
    pub fn rx(&mut self, q: usize, theta: f64) -> &mut Self {
        self.push(Gate::RX(q, theta))
    }
    pub fn ry(&mut self, q: usize, theta: f64) -> &mut Self {
        self.push(Gate::RY(q, theta))
    }
    pub fn rz(&mut self, q: usize, theta: f64) -> &mut Self {
        self.push(Gate::RZ(q, theta))
    }
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.push(Gate::CX(control, target))
    }
    pub fn cz(&mut self, control: usize, target: usize) -> &mut Self {
        self.push(Gate::CZ(control, target))
    }
    pub fn cp(&mut self, control: usize, target: usize, theta: f64) -> &mut Self {
        self.push(Gate::CP(control, target, theta))
    }
    pub fn swap(&mut self, a: usize, b: usize) -> &mut Self {
        self.push(Gate::Swap(a, b))
    }
    pub fn toffoli(&mut self, c1: usize, c2: usize, target: usize) -> &mut Self {
        self.push(Gate::Toffoli(c1, c2, target))
    }
    pub fn custom(&mut self, q: usize, matrix: [[Complex; 2]; 2]) -> &mut Self {
        self.push(Gate::Custom1(q, matrix))
    }

    pub fn to_qasm(&self) -> Result<String, QasmError> {
        if self.qubits == 0 {
            return Err(QasmError::EmptyCircuit);
        }
        for gate in &self.gates {
            for &q in &gate.qubits() {
                if q >= self.qubits {
                    return Err(QasmError::QubitOutOfRange {
                        index: q,
                        qubits: self.qubits,
                    });
                }
            }
        }
        let mut output = String::new();
        let _ = writeln!(output, "OPENQASM 2.0;");
        let _ = writeln!(output, "include \"qelib1.inc\";");
        let _ = writeln!(output, "qreg q[{}];", self.qubits);
        let _ = writeln!(output, "creg c[{}];", self.qubits);
        for gate in &self.gates {
            let _ = writeln!(output, "{}", gate.qasm());
        }
        for q in 0..self.qubits {
            let _ = writeln!(output, "measure q[{q}] -> c[{q}];");
        }
        Ok(output)
    }

    /// Number of sequential gate layers: each gate occupies the layers of every qubit it touches.
    pub fn depth(&self) -> usize {
        self.layers.iter().copied().max().unwrap_or(0)
    }

    /// Count gates grouped by their stable names.
    pub fn gate_counts(&self) -> std::collections::BTreeMap<&'static str, usize> {
        let mut counts = std::collections::BTreeMap::new();
        for gate in &self.gates {
            *counts.entry(gate.name()).or_insert(0) += 1;
        }
        counts
    }

    /// Summary tuple used by CLI and tests: (qubits, gates, depth).
    pub fn summary(&self) -> (usize, usize, usize) {
        (self.qubits, self.gates.len(), self.depth())
    }

    /// Dominant amplitude probability of a quantum state (max |amplitude|^2).
    pub fn state_purity(&self, state: &crate::state::StateVector) -> f64 {
        let mut max_probability = 0.0_f64;
        for amplitude in state.amplitudes() {
            max_probability = max_probability.max(amplitude.norm_sqr());
        }
        max_probability
    }

    /// Build a Bell state on two qubits.
    pub fn bell_pair(&mut self, a: usize, b: usize) -> &mut Self {
        self.h(a).cnot(a, b)
    }

    /// Prepare an equal superposition over `targets`, like the existing W helper.
    pub fn w_state(&mut self, targets: &[usize]) -> Result<&mut Self, QasmError> {
        if targets.is_empty() {
            return Err(QasmError::EmptyCircuit);
        }
        for &target in targets {
            self.check_target(target)?;
        }
        let first = targets[0];
        self.h(first);
        for &target in targets.iter().skip(1) {
            self.cnot(first, target);
        }
        Ok(self)
    }

    fn check_target(&self, target: usize) -> Result<(), QasmError> {
        if target >= self.qubits {
            Err(QasmError::QubitOutOfRange {
                index: target,
                qubits: self.qubits,
            })
        } else {
            Ok(())
        }
    }

    /// Maximum number of gates applied to any qubit.
    /// Minimum gate count across qubits that are touched by at least one gate.
    pub fn min_gates_per_touched_qubit(&self) -> usize {
        let touched = self.distinct_touching_qubits();
        touched
            .iter()
            .filter_map(|&qubit| self.layers.get(qubit).copied())
            .min()
            .unwrap_or(0)
    }

    /// Number of qubits touched by at least one gate.
    pub fn active_qubit_count(&self) -> usize {
        self.distinct_touching_qubits().len()
    }

    /// Number of gates whose definition depends on a continuous angle.
    pub fn parameterized_gate_count(&self) -> usize {
        self.gates
            .iter()
            .filter(|gate| {
                matches!(
                    gate,
                    Gate::RX(_, _) | Gate::RY(_, _) | Gate::RZ(_, _) | Gate::CP(_, _, _)
                )
            })
            .count()
    }

    pub fn max_gates_per_qubit(&self) -> usize {
        let mut counts = vec![0usize; self.qubits];
        for gate in &self.gates {
            let indices = match gate {
                Gate::H(q)
                | Gate::X(q)
                | Gate::Y(q)
                | Gate::Z(q)
                | Gate::S(q)
                | Gate::SDG(q)
                | Gate::T(q)
                | Gate::RX(q, _)
                | Gate::RY(q, _)
                | Gate::RZ(q, _)
                | Gate::Custom1(q, _) => vec![*q],
                Gate::CX(a, b) | Gate::CZ(a, b) | Gate::CP(a, b, _) | Gate::Swap(a, b) => {
                    vec![*a, *b]
                }
                Gate::Toffoli(a, b, c) => vec![*a, *b, *c],
            };
            for index in indices {
                if index < counts.len() {
                    counts[index] += 1;
                }
            }
        }
        counts.into_iter().max().unwrap_or(0)
    }

    pub fn diagram(&self) -> String {
        (0..self.qubits)
            .map(|q| {
                let mut line = format!("q{q}: ");
                for gate in &self.gates {
                    let symbol = if gate.qubits().contains(&q) {
                        gate_symbol(gate)
                    } else {
                        "───"
                    };
                    let _ = write!(line, "{symbol:>3}");
                }
                line
            })
            .collect::<Vec<_>>()
            .join("\n")
    }
}

fn gate_symbol(gate: &Gate) -> &'static str {
    match gate {
        Gate::H(_) => " H ",
        Gate::X(_) | Gate::CX(_, _) => " X ",
        Gate::Y(_) => " Y ",
        Gate::Z(_) | Gate::CZ(_, _) => " Z ",
        Gate::S(_) => " S ",
        Gate::SDG(_) => "SDG",
        Gate::T(_) => " T ",
        Gate::RX(_, _) => "RX ",
        Gate::RY(_, _) => "RY ",
        Gate::RZ(_, _) => "RZ ",
        Gate::CP(_, _, _) => "CP ",
        Gate::Swap(_, _) => "SW ",
        Gate::Toffoli(_, _, _) => "CCX",
        Gate::Custom1(_, _) => " U ",
    }
}
