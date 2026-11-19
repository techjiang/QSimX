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
}

impl Circuit {
    pub fn new(qubits: usize) -> Self {
        Self {
            qubits,
            gates: Vec::new(),
        }
    }

    pub fn qubits(&self) -> usize {
        self.qubits
    }

    pub fn gates(&self) -> &[Gate] {
        &self.gates
    }

    pub fn push(&mut self, gate: Gate) -> &mut Self {
        self.gates.push(gate);
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
        Gate::T(_) => " T ",
        Gate::RX(_, _) => "RX ",
        Gate::RY(_, _) => "RY ",
        Gate::RZ(_, _) => "RZ ",
        Gate::Swap(_, _) => "SW ",
        Gate::Toffoli(_, _, _) => "CCX",
        Gate::Custom1(_, _) => " U ",
    }
}
