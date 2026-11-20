use std::f64::consts::{FRAC_1_SQRT_2, PI};

use crate::state::{Complex, StateVector, StateVectorError};

pub const HADAMARD: [[Complex; 2]; 2] = [
    [
        Complex::new(FRAC_1_SQRT_2, 0.0),
        Complex::new(FRAC_1_SQRT_2, 0.0),
    ],
    [
        Complex::new(FRAC_1_SQRT_2, 0.0),
        Complex::new(-FRAC_1_SQRT_2, 0.0),
    ],
];
pub const PAULI_X: [[Complex; 2]; 2] =
    [[Complex::ZERO, Complex::ONE], [Complex::ONE, Complex::ZERO]];
pub const PAULI_Y: [[Complex; 2]; 2] = [
    [Complex::ZERO, Complex::new(0.0, -1.0)],
    [Complex::new(0.0, 1.0), Complex::ZERO],
];
pub const PAULI_Z: [[Complex; 2]; 2] = [
    [Complex::ONE, Complex::ZERO],
    [Complex::ZERO, Complex::new(-1.0, 0.0)],
];
pub const S_GATE: [[Complex; 2]; 2] = [
    [Complex::ONE, Complex::ZERO],
    [Complex::ZERO, Complex::new(0.0, 1.0)],
];
pub const T_GATE: [[Complex; 2]; 2] = [
    [Complex::ONE, Complex::ZERO],
    [Complex::ZERO, Complex::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2)],
];

#[derive(Clone, Debug)]
pub enum Gate {
    H(usize),
    X(usize),
    Y(usize),
    Z(usize),
    S(usize),
    T(usize),
    RX(usize, f64),
    RY(usize, f64),
    RZ(usize, f64),
    CX(usize, usize),
    CZ(usize, usize),
    CP(usize, usize, f64),
    Swap(usize, usize),
    Toffoli(usize, usize, usize),
    Custom1(usize, [[Complex; 2]; 2]),
}

pub fn rx_matrix(theta: f64) -> [[Complex; 2]; 2] {
    let c = (theta / 2.0).cos();
    let s = (theta / 2.0).sin();
    [
        [Complex::new(c, 0.0), Complex::new(0.0, -s)],
        [Complex::new(0.0, -s), Complex::new(c, 0.0)],
    ]
}

pub fn ry_matrix(theta: f64) -> [[Complex; 2]; 2] {
    let c = (theta / 2.0).cos();
    let s = (theta / 2.0).sin();
    [
        [Complex::new(c, 0.0), Complex::new(-s, 0.0)],
        [Complex::new(s, 0.0), Complex::new(c, 0.0)],
    ]
}

pub fn rz_matrix(theta: f64) -> [[Complex; 2]; 2] {
    [
        [Complex::from_polar(1.0, -theta / 2.0), Complex::ZERO],
        [Complex::ZERO, Complex::from_polar(1.0, theta / 2.0)],
    ]
}

impl Gate {
    /// Stable human-readable name for this gate.
    pub fn name(&self) -> &'static str {
        match self {
            Gate::H(_) => "H",
            Gate::X(_) => "X",
            Gate::Y(_) => "Y",
            Gate::Z(_) => "Z",
            Gate::S(_) => "S",
            Gate::T(_) => "T",
            Gate::RX(_, _) => "RX",
            Gate::RY(_, _) => "RY",
            Gate::RZ(_, _) => "RZ",
            Gate::CX(_, _) => "CX",
            Gate::CZ(_, _) => "CZ",
            Gate::CP(_, _, _) => "CP",
            Gate::Swap(_, _) => "SWAP",
            Gate::Toffoli(_, _, _) => "TOFFOLI",
            Gate::Custom1(_, _) => "U",
        }
    }

    /// QASM gate name accepted by qelib1.inc.
    pub fn qasm_name(&self) -> &'static str {
        match self {
            Gate::H(_) => "h",
            Gate::X(_) => "x",
            Gate::Y(_) => "y",
            Gate::Z(_) => "z",
            Gate::S(_) => "s",
            Gate::T(_) => "t",
            Gate::RX(_, _) => "rx",
            Gate::RY(_, _) => "ry",
            Gate::RZ(_, _) => "rz",
            Gate::CX(_, _) => "cx",
            Gate::CZ(_, _) => "cz",
            Gate::CP(_, _, _) => "cp",
            Gate::Swap(_, _) => "swap",
            Gate::Toffoli(_, _, _) => "ccx",
            Gate::Custom1(_, _) => "u3",
        }
    }

    pub fn qubits(&self) -> Vec<usize> {
        match *self {
            Gate::H(q)
            | Gate::X(q)
            | Gate::Y(q)
            | Gate::Z(q)
            | Gate::S(q)
            | Gate::T(q)
            | Gate::RX(q, _)
            | Gate::RY(q, _)
            | Gate::RZ(q, _)
            | Gate::Custom1(q, _) => vec![q],
            Gate::CX(c, t) | Gate::CZ(c, t) | Gate::Swap(c, t) => vec![c, t],
            Gate::CP(c, t, _) => vec![c, t],
            Gate::Toffoli(c1, c2, t) => vec![c1, c2, t],
        }
    }

    pub fn qasm(&self) -> String {
        let angle = |value: f64| format!("{value:.8}");
        match *self {
            Gate::H(q) => format!("h q[{q}];"),
            Gate::X(q) => format!("x q[{q}];"),
            Gate::Y(q) => format!("y q[{q}];"),
            Gate::Z(q) => format!("z q[{q}];"),
            Gate::S(q) => format!("s q[{q}];"),
            Gate::T(q) => format!("t q[{q}];"),
            Gate::RX(q, theta) => format!("rx({}) q[{q}];", angle(theta)),
            Gate::RY(q, theta) => format!("ry({}) q[{q}];", angle(theta)),
            Gate::RZ(q, theta) => format!("rz({}) q[{q}];", angle(theta)),
            Gate::CX(c, t) => format!("cx q[{c}],q[{t}];"),
            Gate::CZ(c, t) => format!("cz q[{c}],q[{t}];"),
            Gate::CP(c, t, theta) => format!("cp({}) q[{c}],q[{t}];", angle(theta)),
            Gate::Swap(a, b) => format!("swap q[{a}],q[{b}];"),
            Gate::Toffoli(c1, c2, t) => format!("ccx q[{c1}],q[{c2}],q[{t}];"),
            Gate::Custom1(q, _) => format!("u3(0,0,0) q[{q}]; // custom gate"),
        }
    }

    pub fn apply(&self, state: &mut StateVector) -> Result<(), StateVectorError> {
        match *self {
            Gate::H(q) => state.apply_1q(q, HADAMARD),
            Gate::X(q) => state.apply_1q(q, PAULI_X),
            Gate::Y(q) => state.apply_1q(q, PAULI_Y),
            Gate::Z(q) => state.apply_1q(q, PAULI_Z),
            Gate::S(q) => state.apply_1q(q, S_GATE),
            Gate::T(q) => state.apply_1q(q, T_GATE),
            Gate::RX(q, theta) => state.apply_1q(q, rx_matrix(theta)),
            Gate::RY(q, theta) => state.apply_1q(q, ry_matrix(theta)),
            Gate::RZ(q, theta) => state.apply_1q(q, rz_matrix(theta)),
            Gate::CX(c, t) => state.apply_2q(c, t, PAULI_X),
            Gate::CZ(c, t) => state.apply_2q(c, t, PAULI_Z),
            Gate::CP(c, t, theta) => state.apply_cp(c, t, theta),
            Gate::Custom1(q, matrix) => state.apply_1q(q, matrix),
            Gate::Swap(a, b) => {
                state.apply_2q(a, b, PAULI_X)?;
                state.apply_2q(b, a, PAULI_X)?;
                state.apply_2q(a, b, PAULI_X)
            }
            Gate::Toffoli(c1, c2, t) => {
                state.check_qubit(c1)?;
                state.check_qubit(c2)?;
                state.check_qubit(t)?;
                if c1 == c2 || c1 == t || c2 == t {
                    return Err(StateVectorError::SameQubit);
                }
                let c1_bit = 1 << c1;
                let c2_bit = 1 << c2;
                let t_bit = 1 << t;
                for index in 0..state.amplitudes().len() {
                    if index & c1_bit == 0 || index & c2_bit == 0 || index & t_bit != 0 {
                        continue;
                    }
                    let other = index ^ t_bit;
                    let a = state.amplitudes()[index];
                    state.amplitudes_mut()[index] = state.amplitudes()[other];
                    state.amplitudes_mut()[other] = a;
                }
                Ok(())
            }
        }
    }
}

pub(crate) fn _pi() -> f64 {
    PI
}
