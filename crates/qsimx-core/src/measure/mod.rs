use crate::state::{Basis, Complex, StateVector, StateVectorError};

pub fn expectation(
    state: &StateVector,
    basis: Basis,
    qubit: usize,
) -> Result<f64, StateVectorError> {
    let matrix = match basis {
        Basis::Z => [
            [Complex::ONE, Complex::ZERO],
            [Complex::ZERO, Complex::new(-1.0, 0.0)],
        ],
        Basis::X => [[Complex::ZERO, Complex::ONE], [Complex::ONE, Complex::ZERO]],
        Basis::Y => [
            [Complex::ZERO, Complex::new(0.0, -1.0)],
            [Complex::new(0.0, 1.0), Complex::ZERO],
        ],
    };
    state.check_qubit(qubit)?;
    let mut measured = state.clone();
    measured.apply_1q(qubit, matrix)?;
    let expectation: Complex = state
        .amplitudes()
        .iter()
        .zip(measured.amplitudes())
        .map(|(&a, &b)| a.conj() * b)
        .fold(Complex::ZERO, |mut acc, value| {
            acc += value;
            acc
        });
    Ok(expectation.re)
}

pub fn computational_probability(state: &StateVector, basis_state: usize) -> f64 {
    state.probability(basis_state)
}
