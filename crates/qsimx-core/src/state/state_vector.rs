use std::ops::{Add, AddAssign, Mul, MulAssign};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };

    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn from_polar(radius: f64, angle: f64) -> Self {
        Self::new(radius * angle.cos(), radius * angle.sin())
    }

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn abs(self) -> f64 {
        self.norm_sqr().sqrt()
    }

    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{:.6}+{:.6}i", self.re, self.im)
        } else {
            write!(f, "{:.6}{:.6}i", self.re, self.im)
        }
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Basis {
    Z,
    X,
    Y,
}

#[derive(Debug, Error, PartialEq)]
pub enum StateVectorError {
    #[error("quantum state must contain at least one qubit")]
    EmptyState,
    #[error("qubit index {index} is out of range for {qubits} qubits")]
    QubitOutOfRange { index: usize, qubits: usize },
    #[error("state vector has wrong length: expected {expected}, got {actual}")]
    WrongLength { expected: usize, actual: usize },
    #[error("state vector is not normalized: norm squared is {norm_sqr}")]
    NotNormalized { norm_sqr: f64 },
    #[error("control and target qubits must be distinct")]
    SameQubit,
}

#[derive(Clone, Debug)]
pub struct StateVector {
    qubits: usize,
    amplitudes: Vec<Complex>,
}

impl StateVector {
    pub fn zero(qubits: usize) -> Result<Self, StateVectorError> {
        if qubits == 0 {
            return Err(StateVectorError::EmptyState);
        }
        let mut amplitudes = vec![Complex::ZERO; 1 << qubits];
        amplitudes[0] = Complex::ONE;
        Ok(Self { qubits, amplitudes })
    }

    pub fn from_amplitudes(
        qubits: usize,
        amplitudes: Vec<Complex>,
    ) -> Result<Self, StateVectorError> {
        if qubits == 0 {
            return Err(StateVectorError::EmptyState);
        }
        let expected = 1usize
            .checked_shl(qubits as u32)
            .ok_or(StateVectorError::WrongLength {
                expected: usize::MAX,
                actual: amplitudes.len(),
            })?;
        if amplitudes.len() != expected {
            return Err(StateVectorError::WrongLength {
                expected,
                actual: amplitudes.len(),
            });
        }
        let state = Self { qubits, amplitudes };
        if !state.is_normalized() {
            return Err(StateVectorError::NotNormalized {
                norm_sqr: state.norm_sqr(),
            });
        }
        Ok(state)
    }

    pub fn qubits(&self) -> usize {
        self.qubits
    }

    pub fn amplitudes(&self) -> &[Complex] {
        &self.amplitudes
    }

    pub fn amplitudes_mut(&mut self) -> &mut [Complex] {
        &mut self.amplitudes
    }

    pub fn norm_sqr(&self) -> f64 {
        self.amplitudes.iter().map(|a| a.norm_sqr()).sum()
    }

    /// Render one amplitude in `a+bi` form with fixed precision.
    /// Number of computational basis states with probability above `threshold`.
    pub fn significant_basis_count(&self, threshold: f64) -> usize {
        self.amplitudes
            .iter()
            .filter(|amplitude| amplitude.norm_sqr() > threshold)
            .count()
    }

    pub fn amplitude_string(&self, basis_state: usize, precision: usize) -> String {
        match self.amplitudes.get(basis_state) {
            Some(a) => format!("{a:.precision$}"),
            None => format!("{:.precision$}", Complex::ZERO),
        }
    }

    /// Squared overlap |<self|other>|^2.
    pub fn fidelity(&self, other: &StateVector) -> Result<f64, StateVectorError> {
        if self.qubits != other.qubits {
            return Err(StateVectorError::WrongLength {
                expected: self.amplitudes.len(),
                actual: other.amplitudes.len(),
            });
        }
        let overlap: Complex = self
            .amplitudes
            .iter()
            .zip(&other.amplitudes)
            .map(|(&a, &b)| a.conj() * b)
            .fold(Complex::ZERO, |mut acc, value| {
                acc += value;
                acc
            });
        Ok(overlap.norm_sqr())
    }

    pub fn probability(&self, basis_state: usize) -> f64 {
        self.amplitudes
            .get(basis_state)
            .map(|a| a.norm_sqr())
            .unwrap_or(0.0)
    }

    pub fn probabilities(&self) -> Vec<f64> {
        self.amplitudes.iter().map(|a| a.norm_sqr()).collect()
    }

    /// Quantum phase of a basis-state amplitude, in radians.
    pub fn phase(&self, basis_state: usize) -> f64 {
        self.amplitudes
            .get(basis_state)
            .map(|amplitude| amplitude.arg())
            .unwrap_or(0.0)
    }

    /// Probability of the |0> outcome for one measured qubit.
    pub fn zero_probability(&self, qubit: usize) -> Result<f64, StateVectorError> {
        self.check_qubit(qubit)?;
        Ok(self
            .amplitudes
            .iter()
            .enumerate()
            .filter(|(index, _)| index & (1 << qubit) == 0)
            .map(|(_, amplitude)| amplitude.norm_sqr())
            .sum())
    }

    /// Probability of the |1> outcome for one measured qubit.
    pub fn one_probability(&self, qubit: usize) -> Result<f64, StateVectorError> {
        Ok(1.0 - self.zero_probability(qubit)?)
    }

    /// Von Neumann entropy of one qubit, in bits.
    pub fn entanglement_entropy(&self, qubit: usize) -> Result<f64, StateVectorError> {
        let rho = self.reduced_density_matrix(qubit)?;
        let trace = rho[0][0].re + rho[1][1].re;
        if trace <= f64::EPSILON {
            return Ok(0.0);
        }
        let det = rho[0][0].re * rho[1][1].re - rho[0][1].norm_sqr();
        let discriminant = (trace * trace - 4.0 * det).max(0.0).sqrt();
        let lambda1 = ((trace + discriminant) / 2.0).clamp(0.0, 1.0);
        let lambda2 = ((trace - discriminant) / 2.0).clamp(0.0, 1.0);
        let term = |lambda: f64| {
            if lambda > f64::EPSILON {
                -lambda * lambda.log2()
            } else {
                0.0
            }
        };
        Ok(term(lambda1) + term(lambda2))
    }

    /// Purity of one qubit's reduced density matrix.
    pub fn purity(&self, qubit: usize) -> Result<f64, StateVectorError> {
        let rho = self.reduced_density_matrix(qubit)?;
        Ok(rho[0][0].norm_sqr() + 2.0 * rho[0][1].norm_sqr() + rho[1][1].norm_sqr())
    }

    /// Total variation distance between two probability distributions.
    pub fn total_variation_distance(&self, other: &StateVector) -> Result<f64, StateVectorError> {
        if self.qubits != other.qubits {
            return Err(StateVectorError::WrongLength {
                expected: self.amplitudes.len(),
                actual: other.amplitudes.len(),
            });
        }
        Ok(self
            .amplitudes
            .iter()
            .zip(&other.amplitudes)
            .map(|(&a, &b)| (a.norm_sqr() - b.norm_sqr()).abs())
            .sum::<f64>()
            / 2.0)
    }

    /// Difference between the largest and smallest basis-state probability.
    pub fn probability_span(&self) -> f64 {
        let mut min = f64::INFINITY;
        let mut max = 0.0_f64;
        for amplitude in &self.amplitudes {
            let probability = amplitude.norm_sqr();
            min = min.min(probability);
            max = max.max(probability);
        }
        if self.amplitudes.is_empty() {
            0.0
        } else {
            max - min
        }
    }

    /// Bitwise XOR parity of a basis state, useful for parity-check summaries.
    pub fn basis_parity(&self, basis_state: usize) -> Result<u8, StateVectorError> {
        if basis_state >= self.amplitudes.len() {
            return Err(StateVectorError::WrongLength {
                expected: self.amplitudes.len(),
                actual: basis_state + 1,
            });
        }
        Ok((basis_state.count_ones() % 2) as u8)
    }

    /// Probability that a basis state has the requested bit value on one qubit.
    pub fn conditional_bit_probability(
        &self,
        qubit: usize,
        bit: u8,
        fixed_qubit: usize,
        fixed_bit: u8,
    ) -> Result<f64, StateVectorError> {
        self.check_qubit(qubit)?;
        self.check_qubit(fixed_qubit)?;
        if qubit == fixed_qubit {
            return Err(StateVectorError::SameQubit);
        }
        if bit > 1 || fixed_bit > 1 {
            return Err(StateVectorError::WrongLength {
                expected: 2,
                actual: bit.max(fixed_bit) as usize + 1,
            });
        }
        let target = self.amplitudes.iter().enumerate().filter(|(index, _)| {
            let bit_value = ((index >> qubit) & 1) as u8;
            let fixed_value = ((index >> fixed_qubit) & 1) as u8;
            bit_value == bit && fixed_value == fixed_bit
        });
        Ok(target.map(|(_, amplitude)| amplitude.norm_sqr()).sum())
    }

    /// Bitwise AND of all set-bit positions in a basis state.
    pub fn basis_bitwise_and(&self, basis_state: usize) -> Result<usize, StateVectorError> {
        if basis_state >= self.amplitudes.len() {
            return Err(StateVectorError::WrongLength {
                expected: self.amplitudes.len(),
                actual: basis_state + 1,
            });
        }
        let mut result = usize::MAX;
        for bit in 0..self.qubits {
            if basis_state & (1 << bit) != 0 {
                result &= bit;
            }
        }
        Ok(if result == usize::MAX { 0 } else { result })
    }

    /// Probability mass carried by the computational basis states with even parity.
    pub fn even_parity_probability(&self) -> f64 {
        self.amplitudes
            .iter()
            .enumerate()
            .filter(|(index, _)| index.count_ones() % 2 == 0)
            .map(|(_, amplitude)| amplitude.norm_sqr())
            .sum()
    }

    /// Rényi entropy of order 2 over the basis-state distribution, in bits.
    pub fn renyi2_entropy(&self) -> f64 {
        let purity: f64 = self.amplitudes.iter().map(|a| a.norm_sqr().powi(2)).sum();
        if purity <= f64::EPSILON {
            0.0
        } else {
            -purity.log2()
        }
    }

    /// Cumulative probability distribution indexed by basis state.
    /// Jensen-Shannon divergence between two basis-state probability distributions.
    pub fn jensen_shannon_divergence(&self, other: &StateVector) -> Result<f64, StateVectorError> {
        if self.qubits != other.qubits {
            return Err(StateVectorError::WrongLength {
                expected: self.amplitudes.len(),
                actual: other.amplitudes.len(),
            });
        }
        let mix = 0.5;
        let kl = |p: f64, q: f64| {
            if p <= 0.0 { 0.0 } else { p * (p / q).log2() }
        };
        let mut divergence = 0.0;
        for (a, b) in self.amplitudes.iter().zip(other.amplitudes.iter()) {
            let p = a.norm_sqr();
            let q = b.norm_sqr();
            let m = mix * p + mix * q;
            divergence += 0.5 * kl(p, m) + 0.5 * kl(q, m);
        }
        Ok(divergence)
    }

    /// Probability mass on basis states whose popcount is odd.
    pub fn odd_parity_probability(&self) -> f64 {
        let mut total = 0.0;
        for (index, amplitude) in self.amplitudes.iter().enumerate() {
            if index.count_ones() % 2 == 1 {
                total += amplitude.norm_sqr();
            }
        }
        total
    }

    /// Largest absolute difference between consecutive basis-state probabilities.
    pub fn max_probability_step(&self) -> f64 {
        self.amplitudes
            .windows(2)
            .map(|pair| (pair[0].norm_sqr() - pair[1].norm_sqr()).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Probability that two independently prepared copies collapse to the same basis state.
    pub fn collision_probability(&self) -> f64 {
        self.amplitudes
            .iter()
            .map(|amplitude| amplitude.norm_sqr().powi(2))
            .sum()
    }

    /// Effective number of outcomes, defined as the inverse collision probability.
    pub fn participation_ratio(&self) -> f64 {
        let collision = self.collision_probability();
        if collision <= f64::EPSILON {
            0.0
        } else {
            1.0 / collision
        }
    }

    /// Probability mass outside the most likely basis state.
    pub fn tail_probability(&self) -> f64 {
        (1.0 - self.max_probability()).max(0.0)
    }

    /// Shannon entropy of the ideal measurement distribution, in bits.
    pub fn shannon_entropy(&self) -> f64 {
        self.amplitudes
            .iter()
            .map(|amplitude| amplitude.norm_sqr())
            .filter(|probability| *probability > f64::EPSILON)
            .map(|probability| -probability * probability.log2())
            .sum()
    }

    /// Gini coefficient of the ideal measurement distribution.
    pub fn probability_gini(&self) -> f64 {
        let probabilities = self.probabilities();
        let n = probabilities.len() as f64;
        if n <= 1.0 {
            return 0.0;
        }
        let mut sorted = probabilities;
        sorted.sort_by(|a, b| a.total_cmp(b));
        let weighted: f64 = sorted
            .iter()
            .enumerate()
            .map(|(index, probability)| (index as f64 + 1.0) * probability)
            .sum();
        2.0 * weighted / n - (n + 1.0) / n
    }

    /// Mean absolute difference between consecutive basis-state probabilities.
    pub fn probability_smoothness(&self) -> f64 {
        let probabilities = self.probabilities();
        if probabilities.len() < 2 {
            return 0.0;
        }
        probabilities
            .windows(2)
            .map(|p| (p[0] - p[1]).abs())
            .sum::<f64>()
            / (probabilities.len() - 1) as f64
    }

    /// Effective branch count from participation ratio, capped by the state dimension.
    pub fn effective_outcomes(&self) -> usize {
        let dimension = self.amplitudes.len().max(1);
        (self.participation_ratio().round() as usize).clamp(1, dimension)
    }

    /// Shannon entropy per active qubit, normalized to `[0, 1]`.
    pub fn entropy_per_qubit(&self) -> f64 {
        let qubits = self.qubits.max(1) as f64;
        (self.shannon_entropy() / qubits).min(1.0)
    }

    /// Largest basis-state probability.
    pub fn max_probability(&self) -> f64 {
        self.amplitudes
            .iter()
            .map(|amplitude| amplitude.norm_sqr())
            .fold(0.0_f64, f64::max)
    }

    /// Peak concentration of the measurement distribution: `max(p) * dimension`.
    pub fn peak_concentration(&self) -> f64 {
        let dimension = self.amplitudes.len();
        if dimension == 0 {
            return 0.0;
        }
        self.max_probability() * dimension as f64
    }

    /// Effective dimension from the largest basis-state probability.
    pub fn inverse_max_probability(&self) -> f64 {
        let largest = self.max_probability();
        if largest <= f64::EPSILON {
            0.0
        } else {
            1.0 / largest
        }
    }

    pub fn cumulative_probabilities(&self) -> Vec<f64> {
        let mut cumulative = 0.0;
        self.amplitudes
            .iter()
            .map(|a| {
                cumulative += a.norm_sqr();
                cumulative
            })
            .collect()
    }

    pub fn normalize(&mut self) {
        let norm_sqr = self.norm_sqr();
        if norm_sqr > f64::EPSILON {
            let scale = 1.0 / norm_sqr.sqrt();
            for amplitude in &mut self.amplitudes {
                amplitude.re *= scale;
                amplitude.im *= scale;
            }
        }
    }

    pub fn is_normalized(&self) -> bool {
        (self.norm_sqr() - 1.0).abs() < 1e-9
    }

    pub fn apply_1q(
        &mut self,
        target: usize,
        matrix: [[Complex; 2]; 2],
    ) -> Result<(), StateVectorError> {
        self.check_qubit(target)?;
        let stride = 1 << target;
        let block = stride << 1;
        for base in (0..self.amplitudes.len()).step_by(block) {
            for offset in 0..stride {
                let i0 = base + offset;
                let i1 = i0 + stride;
                let a0 = self.amplitudes[i0];
                let a1 = self.amplitudes[i1];
                self.amplitudes[i0] = matrix[0][0] * a0 + matrix[0][1] * a1;
                self.amplitudes[i1] = matrix[1][0] * a0 + matrix[1][1] * a1;
            }
        }
        Ok(())
    }

    pub fn apply_2q(
        &mut self,
        control: usize,
        target: usize,
        matrix: [[Complex; 2]; 2],
    ) -> Result<(), StateVectorError> {
        self.check_qubit(control)?;
        self.check_qubit(target)?;
        if control == target {
            return Err(StateVectorError::SameQubit);
        }
        let target_bit = 1 << target;
        let control_bit = 1 << control;
        for index in 0..self.amplitudes.len() {
            if index & control_bit == 0 || index & target_bit != 0 {
                continue;
            }
            let i0 = index;
            let i1 = index | target_bit;
            let a0 = self.amplitudes[i0];
            let a1 = self.amplitudes[i1];
            self.amplitudes[i0] = matrix[0][0] * a0 + matrix[0][1] * a1;
            self.amplitudes[i1] = matrix[1][0] * a0 + matrix[1][1] * a1;
        }
        Ok(())
    }

    pub fn apply_cp(
        &mut self,
        control: usize,
        target: usize,
        theta: f64,
    ) -> Result<(), StateVectorError> {
        self.check_qubit(control)?;
        self.check_qubit(target)?;
        if control == target {
            return Err(StateVectorError::SameQubit);
        }
        let target_bit = 1 << target;
        let control_bit = 1 << control;
        let phase = Complex::from_polar(1.0, theta);
        for index in 0..self.amplitudes.len() {
            if index & control_bit != 0 && index & target_bit != 0 {
                self.amplitudes[index] *= phase;
            }
        }
        Ok(())
    }

    pub fn reduced_density_matrix(
        &self,
        qubit: usize,
    ) -> Result<[[Complex; 2]; 2], StateVectorError> {
        self.check_qubit(qubit)?;
        let bit = 1 << qubit;
        let mut rho = [[Complex::ZERO; 2]; 2];
        for index in 0..self.amplitudes.len() {
            if index & bit != 0 {
                continue;
            }
            let a0 = self.amplitudes[index];
            let a1 = self.amplitudes[index | bit];
            rho[0][0] += a0 * a0.conj();
            rho[0][1] += a0 * a1.conj();
            rho[1][0] += a1 * a0.conj();
            rho[1][1] += a1 * a1.conj();
        }
        Ok(rho)
    }

    /// Expectation value of the Pauli Z operator on one qubit.
    pub fn expectation_z(&self, qubit: usize) -> Result<f64, StateVectorError> {
        let rho = self.reduced_density_matrix(qubit)?;
        Ok(rho[0][0].re - rho[1][1].re)
    }

    pub fn bloch_vector(&self, qubit: usize) -> Result<[f64; 3], StateVectorError> {
        let rho = self.reduced_density_matrix(qubit)?;
        Ok([
            2.0 * rho[0][1].re,
            -2.0 * rho[0][1].im,
            rho[0][0].re - rho[1][1].re,
        ])
    }

    pub fn measure(&self, rng: &mut StdRng) -> usize {
        let threshold: f64 = rng.random();
        let mut cumulative = 0.0;
        for (index, amplitude) in self.amplitudes.iter().enumerate() {
            cumulative += amplitude.norm_sqr();
            if threshold < cumulative {
                return index;
            }
        }
        self.amplitudes.len() - 1
    }

    pub fn sample(&self, shots: usize, seed: u64) -> Vec<usize> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..shots).map(|_| self.measure(&mut rng)).collect()
    }

    pub fn check_qubit(&self, index: usize) -> Result<(), StateVectorError> {
        if index >= self.qubits {
            Err(StateVectorError::QubitOutOfRange {
                index,
                qubits: self.qubits,
            })
        } else {
            Ok(())
        }
    }
}
