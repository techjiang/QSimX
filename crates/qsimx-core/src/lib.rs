pub mod algorithm;
pub mod circuit;
pub mod gate;
pub mod measure;
pub mod state;

pub mod prelude {
    pub use crate::VERSION;
    pub use crate::algorithm::{
        GroverAnalysisResult, GroverError, GroverOracle, GroverResult, QftResult,
        QftRoundTripError, grover_amplification_gain, grover_success_roughness, grover_ascending_steps,
        grover_first_success_iteration, grover_iterations, grover_peak_delta, grover_peak_drop,
        grover_peak_iteration, grover_peak_step, grover_peak_width, grover_peak_width_ratio,
        grover_recommended_gain, grover_search, grover_success_area, grover_success_cumulative,
        grover_success_curve, grover_success_kurtosis, grover_success_path_length,
        grover_success_skewness, grover_success_variance, inverse_qft, qft, qft_adjacent_phase_gap,
        qft_adjacent_phase_step, qft_average_probability, qft_density_balance, qft_density_gap,
        qft_gate_balance, qft_gate_density, qft_gate_fidelity, qft_gate_symmetry,
        qft_inverse_distribution_gap, qft_pair_is_reversible, qft_phase_balance,
        qft_phase_circular_skewness, qft_phase_circular_variance, qft_phase_resultant_length,
        qft_probability_imbalance, qft_roundtrip_fidelity, qft_sparsity, qft_uniform_deviation,
        recommended_grover_iterations,
    };
    pub use crate::circuit::{
        Circuit, QasmError, SimulationResult, SimulatorError, StateVectorSimulator, apply_gate,
    };
    pub use crate::gate::Gate;
    pub use crate::measure::{computational_probability, expectation};
    pub use crate::state::{Basis, Complex, StateVector, StateVectorError};
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
