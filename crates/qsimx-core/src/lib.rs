pub mod algorithm;
pub mod circuit;
pub mod gate;
pub mod measure;
pub mod state;

pub mod prelude {
    pub use crate::VERSION;
    pub use crate::algorithm::{
        GroverOracle, GroverResult, QftResult, QftRoundTripError, grover_amplification_gain,
        grover_iterations, grover_search, grover_success_cumulative, grover_success_curve,
        inverse_qft, qft, qft_average_probability, qft_density_gap, qft_gate_density, qft_gate_symmetry,
        qft_roundtrip_fidelity, qft_sparsity, recommended_grover_iterations,
    };
    pub use crate::circuit::{
        Circuit, QasmError, SimulationResult, SimulatorError, StateVectorSimulator, apply_gate,
    };
    pub use crate::gate::Gate;
    pub use crate::measure::{computational_probability, expectation};
    pub use crate::state::{Basis, Complex, StateVector, StateVectorError};
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
