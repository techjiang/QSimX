pub mod algorithm;
pub mod circuit;
pub mod gate;
pub mod measure;
pub mod state;

pub mod prelude {
    pub use crate::VERSION;
    pub use crate::algorithm::{
        GroverOracle, GroverResult, QftResult, grover_search, inverse_qft, qft,
    };
    pub use crate::circuit::{
        Circuit, QasmError, SimulationResult, SimulatorError, StateVectorSimulator, apply_gate,
    };
    pub use crate::gate::Gate;
    pub use crate::measure::{computational_probability, expectation};
    pub use crate::state::{Basis, Complex, StateVector, StateVectorError};
pub use crate::algorithm::{grover_iterations, recommended_grover_iterations};
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
