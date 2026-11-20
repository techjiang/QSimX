pub mod algorithm;
pub mod circuit;
pub mod gate;
pub mod measure;
pub mod state;

pub mod prelude {
    pub use crate::algorithm::{
        grover_search, inverse_qft, qft, GroverOracle, GroverResult, QftResult,
    };
    pub use crate::circuit::{
        Circuit, QasmError, SimulationResult, SimulatorError, StateVectorSimulator,
    };
    pub use crate::gate::Gate;
    pub use crate::measure::{computational_probability, expectation};
    pub use crate::state::{Basis, Complex, StateVector, StateVectorError};
    pub use crate::VERSION;
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
