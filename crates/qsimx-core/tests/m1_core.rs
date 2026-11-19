use qsimx::prelude::*;

#[test]
fn bell_state_has_only_00_and_11() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.probability(0b00) - 0.5).abs() < 1e-12);
    assert!(state.probability(0b01) < 1e-12);
    assert!(state.probability(0b10) < 1e-12);
    assert!((state.probability(0b11) - 0.5).abs() < 1e-12);
    assert!(state.is_normalized());
}

#[test]
fn ghz_state_has_two_equal_outcomes() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let result = StateVectorSimulator::new().run(&circuit, 3000).unwrap();
    let counts = result.counts();
    assert_eq!(counts.len(), 2);
    assert_eq!(counts[&0b000] + counts[&0b111], 3000);
    assert!((counts[&0b000] as f64 / 3000.0 - 0.5).abs() < 0.05);
}

#[test]
fn phase_kickback_ghz_is_correct() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).z(0).cnot(0, 1).cnot(0, 2);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.probability(0b000) - 0.5).abs() < 1e-12);
    assert!((state.probability(0b111) - 0.5).abs() < 1e-12);
}

#[test]
fn bloch_vector_tracks_pauli_states() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    let zero = sim.run_state(&circuit).unwrap();
    assert_eq!(zero.bloch_vector(0).unwrap(), [0.0, 0.0, 1.0]);

    circuit.x(0);
    let one = sim.run_state(&circuit).unwrap();
    assert_eq!(one.bloch_vector(0).unwrap(), [0.0, 0.0, -1.0]);

    let mut circuit = Circuit::new(1);
    circuit.h(0);
    let plus = sim.run_state(&circuit).unwrap();
    let [x, y, z] = plus.bloch_vector(0).unwrap();
    assert!((x - 1.0).abs() < 1e-12);
    assert!(y.abs() < 1e-12);
    assert!(z.abs() < 1e-12);
}

#[test]
fn toffoli_and_swap_work() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(3);
    circuit.x(0).x(1).toffoli(0, 1, 2);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.probability(0b111) - 1.0).abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.x(0).swap(0, 1);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.probability(0b10) - 1.0).abs() < 1e-12);
}

#[test]
fn qasm_export_contains_full_circuit() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let qasm = circuit.to_qasm().unwrap();
    assert!(qasm.contains("OPENQASM 2.0;"));
    assert!(qasm.contains("qreg q[2];"));
    assert!(qasm.contains("h q[0];"));
    assert!(qasm.contains("cx q[0],q[1];"));
    assert!(qasm.contains("measure q[1] -> c[1];"));
}

#[test]
fn rotation_circuit_is_unitary() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit
        .rx(0, 0.7)
        .ry(0, -0.3)
        .rz(0, 1.2)
        .rz(0, -1.2)
        .ry(0, 0.3)
        .rx(0, -0.7);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.probability(0) - 1.0).abs() < 1e-10);
}

#[test]
fn measurement_basis_expectations_are_right() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.h(0);
    let state = sim.run_state(&circuit).unwrap();
    assert!((expectation(&state, Basis::X, 0).unwrap() - 1.0).abs() < 1e-12);
    assert!(expectation(&state, Basis::Y, 0).unwrap().abs() < 1e-12);
    assert!(expectation(&state, Basis::Z, 0).unwrap().abs() < 1e-12);
}


#[test]
fn fidelity_detects_identical_and_orthogonal_states() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    let zero = sim.run_state(&circuit).unwrap();
    assert!((zero.fidelity(&zero).unwrap() - 1.0).abs() < 1e-12);

    circuit.x(0);
    let one = sim.run_state(&circuit).unwrap();
    assert!(zero.fidelity(&one).unwrap() < 1e-12);
}
