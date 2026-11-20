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


#[test]
fn iteration7_helpers_work() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let state = sim.run_state(&circuit).unwrap();

    let cumulative = state.cumulative_probabilities();
    assert_eq!(cumulative.len(), 8);
    assert!((cumulative[0] - 0.5).abs() < 1e-12);
    assert!((cumulative[7] - 1.0).abs() < 1e-12);

    assert_eq!(circuit.depth(), 3);
    assert_eq!(circuit.gate_counts().get("H"), Some(&1));
    assert_eq!(circuit.gate_counts().get("CX"), Some(&2));

    let result = sim.run(&circuit, 1024).unwrap();
    assert!((result.shannon_entropy() - 1.0).abs() < 0.05);
    assert_eq!(result.top_counts(2).len(), 2);
    assert!(state.amplitude_string(0, 6).contains('+'));
}

#[test]
fn iteration8_state_and_circuit_helpers_work() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.h(0);
    let state = sim.run_state(&circuit).unwrap();

    assert!((state.zero_probability(0).unwrap() - 0.5).abs() < 1e-12);
    assert!(state.phase(0).abs() < 1e-12);
    assert_eq!(circuit.summary(), (1, 1, 1));
    assert_eq!(circuit.gates()[0].qasm_name(), "h");
    assert!((state.fidelity(&state).unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration8_result_fidelity_matches_state() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let result = sim.run(&circuit, 128).unwrap();
    let reference = sim.run_state(&circuit).unwrap();
    assert!((result.state_fidelity(&reference).unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration9_one_probability_complements_zero_probability() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.h(0);
    let state = sim.run_state(&circuit).unwrap();
    let zero = state.zero_probability(0).unwrap();
    let one = state.one_probability(0).unwrap();
    assert!((zero - 0.5).abs() < 1e-12);
    assert!((one - 0.5).abs() < 1e-12);
    assert!((zero + one - 1.0).abs() < 1e-12);
}

#[test]
fn iteration9_bell_pair_builder_creates_expected_state() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.bell_pair(0, 1);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.probability(0) - 0.5).abs() < 1e-12);
    assert!((state.probability(3) - 0.5).abs() < 1e-12);
}

#[test]
fn iteration9_most_likely_reports_dominant_outcome() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.x(0);
    let result = sim.run(&circuit, 256).unwrap();
    assert_eq!(result.most_likely(), Some((1, 1.0)));
}

#[test]
fn iteration10_result_probability_and_min_entropy_work() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.h(0);
    let result = sim.run(&circuit, 2048).unwrap();
    let p0 = result.probability_by_state(0);
    let p1 = result.probability_by_state(1);
    assert!((p0 + p1 - 1.0).abs() < 1e-12);
    assert!((p0 - 0.5).abs() < 0.05);
    assert!((result.min_entropy() - 1.0).abs() < 0.05);

    let empty = StateVectorSimulator::new().run(&circuit, 0).unwrap();
    assert_eq!(empty.probability_by_state(0), 0.0);
    assert_eq!(empty.min_entropy(), 0.0);
}

#[test]
fn iteration10_expectation_z_matches_bloch_z() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.ry(0, 0.7);
    let state = sim.run_state(&circuit).unwrap();
    let z = state.expectation_z(0).unwrap();
    assert!((z - state.bloch_vector(0).unwrap()[2]).abs() < 1e-12);
    assert!((z - 0.7_f64.cos()).abs() < 1e-12);
}

#[test]
fn iteration10_sdg_inverts_s() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(1);
    circuit.s(0).sdg(0);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.probability(0) - 1.0).abs() < 1e-12);

    let mut circuit = Circuit::new(1);
    circuit.sdg(0);
    let qasm = circuit.to_qasm().unwrap();
    assert!(qasm.contains("sdg q[0];"));
    assert_eq!(circuit.gates()[0].name(), "SDG");
}
