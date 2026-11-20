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
fn qft_maps_computational_basis_to_product_phases() {
    let sim = StateVectorSimulator::new();
    let targets = [0, 1];
    let q = qft(&targets, 2).unwrap();
    assert_eq!(q.circuit.gates().len(), 4);

    let mut input = Circuit::new(2);
    input.x(0);
    let basis = sim.run_state(&input).unwrap();
    let mut expected = basis.clone();
    for gate in q.circuit.gates() {
        apply_gate(gate, &mut expected).unwrap();
    }

    let expected_amplitudes = expected.amplitudes();
    assert!((expected_amplitudes[0].abs() - 0.5).abs() < 1e-10);
    assert!((expected_amplitudes[1].abs() - 0.5).abs() < 1e-10);
    assert!((expected_amplitudes[2].abs() - 0.5).abs() < 1e-10);
    assert!((expected_amplitudes[3].abs() - 0.5).abs() < 1e-10);
}

#[test]
fn inverse_qft_undoes_qft() {
    let sim = StateVectorSimulator::new();
    let targets = [0, 1, 2];
    let mut input = Circuit::new(3);
    input.x(0).h(1).s(2);
    let original = sim.run_state(&input).unwrap();

    let forward = qft(&targets, 3).unwrap();
    let mut transformed = original.clone();
    for gate in forward.circuit.gates() {
        apply_gate(gate, &mut transformed).unwrap();
    }

    let inverse = inverse_qft(&targets, 3).unwrap();
    for gate in inverse.circuit.gates() {
        apply_gate(gate, &mut transformed).unwrap();
    }
    assert!((transformed.fidelity(&original).unwrap() - 1.0).abs() < 1e-10);
}

#[test]
fn grover_amplifies_the_marked_state() {
    for qubits in 1..=3 {
        let marked = (1usize << qubits) - 1;
        let iterations = recommended_grover_iterations(qubits);
        let result = grover_search(qubits, &GroverOracle::marked(marked), iterations).unwrap();
        assert!(
            result.marked_probability >= (1.0 / (1usize << qubits) as f64) - 1e-12,
            "qubits={qubits}"
        );
        assert_eq!(result.iterations, iterations);
    }
}

#[test]
fn grover_rejects_invalid_input() {
    assert!(grover_search(0, &GroverOracle::marked(0), 1).is_err());
    assert!(grover_search(2, &GroverOracle::marked(4), 1).is_err());
    assert!(qft(&[2], 2).is_err());
    assert!(qft(&[0, 0], 2).is_err());
}

#[test]
fn grover_iterations_avoid_small_space_overshoot() {
    assert_eq!(recommended_grover_iterations(1), 1);
    assert_eq!(recommended_grover_iterations(2), 1);
    assert_eq!(recommended_grover_iterations(3), 2);
    assert_eq!(recommended_grover_iterations(4), 3);
    assert_eq!(grover_iterations(10, 0), 1);
    assert_eq!(grover_iterations(0, 1), 1);
}

#[test]
fn grover_search_with_custom_oracle_reports_peak_probability() {
    let mut oracle = Circuit::new(2);
    oracle.x(0).x(1).z(0).z(1).x(0).x(1);
    let result = grover_search(2, &GroverOracle::circuit(oracle), 1).unwrap();
    assert!((result.marked_probability - 0.25).abs() < 1e-12);
    assert_eq!(result.iterations, 1);
}

#[test]
fn iteration10_analysis_helpers_work() {
    let sim = StateVectorSimulator::new();

    let mut w = Circuit::new(3);
    w.w_state(&[0, 1, 2]).unwrap();
    let state = sim.run_state(&w).unwrap();
    assert_eq!(state.significant_basis_count(1e-12), 2);
    assert!((state.probabilities().iter().sum::<f64>() - 1.0).abs() < 1e-10);

    let mut bell = Circuit::new(2);
    bell.bell_pair(0, 1);
    let result = sim.run(&bell, 1024).unwrap();
    assert!((result.effective_outcomes() - 2.0).abs() < 0.1);
}

#[test]
fn iteration10_grover_iteration_helper_is_bounded() {
    assert_eq!(recommended_grover_iterations(0), 1);
    assert_eq!(recommended_grover_iterations(1), 1);
    assert_eq!(recommended_grover_iterations(2), 1);
    assert!(recommended_grover_iterations(8) >= 1);
    assert_eq!(grover_iterations(8, 0), 1);
}

#[test]
fn iteration11_entanglement_and_purity_are_complementary() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.bell_pair(0, 1);
    let state = sim.run_state(&circuit).unwrap();
    assert!((state.entanglement_entropy(0).unwrap() - 1.0).abs() < 1e-12);
    assert!((state.purity(0).unwrap() - 0.5).abs() < 1e-12);

    let mut separable = Circuit::new(2);
    separable.h(0);
    let product = sim.run_state(&separable).unwrap();
    assert!(product.entanglement_entropy(0).unwrap().abs() < 1e-12);
    assert!((product.purity(0).unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration11_grover_success_curve_starts_flat() {
    let oracle = GroverOracle::marked(0b11);
    let curve = grover_success_curve(2, &oracle, 2).unwrap();
    assert_eq!(curve.len(), 3);
    assert!((curve[0] - 0.25).abs() < 1e-12);
    assert!(curve[1] >= curve[0]);
}

#[test]
fn iteration11_qft_roundtrip_preserves_state() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).s(2);
    let state = sim.run_state(&circuit).unwrap();
    let fidelity = qft_roundtrip_fidelity(&[0, 1, 2], 3, &state).unwrap();
    assert!((fidelity - 1.0).abs() < 1e-10);
}

#[test]
fn iteration11_entropy_stability_is_nonnegative() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let result = sim.run(&circuit, 1024).unwrap();
    assert!(result.entropy_stability() >= 0.0);
}

#[test]
fn iteration12_total_variation_distance_bounds() {
    let zero = StateVector::zero(2).unwrap();
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let bell = StateVectorSimulator::new().run_state(&circuit).unwrap();
    let distance = zero.total_variation_distance(&bell).unwrap();
    assert!((distance - 0.5).abs() < 1e-12);
    assert_eq!(zero.total_variation_distance(&zero).unwrap(), 0.0);
}

#[test]
fn iteration12_state_purity_reports_max_probability() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((Circuit::state_purity(&circuit, &state) - 0.5).abs() < 1e-12);
}

#[test]
fn iteration12_qft_sparsity_counts_modes() {
    let result = qft(&[0, 1, 2], 3).unwrap();
    assert_eq!(qft_sparsity(&result, 1e-12), 8);
}

#[test]
fn iteration13_probability_span_is_bounded() {
    let zero = StateVector::zero(2).unwrap();
    assert!((zero.probability_span() - 1.0).abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let bell = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((bell.probability_span() - 0.5).abs() < 1e-12);
}

#[test]
fn iteration13_histogram_text_is_scaled() {
    let mut circuit = Circuit::new(1);
    circuit.x(0);
    let result = StateVectorSimulator::new().run(&circuit, 8).unwrap();
    let text = result.histogram_text(1, 4);
    assert!(text.contains("1 | #### 8"), "histogram was: {text}");
    let empty_bar = result.histogram_text(1, 0);
    assert!(
        empty_bar.contains("1 |  8"),
        "empty histogram was: {empty_bar}"
    );
}

#[test]
fn iteration13_layered_circuit_tracks_active_qubits() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).add_layered(Gate::X(2));
    assert_eq!(circuit.distinct_touching_qubits(), vec![0, 2]);
    assert_eq!(circuit.depth(), 2);
}

#[test]
fn iteration13_qft_density_and_grover_gain() {
    let qft = qft(&[0, 1, 2], 3).unwrap();
    assert!((qft_gate_density(&qft) - 7.0 / 3.0).abs() < 1e-12);
    assert!((grover_amplification_gain(2, 0.25) - 1.0).abs() < 1e-12);
    assert!((grover_amplification_gain(3, 0.5) - 4.0).abs() < 1e-12);
}

#[test]
fn iteration14_basis_parity_and_parity_counts_are_consistent() {
    let mut circuit = Circuit::new(2);
    circuit.x(1);
    let result = StateVectorSimulator::new().run(&circuit, 8).unwrap();
    assert_eq!(result.parity_counts(), [0, 8]);
    assert_eq!(result.state().basis_parity(0b01).unwrap(), 1);
    assert_eq!(result.state().basis_parity(0b10).unwrap(), 1);
    assert!(result.state().basis_parity(4).is_err());
}

#[test]
fn iteration14_grover_cumulative_is_monotonic() {
    let oracle = GroverOracle::marked(0b111);
    let cumulative = grover_success_cumulative(3, &oracle, 3).unwrap();
    assert_eq!(cumulative.len(), 4);
    assert!(cumulative.windows(2).all(|pair| pair[1] >= pair[0]));
}

#[test]
fn iteration14_qft_average_probability_is_uniform() {
    let result = qft(&[0, 1, 2], 3).unwrap();
    assert!((qft_average_probability(&result) - 1.0 / 8.0).abs() < 1e-12);
}

#[test]
fn iteration14_w_state_helper_builds_expected_ghz_state() {
    let mut circuit = Circuit::new(3);
    let targets = [0usize, 1, 2];
    circuit.w_state(&targets).unwrap();
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.probability(0b000) - 0.5).abs() < 1e-12);
    assert!((state.probability(0b111) - 0.5).abs() < 1e-12);
}

#[test]
fn iteration14_w_state_rejects_empty_and_duplicate_targets() {
    let mut circuit = Circuit::new(2);
    assert!(circuit.w_state(&[]).is_err());
    assert!(circuit.w_state(&[2]).is_err());
}
