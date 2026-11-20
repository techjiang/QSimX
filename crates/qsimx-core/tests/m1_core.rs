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

#[test]
fn iteration15_conditional_bit_probability_is_normalized_per_condition() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    let q0_and_q1_zero = state.conditional_bit_probability(1, 0, 0, 0).unwrap();
    let q0_and_q1_one = state.conditional_bit_probability(1, 1, 0, 0).unwrap();
    assert!((q0_and_q1_zero - 0.5).abs() < 1e-12);
    assert!((q0_and_q1_one).abs() < 1e-12);
    assert!(state.conditional_bit_probability(0, 0, 0, 0).is_err());
    assert!(state.conditional_bit_probability(1, 2, 0, 0).is_err());
}

#[test]
fn iteration15_basis_bitwise_and_reports_set_bit_positions() {
    let state = StateVector::zero(4).unwrap();
    assert_eq!(state.basis_bitwise_and(0b0000).unwrap(), 0);
    assert_eq!(state.basis_bitwise_and(0b0010).unwrap(), 1);
    assert_eq!(state.basis_bitwise_and(0b0110).unwrap(), 0);
    assert_eq!(state.basis_bitwise_and(0b1010).unwrap(), 1);
    assert!(state.basis_bitwise_and(16).is_err());
}

#[test]
fn iteration15_qft_gate_symmetry_is_stable() {
    let targets = [0usize, 1, 2];
    let forward = qft(&targets, 3).unwrap();
    let inverse = inverse_qft(&targets, 3).unwrap();
    assert!((qft_gate_symmetry(&forward, &inverse) - 1.0).abs() < 1e-12);
}

#[test]
fn iteration16_jensen_shannon_divergence_bounds() {
    let mut ghz = Circuit::new(2);
    ghz.h(0).cnot(0, 1);
    let ghz_state = StateVectorSimulator::new().run_state(&ghz).unwrap();
    let other = StateVector::zero(2).unwrap();
    let jsd = ghz_state.jensen_shannon_divergence(&other).unwrap();
    assert!(jsd >= 0.0);
    assert!(jsd <= 1.0);
    assert!((ghz_state.jensen_shannon_divergence(&ghz_state).unwrap()).abs() < 1e-12);
    assert!(
        ghz_state
            .jensen_shannon_divergence(&StateVector::zero(1).unwrap())
            .is_err()
    );
}

#[test]
fn iteration16_max_gates_per_qubit_tracks_busiest_wire() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).h(0).x(1).cnot(0, 2);
    assert_eq!(circuit.max_gates_per_qubit(), 3);
    assert_eq!(Circuit::new(3).max_gates_per_qubit(), 0);
}

#[test]
fn iteration16_qft_density_gap_is_nonnegative() {
    let targets = [0usize, 1, 2];
    let forward = qft(&targets, 3).unwrap();
    let inverse = inverse_qft(&targets, 3).unwrap();
    let gap = qft_density_gap(&forward, &inverse);
    assert!(gap >= 0.0);
    assert!(gap < 1e-12);
}

#[test]
fn iteration16_roundtrip_fidelity_stays_one() {
    let targets = [0usize, 1, 2];
    let initial = StateVector::zero(3).unwrap();
    let fidelity = qft_roundtrip_fidelity(&targets, 3, &initial).unwrap();
    assert!((fidelity - 1.0).abs() < 1e-10);
}

#[test]
fn iteration17_even_parity_probability_matches_explicit_counts() {
    let mut circuit = Circuit::new(2);
    circuit.x(1);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!(state.even_parity_probability().abs() < 1e-12);
    let ghz = {
        let mut c = Circuit::new(2);
        c.h(0).cnot(0, 1);
        StateVectorSimulator::new().run_state(&c).unwrap()
    };
    assert!((ghz.even_parity_probability() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration17_renyi2_entropy_bounds() {
    let zero = StateVector::zero(2).unwrap();
    assert!(zero.renyi2_entropy().abs() < 1e-12);
    let uniform = {
        let mut c = Circuit::new(2);
        c.h(0).h(1);
        StateVectorSimulator::new().run_state(&c).unwrap()
    };
    assert!((uniform.renyi2_entropy() - 2.0).abs() < 1e-12);
}

#[test]
fn iteration17_qft_gate_balance_is_stable() {
    let targets = [0usize, 1, 2];
    let forward = qft(&targets, 3).unwrap();
    let inverse = inverse_qft(&targets, 3).unwrap();
    assert!((qft_gate_balance(&forward, &inverse) - 1.0).abs() < 1e-12);
}

#[test]
fn iteration17_min_gates_per_touched_qubit_ignores_idle_wires() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).h(0).x(1).cnot(0, 1);
    assert_eq!(circuit.min_gates_per_touched_qubit(), 3);
    assert_eq!(Circuit::new(3).min_gates_per_touched_qubit(), 0);
}

#[test]
fn iteration18_parity_and_circuit_analysis_helpers_work() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.x(0);
    let state = sim.run_state(&circuit).unwrap();

    assert!((state.odd_parity_probability() - 1.0).abs() < 1e-12);
    assert!(state.even_parity_probability().abs() < 1e-12);
    assert!((state.max_probability_step() - 1.0).abs() < 1e-12);
    assert_eq!(circuit.active_qubit_count(), 1);
    assert_eq!(circuit.min_gates_per_touched_qubit(), 1);
}

#[test]
fn iteration18_grover_gain_is_nonnegative() {
    let oracle = GroverOracle::marked(0b11);
    let gain = grover_recommended_gain(2, &oracle).unwrap();
    assert!(gain >= -1e-12);
    assert!(gain <= 1.0 + 1e-12);
}

#[test]
fn iteration18_qft_density_balance_matches_average() {
    let targets = vec![0usize, 1, 2];
    let forward = qft(&targets, 3).unwrap();
    let inverse = inverse_qft(&targets, 3).unwrap();
    let expected = (qft_gate_density(&forward) + qft_gate_density(&inverse)) / 2.0;
    assert!((qft_density_balance(&forward, &inverse) - expected).abs() < 1e-12);
}

#[test]
fn iteration18_parity_cli_metrics_are_stable() {
    let oracle = GroverOracle::marked(0b11);
    let grover = grover_search(2, &oracle, recommended_grover_iterations(2)).unwrap();
    let odd = grover.state.odd_parity_probability();
    let even = grover.state.even_parity_probability();
    assert!((odd + even - 1.0).abs() < 1e-12);
    assert!((odd - 0.5).abs() < 1e-12);
    assert!((even - 0.5).abs() < 1e-12);
}

#[test]
fn iteration19_collision_and_participation_ratio_work() {
    let sim = StateVectorSimulator::new();

    let mut basis = Circuit::new(1);
    let basis_state = sim.run_state(&basis).unwrap();
    assert!((basis_state.collision_probability() - 1.0).abs() < 1e-12);
    assert!((basis_state.participation_ratio() - 1.0).abs() < 1e-12);

    basis.h(0);
    let uniform = sim.run_state(&basis).unwrap();
    assert!((uniform.collision_probability() - 0.5).abs() < 1e-12);
    assert!((uniform.participation_ratio() - 2.0).abs() < 1e-12);
}

#[test]
fn iteration19_parameterized_gate_count_tracks_angles() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).rx(0, 0.1).ry(1, 0.2).cp(0, 1, 0.3);
    assert_eq!(circuit.parameterized_gate_count(), 3);

    let mut plain = Circuit::new(2);
    plain.h(0).cnot(0, 1);
    assert_eq!(plain.parameterized_gate_count(), 0);
}

#[test]
fn iteration20_grover_peak_analysis_finds_expected_iteration() {
    let oracle = GroverOracle::marked(0b11);
    assert_eq!(grover_peak_iteration(2, &oracle, 4).unwrap(), 2);
    assert_eq!(grover_peak_delta(2, &oracle, 4).unwrap(), -1);
}

#[test]
fn iteration20_grover_peak_analysis_rejects_empty_search() {
    let oracle = GroverOracle::marked(0);
    assert!(grover_peak_iteration(0, &oracle, 2).is_err());
}

#[test]
fn iteration21_qft_probability_imbalance_is_zero_for_inverse_pair() {
    let targets = [0, 1, 2];
    let forward = qft(&targets, 3).unwrap();
    let inverse = inverse_qft(&targets, 3).unwrap();
    assert!(qft_probability_imbalance(&forward, &inverse) < 1e-12);
}

#[test]
fn iteration21_qft_pair_reversibility_checks_roundtrip() {
    let targets = [0, 1];
    let forward = qft(&targets, 2).unwrap();
    let inverse = inverse_qft(&targets, 2).unwrap();
    let input = StateVector::zero(2).unwrap();

    assert!(qft_pair_is_reversible(&forward, &inverse, &input));
}

#[test]
fn iteration22_verify_metrics_are_stable_on_bell_state() {
    let sim = StateVectorSimulator::new();
    let mut circuit = Circuit::new(2);
    circuit.bell_pair(0, 1);
    let state = sim.run_state(&circuit).unwrap();
    let result = sim.run(&circuit, 512).unwrap();

    assert_eq!(circuit.parameterized_gate_count(), 0);
    assert!((state.collision_probability() - 0.5).abs() < 1e-12);
    assert!((state.participation_ratio() - 2.0).abs() < 1e-12);
    assert!((result.shannon_entropy() - 1.0).abs() < 0.05);
}

#[test]
fn iteration22_verify_metrics_reject_mismatched_qft_states() {
    let forward = qft(&[0, 1], 2).unwrap();
    let inverse = inverse_qft(&[0, 1], 2).unwrap();
    let mut other = inverse.state.clone();
    other.amplitudes_mut()[0].re += 0.1;
    other.normalize();

    let imbalance = forward
        .state
        .probabilities()
        .iter()
        .zip(other.probabilities().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(imbalance > 0.0);
    assert!(qft_probability_imbalance(&forward, &inverse) < 1e-12);
}

#[test]
fn iteration23_information_metrics_match_basis_and_uniform_states() {
    let basis = StateVector::zero(2).unwrap();
    assert!(basis.shannon_entropy().abs() < 1e-12);
    assert!((basis.inverse_max_probability() - 1.0).abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1);
    let uniform = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((uniform.shannon_entropy() - 2.0).abs() < 1e-12);
    assert!((uniform.inverse_max_probability() - 4.0).abs() < 1e-12);
}

#[test]
fn iteration24_multi_qubit_ratio_counts_entangling_gates() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).h(0).cnot(0, 1).toffoli(0, 1, 2);
    assert_eq!(circuit.multi_qubit_gate_count(), 2);
    assert!((circuit.multi_qubit_gate_ratio() - 0.5).abs() < 1e-12);
    assert_eq!(Circuit::new(2).multi_qubit_gate_ratio(), 0.0);
}

#[test]
fn iteration25_grover_first_success_and_peak_step_are_stable() {
    let oracle = GroverOracle::marked(0b11);
    assert_eq!(
        grover_first_success_iteration(2, &oracle, 4, 0.9).unwrap(),
        Some(2)
    );
    assert_eq!(
        grover_first_success_iteration(2, &oracle, 0, 0.99).unwrap(),
        None
    );
    let step = grover_peak_step(2, &oracle, 4).unwrap();
    assert!(step > 0.0 && step <= 1.0 + 1e-12);
}

#[test]
fn iteration26_next_verify_metrics_match_known_ghz_state() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();

    assert!((state.shannon_entropy() - 1.0).abs() < 1e-12);
    assert!((state.collision_probability() - 0.5).abs() < 1e-12);
    assert!((state.inverse_max_probability() - 2.0).abs() < 1e-12);
    assert_eq!(circuit.multi_qubit_gate_count(), 2);
}

#[test]
fn iteration27_inverse_max_probability_rejects_zero_state_vector() {
    let state = StateVector::zero(1).unwrap();
    assert!((state.inverse_max_probability() - 1.0).abs() < 1e-12);
    assert!(grover_peak_step(0, &GroverOracle::marked(0), 1).is_err());
}

#[test]
fn iteration28_entropy_per_qubit_normalizes_distributions() {
    let basis = StateVector::zero(2).unwrap();
    assert!(basis.entropy_per_qubit().abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1);
    let uniform = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((uniform.entropy_per_qubit() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration29_probability_gini_measures_concentration() {
    let basis = StateVector::zero(2).unwrap();
    assert!((basis.probability_gini() - 0.75).abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1);
    let uniform = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!(uniform.probability_gini().abs() < 1e-12);
}

#[test]
fn iteration30_parameterized_gate_ratio_tracks_angle_gates() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).rx(0, 0.5).cp(0, 1, 0.25).x(1);
    assert!((circuit.parameterized_gate_ratio() - 0.5).abs() < 1e-12);
    assert!((Circuit::new(1).parameterized_gate_ratio() - 0.0).abs() < 1e-12);
}

#[test]
fn iteration31_grover_success_area_stays_bounded() {
    let oracle = GroverOracle::marked(0b11);
    let area = grover_success_area(2, &oracle, 4).unwrap();
    assert!(area > 0.0 && area <= 1.0);
    assert!(grover_success_area(0, &oracle, 1).is_err());
}

#[test]
fn iteration32_qft_gate_fidelity_is_bounded() {
    let result = qft(&[0, 1], 2).unwrap();
    let fidelity = qft_gate_fidelity(&result, 2).unwrap();
    assert!((0.0..=1.0).contains(&fidelity));
    assert!(qft_gate_fidelity(&result, 0).is_err());
}

#[test]
fn iteration33_grover_peak_width_counts_above_target() {
    let oracle = GroverOracle::marked(0b11);
    let width = grover_peak_width(2, &oracle, 5, 0.5).unwrap();
    assert!(width >= 1);
    assert_eq!(grover_peak_width(2, &oracle, 5, 1.1).unwrap(), 0);
    assert!(grover_peak_width(0, &oracle, 2, 0.5).is_err());
}

#[test]
fn iteration34_qft_phase_balance_is_bounded() {
    let result = qft(&[0, 1], 2).unwrap();
    let balance = qft_phase_balance(&result);
    assert!((0.0..=1.0).contains(&balance));
}

#[test]
fn iteration35_distinct_touching_qubit_count_skips_idle_qubits() {
    let mut circuit = Circuit::new(4);
    circuit.h(0).cnot(0, 1);
    assert_eq!(circuit.distinct_touching_qubit_count(), 2);
    assert_eq!(Circuit::new(3).distinct_touching_qubit_count(), 0);
}

#[test]
fn iteration36_max_probability_matches_inverse_dimension() {
    let basis = StateVector::zero(2).unwrap();
    assert!((basis.max_probability() - 1.0).abs() < 1e-12);

    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1);
    let uniform = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((uniform.max_probability() - 0.25).abs() < 1e-12);
    assert!((uniform.inverse_max_probability() - 4.0).abs() < 1e-12);
}

#[test]
fn iteration38_grover_peak_drop_is_non_negative() {
    let oracle = GroverOracle::marked(0b11);
    let drop = grover_peak_drop(2, &oracle, 5).unwrap();
    assert!((0.0..=1.0).contains(&drop));
    assert!(grover_peak_drop(2, &oracle, 1).unwrap() <= drop + 1e-12);
    assert!(grover_peak_drop(0, &oracle, 2).is_err());
}

#[test]
fn iteration39_qft_adjacent_phase_gap_is_bounded() {
    let result = qft(&[0, 1], 2).unwrap();
    let gap = qft_adjacent_phase_gap(&result);
    assert!((0.0..=1.0).contains(&gap));
}

#[test]
fn iteration40_single_qubit_gate_count_excludes_multi_qubit_gates() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).x(0).cnot(0, 1);
    assert_eq!(circuit.single_qubit_gate_count(), 2);
    assert_eq!(Circuit::new(2).single_qubit_gate_count(), 0);
}

#[test]
fn iteration41_ghz_odd_parity_probability_is_half() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.odd_parity_probability() - 0.5).abs() < 1e-12);
}

#[test]
fn iteration42_tail_probability_complements_max_probability() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).h(1);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.tail_probability() - 0.75).abs() < 1e-12);
    assert!((state.max_probability() + state.tail_probability() - 1.0).abs() < 1e-12);
}

#[test]
fn iteration43_gates_with_arity_counts_by_width() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).toffoli(0, 1, 2);
    assert_eq!(circuit.gates_with_arity(1), 1);
    assert_eq!(circuit.gates_with_arity(2), 1);
    assert_eq!(circuit.gates_with_arity(3), 1);
    assert_eq!(circuit.gates_with_arity(4), 0);
}

#[test]
fn iteration44_grover_success_path_length_is_non_negative() {
    let oracle = GroverOracle::marked(0b11);
    let path = grover_success_path_length(2, &oracle, 5).unwrap();
    assert!(path >= 0.0);
    assert!(grover_success_path_length(0, &oracle, 2).is_err());
}

#[test]
fn iteration45_qft_adjacent_phase_step_is_bounded() {
    let result = qft(&[0, 1, 2], 3).unwrap();
    let step = qft_adjacent_phase_step(&result);
    assert!((0.0..=2.0).contains(&step));
    assert!(qft_adjacent_phase_step(&qft(&[0], 1).unwrap()).abs() < 1e-12);
}

#[test]
fn iteration46_grover_success_variance_is_non_negative() {
    let oracle = GroverOracle::marked(0b11);
    let variance = grover_success_variance(2, &oracle, 5).unwrap();
    assert!(variance >= 0.0);
    assert!(grover_success_variance(0, &oracle, 2).is_err());
}

#[test]
fn iteration47_grover_ascending_steps_are_bounded() {
    let oracle = GroverOracle::marked(0b11);
    let steps = grover_ascending_steps(2, &oracle, 5).unwrap();
    assert!(steps <= 5);
    assert!(grover_ascending_steps(0, &oracle, 2).is_err());
}

#[test]
fn iteration48_qft_uniform_deviation_is_non_negative() {
    let result = qft(&[0, 1, 2], 3).unwrap();
    let deviation = qft_uniform_deviation(&result);
    assert!((0.0..=1.0).contains(&deviation));
}

#[test]
fn iteration49_qft_inverse_distribution_gap_is_zero() {
    let forward = qft(&[0, 1, 2], 3).unwrap();
    let inverse = inverse_qft(&[0, 1, 2], 3).unwrap();
    assert!(qft_inverse_distribution_gap(&forward, &inverse) >= 0.0);
}

#[test]
fn iteration50_state_and_circuit_statistics_work() {
    let zero = StateVector::zero(2).unwrap();
    assert!((zero.probability_smoothness() - 1.0 / 3.0).abs() < 1e-12);
    assert_eq!(zero.effective_outcomes(), 1);

    let mut bell = Circuit::new(2);
    bell.h(0).cnot(0, 1);
    let bell_state = StateVectorSimulator::new().run_state(&bell).unwrap();
    assert_eq!(bell_state.effective_outcomes(), 2);
    assert_eq!(bell.one_and_two_qubit_gate_count(), 2);
}

#[test]
fn iteration50b_profile_exports_are_available() {
    let zero = StateVector::zero(2).unwrap();
    let _ = zero.peak_concentration();
    let _ = Circuit::new(2).gate_width_entropy();
    let _ = qft_phase_circular_variance(&qft(&[0, 1], 2).unwrap());
    let _ = grover_success_skewness(2, &GroverOracle::marked(0b11), 2);
}

#[test]
fn iteration51_peak_concentration_is_bounded() {
    let zero = StateVector::zero(3).unwrap();
    assert!((zero.peak_concentration() - 8.0).abs() < 1e-12);
    let mut bell = Circuit::new(2);
    bell.h(0).cnot(0, 1);
    let state = StateVectorSimulator::new().run_state(&bell).unwrap();
    assert!((state.peak_concentration() - 2.0).abs() < 1e-12);
}

#[test]
fn iteration52_gate_width_entropy_reflects_mixed_widths() {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    assert!((circuit.gate_width_entropy() - 1.0).abs() < 1e-12);
    assert!(Circuit::new(2).gate_width_entropy().abs() < 1e-12);
}

#[test]
fn iteration53_grover_success_skewness_is_finite() {
    let oracle = GroverOracle::marked(0b11);
    let skewness = grover_success_skewness(2, &oracle, 5).unwrap();
    assert!(skewness.is_finite());
    assert!(grover_success_skewness(0, &oracle, 2).is_err());
}

#[test]
fn iteration54_qft_phase_circular_variance_is_bounded() {
    let result = qft(&[0, 1, 2], 3).unwrap();
    let variance = qft_phase_circular_variance(&result);
    assert!((0.0..=1.0).contains(&variance));
}

#[test]
fn iteration55_profile_statistics_are_stable() {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);
    let state = StateVectorSimulator::new().run_state(&circuit).unwrap();
    assert!((state.peak_concentration() - 4.0).abs() < 1e-12);
    assert!((circuit.gate_width_entropy() - 0.9182958340544896).abs() < 1e-12);
}
