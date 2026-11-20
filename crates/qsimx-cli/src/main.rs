use std::env;

use qsimx::prelude::{
    Circuit, Gate, GroverOracle, StateVector, StateVectorSimulator, grover_amplification_gain,
    grover_ascending_steps, grover_first_success_iteration, grover_peak_delta, grover_peak_drop,
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

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        return;
    }
    match args[1].as_str() {
        "bell" => run_bell(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1024)),
        "ghz" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_ghz(
                qubits,
                args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024),
            );
        }
        "grover" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_grover(
                qubits,
                args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024),
            );
        }
        "qft" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_qft(
                qubits,
                args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024),
            );
        }
        "standard" => run_standard(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3)),
        "inspect" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_inspect(qubits);
        }
        "summary9" => run_summary9(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2)),
        "grover-m2" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_grover_m2(qubits);
        }
        "qft-m2" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_qft_m2(qubits);
        }
        "verify-m2" => run_verify_m2(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3)),
        "verify-m2-plus" => {
            run_verify_m2_plus(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-final" => {
            run_verify_m2_final(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-compact" => {
            run_verify_m2_compact(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-report" => {
            run_verify_m2_report(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-cond" => {
            run_verify_m2_cond(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-dist" => {
            run_verify_m2_dist(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-balance" => {
            run_verify_m2_balance(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-parity" => {
            run_verify_m2_parity(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-diagnostics" => {
            run_verify_m2_diagnostics(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-next" => {
            run_verify_m2_next(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-advanced" => {
            run_verify_m2_advanced(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-shape" => {
            run_verify_m2_shape(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-curve" => {
            run_verify_m2_curve(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-motion" => {
            run_verify_m2_motion(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-stats" => {
            run_verify_m2_stats(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-profile" => {
            run_verify_m2_profile(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-signature" => {
            run_verify_m2_signature(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "verify-m2-fingerprint" => {
            run_verify_m2_fingerprint(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3))
        }
        "version" => println!("QSimX {}", qsimx::VERSION),
        "qasm" => {
            let mut circuit = Circuit::new(3);
            circuit.h(0).cnot(0, 1).cnot(0, 2);
            match circuit.to_qasm() {
                Ok(qasm) => print!("{qasm}"),
                Err(err) => eprintln!("error: {err}"),
            }
        }
        "help" | "--help" | "-h" => print_usage(),
        other => eprintln!("unknown command: {other}"),
    }
}

fn run_bell(shots: usize) {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);
    run(&circuit, shots);
}

fn run_ghz(qubits: usize, shots: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let mut circuit = Circuit::new(qubits);
    circuit.h(0);
    for target in 1..qubits {
        circuit.cnot(0, target);
    }
    run(&circuit, shots);
}

fn run_grover(qubits: usize, shots: usize) {
    if qubits < 2 {
        eprintln!("error: grover needs at least 2 qubits");
        return;
    }
    let mut circuit = Circuit::new(qubits);
    for q in 0..qubits {
        circuit.h(q);
    }
    // Oracle marks |11...1>.
    if qubits == 2 {
        circuit.cz(0, 1);
    } else {
        circuit.toffoli(0, 1, qubits - 1);
    }
    // Diffuser: H -> X -> multi-controlled Z -> X -> H.
    for q in 0..qubits {
        circuit.h(q).x(q);
    }
    if qubits == 2 {
        circuit.cz(0, 1);
    } else {
        circuit.toffoli(0, 1, qubits - 1);
    }
    for q in 0..qubits {
        circuit.x(q).h(q);
    }
    run(&circuit, shots);
}

fn run_qft(qubits: usize, shots: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    use std::f64::consts::PI;
    let mut circuit = Circuit::new(qubits);
    for target in 0..qubits {
        circuit.h(target);
        for control in target + 1..qubits {
            let angle = PI / (1usize << (control - target)) as f64;
            circuit.cp(control, target, angle);
        }
    }
    for q in 0..qubits / 2 {
        circuit.swap(q, qubits - 1 - q);
    }
    run(&circuit, shots);
}

fn run_standard(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let mut circuit = Circuit::new(qubits);
    circuit.h(0);
    for target in 1..qubits {
        circuit.cnot(0, target);
    }
    match StateVectorSimulator::new().run(&circuit, 1024) {
        Ok(result) => {
            let reference = StateVectorSimulator::new().run_state(&circuit).unwrap();
            println!("standard summary: {:?}", circuit.summary());
            println!(
                "zero probability q0: {:.6}",
                result.state().zero_probability(0).unwrap()
            );
            println!("phase |0>: {:.6}", result.state().phase(0));
            println!(
                "state fidelity: {:.6}",
                result.state_fidelity(&reference).unwrap()
            );
            println!(
                "gate names: {:?}",
                circuit
                    .gates()
                    .iter()
                    .map(|gate| gate.name())
                    .collect::<Vec<_>>()
            );
            println!(
                "qasm names: {:?}",
                circuit
                    .gates()
                    .iter()
                    .map(|gate| gate.qasm_name())
                    .collect::<Vec<_>>()
            );
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn run_inspect(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let mut circuit = Circuit::new(qubits);
    circuit.h(0);
    for target in 1..qubits {
        circuit.cnot(0, target);
    }
    match StateVectorSimulator::new().run(&circuit, 1024) {
        Ok(result) => {
            println!("circuit depth: {}", circuit.depth());
            println!("gate counts: {:?}", circuit.gate_counts());
            println!(
                "nonzero amplitudes: {}",
                result
                    .state()
                    .amplitudes()
                    .iter()
                    .filter(|a| a.norm_sqr() > 1e-12)
                    .count()
            );
            let probs = result.state().probabilities();
            let max = probs.iter().copied().fold(0.0f64, f64::max);
            println!("max probability: {max:.6}");
            println!("shannon entropy: {:.6} bits", result.shannon_entropy());
            println!("top counts: {:?}", result.top_counts(4));
            for (state, _) in result.top_counts(4) {
                println!(
                    "  {state:0width$b}: {}",
                    result.state().amplitude_string(state, 6),
                    width = qubits
                );
            }
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn run_summary9(qubits: usize) {
    if qubits < 2 {
        eprintln!("error: summary9 needs at least 2 qubits");
        return;
    }
    let mut circuit = Circuit::new(qubits);
    circuit.bell_pair(0, 1);
    match StateVectorSimulator::new().run(&circuit, 1024) {
        Ok(result) => {
            println!("summary: {:?}", circuit.summary());
            println!(
                "one probability q0: {:.6}",
                result.state().one_probability(0).unwrap()
            );
            println!("most likely: {:?}", result.most_likely());
            println!(
                "qasm names: {:?}",
                circuit
                    .gates()
                    .iter()
                    .map(|g| g.qasm_name())
                    .collect::<Vec<_>>()
            );
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn run_grover_m2(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let marked = (1usize << qubits) - 1;
    let iterations = recommended_grover_iterations(qubits);
    match grover_search(qubits, &GroverOracle::marked(marked), iterations) {
        Ok(result) => {
            println!("algorithm: grover");
            println!("qubits: {qubits}");
            println!("marked: {marked:0width$b}", width = qubits);
            println!("iterations: {}", result.iterations);
            println!("gates: {}", result.circuit.gates().len());
            println!("marked_probability: {:.6}", result.marked_probability);
            println!(
                "
circuit:\n{}",
                result.circuit.diagram()
            );
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn run_qft_m2(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    match qft(&targets, qubits) {
        Ok(result) => {
            println!("algorithm: qft");
            println!("qubits: {qubits}");
            println!("gates: {}", result.circuit.gates().len());
            println!("\ncircuit:\n{}", result.circuit.diagram());
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn run_verify_m2(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let marked = (1usize << qubits) - 1;
    let iterations = recommended_grover_iterations(qubits);
    let grover = grover_search(qubits, &GroverOracle::marked(marked), iterations);
    let targets: Vec<usize> = (0..qubits).collect();
    let transformed = qft(&targets, qubits);

    let mut w = Circuit::new(qubits);
    w.w_state(&targets).unwrap();
    let w_state = StateVectorSimulator::new().run_state(&w);

    let g = match grover {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let q = match transformed {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let ws = match w_state {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover iterations: {}", g.iterations);
    println!("grover gates: {}", g.circuit.gates().len());
    println!("grover marked probability: {:.6}", g.marked_probability);
    println!("qft gates: {}", q.circuit.gates().len());
    println!(
        "w significant amplitudes: {}",
        ws.significant_basis_count(1e-12)
    );
    println!("qft normalized: {}", q.state.is_normalized());
}

fn run_verify_m2_plus(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let iterations = recommended_grover_iterations(qubits);
    let oracle = GroverOracle::marked(marked);

    let mut circuit = Circuit::new(qubits);
    circuit.h(0);
    for target in 1..qubits {
        circuit.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&circuit) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let curve = match grover_success_curve(qubits, &oracle, iterations) {
        Ok(curve) => curve,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let fidelity = match qft_roundtrip_fidelity(&targets, qubits, &state) {
        Ok(fidelity) => fidelity,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("qft roundtrip fidelity: {fidelity:.6}");
    println!(
        "entanglement entropy q0: {:.6}",
        state.entanglement_entropy(0).unwrap()
    );
    println!("purity q0: {:.6}", state.purity(0).unwrap());
    for (iteration, probability) in curve.iter().enumerate() {
        println!("grover iteration {iteration}: {probability:.6}");
    }
}

fn run_verify_m2_final(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let iterations = recommended_grover_iterations(qubits);
    let oracle = GroverOracle::marked(marked);

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let ghz_state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let zero = match StateVector::zero(qubits) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let qft = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let grover = match grover_search(qubits, &oracle, iterations) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let curve = match grover_success_curve(qubits, &oracle, iterations) {
        Ok(curve) => curve,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let fidelity = match qft_roundtrip_fidelity(&targets, qubits, &zero) {
        Ok(fidelity) => fidelity,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!(
        "tvd(|0>, GHZ): {:.6}",
        zero.total_variation_distance(&ghz_state).unwrap()
    );
    println!(
        "ghz max amplitude probability: {:.6}",
        Circuit::state_purity(&ghz, &ghz_state)
    );
    println!("qft significant modes: {}", qft_sparsity(&qft, 1e-12));
    println!("qft roundtrip fidelity: {fidelity:.6}");
    println!("grover iterations: {}", grover.iterations);
    println!(
        "grover marked probability: {:.6}",
        grover.marked_probability
    );
    println!(
        "grover peak probability: {:.6}",
        curve.iter().copied().fold(0.0_f64, f64::max)
    );
}

fn run_verify_m2_compact(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let iterations = recommended_grover_iterations(qubits);
    let oracle = GroverOracle::marked(marked);
    let grover = match grover_search(qubits, &oracle, iterations) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let qft_result = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut demo = Circuit::new(qubits);
    demo.h(0);
    for target in 1..qubits {
        demo.cnot(0, target);
    }
    let result = match StateVectorSimulator::new().run(&demo, 16) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut layered = Circuit::new(qubits);
    layered.h(0).add_layered(Gate::X(qubits - 1));

    println!("qubits: {qubits}");
    println!("probability span: {:.6}", grover.state.probability_span());
    println!("qft gate density: {:.6}", qft_gate_density(&qft_result));
    println!(
        "grover amplification: {:.6}",
        grover_amplification_gain(qubits, grover.marked_probability)
    );
    println!("layered depth: {}", layered.depth());
    println!("touching qubits: {:?}", layered.distinct_touching_qubits());
    println!("histogram:");
    println!("{}", result.histogram_text(qubits, 32));
}

fn run_verify_m2_report(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let oracle = GroverOracle::marked(marked);
    let grover = match grover_search(qubits, &oracle, recommended_grover_iterations(qubits)) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let qft_result = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let cumulative = match grover_success_cumulative(qubits, &oracle, grover.iterations) {
        Ok(values) => values,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut w_circuit = Circuit::new(qubits);
    if let Err(err) = w_circuit.w_state(&targets) {
        eprintln!("error: {err}");
        return;
    }
    let w_result = match StateVectorSimulator::new().run(&w_circuit, 16) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let state = StateVectorSimulator::new().run_state(&w_circuit).unwrap();

    println!("qubits: {qubits}");
    println!(
        "grover marked probability: {:.6}",
        grover.marked_probability
    );
    println!(
        "grover cumulative peak: {:.6}",
        cumulative.last().copied().unwrap_or(0.0)
    );
    println!(
        "qft average probability: {:.8}",
        qft_average_probability(&qft_result)
    );
    println!("qft gate density: {:.6}", qft_gate_density(&qft_result));
    println!("w-state parity counts: {:?}", w_result.parity_counts());
    println!("basis 0 parity: {:?}", state.basis_parity(0));
}

fn run_verify_m2_cond(qubits: usize) {
    if qubits < 2 {
        eprintln!("error: qubits must be at least 2");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let mut circuit = Circuit::new(qubits);
    circuit.h(0);
    for target in 1..qubits {
        circuit.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&circuit) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let last = qubits - 1;
    let conditional = match state.conditional_bit_probability(last, 0, 0, 0) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    println!("qubits: {qubits}");
    println!("conditional p(q{last}=0 | q0=0): {conditional:.6}");
    println!(
        "basis 0b01 bitwise and: {:?}",
        state.basis_bitwise_and(0b01)
    );
    println!(
        "basis 0b11 bitwise and: {:?}",
        state.basis_bitwise_and(0b11)
    );
    println!(
        "qft gate symmetry: {:.6}",
        qft_gate_symmetry(&forward, &inverse)
    );
    println!(
        "qft average probability: {:.8}",
        qft_average_probability(&forward)
    );
}

fn run_verify_m2_dist(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let oracle = GroverOracle::marked(marked);
    let grover = match grover_search(qubits, &oracle, recommended_grover_iterations(qubits)) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut oracle_circuit = Circuit::new(qubits);
    oracle_circuit.h(0);
    for target in 1..qubits {
        oracle_circuit.cnot(0, target);
    }
    let oracle_state = match StateVectorSimulator::new().run_state(&oracle_circuit) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let jsd = match grover.state.jensen_shannon_divergence(&oracle_state) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let roundtrip = match qft_roundtrip_fidelity(&targets, qubits, &grover.state) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    println!("qubits: {qubits}");
    println!(
        "grover marked probability: {:.6}",
        grover.marked_probability
    );
    println!("grover iterations: {}", grover.iterations);
    println!("js divergence(grover, ghz): {:.6}", jsd);
    println!(
        "total variation distance: {:.6}",
        grover
            .state
            .total_variation_distance(&oracle_state)
            .unwrap_or(1.0)
    );
    println!("qft roundtrip fidelity: {:.6}", roundtrip);
    println!(
        "qft density gap: {:.6}",
        qft_density_gap(&forward, &inverse)
    );
    println!(
        "max gates per qubit: {}",
        grover.circuit.max_gates_per_qubit()
    );
    println!("grover circuit depth: {}", grover.circuit.depth());
}

fn run_verify_m2_balance(qubits: usize) {
    if qubits < 2 {
        eprintln!("error: verify-m2-balance needs at least 2 qubits");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let marked = (1usize << qubits) - 1;
    let oracle = GroverOracle::marked(marked);
    let grover = match grover_search(qubits, &oracle, recommended_grover_iterations(qubits)) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!(
        "grover marked probability: {:.6}",
        grover.marked_probability
    );
    println!(
        "grover even parity probability: {:.6}",
        grover.state.even_parity_probability()
    );
    println!(
        "grover renyi2 entropy: {:.6}",
        grover.state.renyi2_entropy()
    );
    println!(
        "qft gate balance: {:.6}",
        qft_gate_balance(&forward, &inverse)
    );
    println!(
        "max gates per qubit: {}",
        grover.circuit.max_gates_per_qubit()
    );
    println!(
        "min gates per touched qubit: {}",
        grover.circuit.min_gates_per_touched_qubit()
    );
}

fn run_verify_m2_parity(qubits: usize) {
    if qubits < 2 {
        eprintln!("error: verify-m2-parity needs at least 2 qubits");
        return;
    }
    let targets: Vec<usize> = (0..qubits).collect();
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let grover = match grover_search(qubits, &oracle, recommended_grover_iterations(qubits)) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let result = match StateVectorSimulator::new().run(&grover.circuit, 1024) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    println!("qubits: {qubits}");
    println!(
        "odd parity probability: {:.6}",
        grover.state.odd_parity_probability()
    );
    println!(
        "even parity probability: {:.6}",
        grover.state.even_parity_probability()
    );
    println!("measurement parity count: {:?}", result.parity_counts());
    println!(
        "max probability step: {:.6}",
        grover.state.max_probability_step()
    );
    println!("active qubits: {}", grover.circuit.active_qubit_count());
    println!(
        "grover recommended gain: {:.6}",
        match grover_recommended_gain(qubits, &oracle) {
            Ok(gain) => gain,
            Err(err) => {
                eprintln!("error: {err}");
                return;
            }
        }
    );
    println!(
        "qft density balance: {:.6}",
        qft_density_balance(&forward, &inverse)
    );
    println!("circuit depth: {}", grover.circuit.depth());
}

fn run_verify_m2_diagnostics(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: verify-m2-diagnostics needs at least 1 qubit");
        return;
    }
    let marked = (1usize << qubits) - 1;
    let oracle = GroverOracle::marked(marked);
    let targets: Vec<usize> = (0..qubits).collect();

    let peak = match grover_peak_iteration(qubits, &oracle, 4) {
        Ok(peak) => peak,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let peak_delta = match grover_peak_delta(qubits, &oracle, 4) {
        Ok(delta) => delta,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let input = match StateVector::zero(qubits) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut bell = Circuit::new(2);
    bell.bell_pair(0, 1);
    let bell_state = match StateVectorSimulator::new().run_state(&bell) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let result = match StateVectorSimulator::new().run(&bell, 512) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover peak iteration: {peak}");
    println!("grover peak delta: {peak_delta}");
    println!(
        "qft probability imbalance: {:.6}",
        qft_probability_imbalance(&forward, &inverse)
    );
    println!(
        "qft reversible pair: {}",
        qft_pair_is_reversible(&forward, &inverse, &input)
    );
    println!(
        "bell parameterized gates: {}",
        bell.parameterized_gate_count()
    );
    println!(
        "bell collision probability: {:.6}",
        bell_state.collision_probability()
    );
    println!(
        "bell participation ratio: {:.6}",
        bell_state.participation_ratio()
    );
    println!("bell shannon entropy: {:.6} bits", result.shannon_entropy());
}

fn run_verify_m2_next(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: verify-m2-next needs at least 1 qubit");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let peak = match grover_first_success_iteration(qubits, &oracle, 2 * qubits + 2, 0.9) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let step = match grover_peak_step(qubits, &oracle, 2 * qubits + 2) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let targets = (0..qubits).collect::<Vec<_>>();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let input = match StateVector::zero(qubits) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover first success >= 0.9: {peak:?}");
    println!("grover peak step: {step:.6}");
    println!(
        "qft reversible pair: {}",
        qft_pair_is_reversible(&forward, &inverse, &input)
    );
    println!(
        "qft probability imbalance: {:.6}",
        qft_probability_imbalance(&forward, &inverse)
    );
    println!(
        "ghz ideal shannon entropy: {:.6} bits",
        state.shannon_entropy()
    );
    println!(
        "ghz collision probability: {:.6}",
        state.collision_probability()
    );
    println!(
        "ghz inverse max probability: {:.6}",
        state.inverse_max_probability()
    );
    println!("ghz multi-qubit gates: {}", ghz.multi_qubit_gate_count());
    println!(
        "ghz multi-qubit gate ratio: {:.6}",
        ghz.multi_qubit_gate_ratio()
    );
}

fn run(circuit: &Circuit, shots: usize) {
    match StateVectorSimulator::new().run(circuit, shots) {
        Ok(result) => {
            println!("circuit:\n{}\n", circuit.diagram());
            println!("probabilities: {:?}", result.probabilities());
            println!("histogram: {}", result.histogram_json(circuit.qubits()));
            println!("counts:");
            for (state, count) in result.sorted_counts() {
                println!("  {state:0width$b}: {count}", width = circuit.qubits());
            }
            let max = result
                .sorted_counts()
                .iter()
                .map(|(_, c)| *c)
                .max()
                .unwrap_or(1);
            println!("histogram:");
            for (state, count) in result.sorted_counts() {
                let filled = (count * 32) / max;
                let bar: String = std::iter::repeat('#').take(filled).collect();
                println!(
                    "  {state:0width$b} | {bar:<32} {count}",
                    width = circuit.qubits()
                );
            }
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn print_usage() {
    println!(
        "QSimX CLI\n\nUsage:\n  qsimx-cli bell [shots]\n  qsimx-cli ghz [qubits] [shots]\n  qsimx-cli grover [qubits] [shots]\n  qsimx-cli qft [qubits] [shots]\n  qsimx-cli qasm\n  qsimx-cli standard [qubits]\n  qsimx-cli inspect [qubits]\n  qsimx-cli verify-m2 [qubits]\n  qsimx-cli verify-m2-plus [qubits]\n  qsimx-cli verify-m2-final [qubits]\n  qsimx-cli verify-m2-compact [qubits]\n  qsimx-cli verify-m2-report [qubits]\n  qsimx-cli verify-m2-cond [qubits]\n  qsimx-cli verify-m2-dist [qubits]\n  qsimx-cli verify-m2-balance [qubits]\n  qsimx-cli verify-m2-parity [qubits]\n  qsimx-cli verify-m2-diagnostics [qubits]\n  qsimx-cli verify-m2-next [qubits]\n  qsimx-cli verify-m2-advanced [qubits]\n  qsimx-cli verify-m2-shape [qubits]\n  qsimx-cli verify-m2-curve [qubits]\n  qsimx-cli verify-m2-motion [qubits]\n  qsimx-cli verify-m2-stats [qubits]\n  qsimx-cli verify-m2-profile [qubits]\n  qsimx-cli verify-m2-signature [qubits]\n  qsimx-cli verify-m2-fingerprint [qubits]\n  qsimx-cli version\n"
    );
}

fn run_verify_m2_advanced(qubits: usize) {
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let area = match grover_success_area(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut angle_circuit = Circuit::new(qubits);
    angle_circuit.h(0).rx(0, 0.5).cp(0, 1 % qubits, 0.25);

    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let fidelity = match qft_gate_fidelity(&forward, qubits) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover success area: {area:.6}");
    println!("ghz entropy per qubit: {:.6}", state.entropy_per_qubit());
    println!("ghz probability gini: {:.6}", state.probability_gini());
    println!(
        "angle circuit parameterized ratio: {:.6}",
        angle_circuit.parameterized_gate_ratio()
    );
    println!("qft gate fidelity: {fidelity:.6}");
}

fn run_verify_m2_shape(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let width = match grover_peak_width(qubits, &oracle, 6, 0.5) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover peak width: {width}");
    println!("qft phase balance: {:.6}", qft_phase_balance(&forward));
    println!("ghz active qubits: {}", ghz.distinct_touching_qubit_count());
    println!("ghz max probability: {:.6}", state.max_probability());
}

fn run_verify_m2_curve(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let drop = match grover_peak_drop(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover peak drop: {drop:.6}");
    println!(
        "qft adjacent phase gap: {:.6}",
        qft_adjacent_phase_gap(&forward)
    );
    println!(
        "ghz odd parity probability: {:.6}",
        state.odd_parity_probability()
    );
    println!("ghz single-qubit gates: {}", ghz.single_qubit_gate_count());
}

fn run_verify_m2_motion(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let path_length = match grover_success_path_length(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover success path length: {path_length:.6}");
    println!(
        "qft adjacent phase step: {:.6}",
        qft_adjacent_phase_step(&forward)
    );
    println!("ghz tail probability: {:.6}", state.tail_probability());
    println!("ghz two-qubit gates: {}", ghz.gates_with_arity(2));
}

fn run_verify_m2_stats(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let variance = match grover_success_variance(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let ascending_steps = match grover_ascending_steps(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let inverse = match inverse_qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover success variance: {variance:.6}");
    println!("grover ascending steps: {ascending_steps}");
    println!(
        "qft uniform deviation: {:.6}",
        qft_uniform_deviation(&forward)
    );
    println!(
        "qft inverse distribution gap: {:.6}",
        qft_inverse_distribution_gap(&forward, &inverse)
    );
    println!(
        "ghz probability smoothness: {:.6}",
        state.probability_smoothness()
    );
    println!("ghz effective outcomes: {}", state.effective_outcomes());
    println!("ghz simple gates: {}", ghz.one_and_two_qubit_gate_count());
}

fn run_verify_m2_profile(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let skewness = match grover_success_skewness(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };

    println!("qubits: {qubits}");
    println!("grover success skewness: {skewness:.6}");
    println!(
        "qft phase circular variance: {:.6}",
        qft_phase_circular_variance(&forward)
    );
    println!("ghz peak concentration: {:.6}", state.peak_concentration());
    println!("ghz gate width entropy: {:.6}", ghz.gate_width_entropy());
}

fn run_verify_m2_signature(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let kurtosis = match grover_success_kurtosis(qubits, &oracle, 6) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    println!("qubits: {qubits}");
    println!("grover success kurtosis: {kurtosis:.6}");
    println!(
        "qft phase resultant length: {:.6}",
        qft_phase_resultant_length(&forward)
    );
    println!(
        "ghz uniformity divergence: {:.6}",
        state.uniformity_divergence()
    );
    println!("ghz controlled gates: {}", ghz.controlled_gate_count());
}

fn run_verify_m2_fingerprint(qubits: usize) {
    if qubits == 0 {
        eprintln!("error: qubits must be at least 1");
        return;
    }
    let oracle = GroverOracle::marked((1usize << qubits) - 1);
    let width_ratio = match grover_peak_width_ratio(qubits, &oracle, 6, 0.5) {
        Ok(value) => value,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let targets: Vec<usize> = (0..qubits).collect();
    let forward = match qft(&targets, qubits) {
        Ok(result) => result,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    let mut ghz = Circuit::new(qubits);
    ghz.h(0);
    for target in 1..qubits {
        ghz.cnot(0, target);
    }
    let state = match StateVectorSimulator::new().run_state(&ghz) {
        Ok(state) => state,
        Err(err) => {
            eprintln!("error: {err}");
            return;
        }
    };
    println!("qubits: {qubits}");
    println!("grover peak width ratio: {width_ratio:.6}");
    println!(
        "qft phase circular skewness: {:.6}",
        qft_phase_circular_skewness(&forward)
    );
    println!(
        "ghz active uniformity divergence: {:.6}",
        state.active_uniformity_divergence()
    );
    println!(
        "ghz entangling gate density: {:.6}",
        ghz.entangling_gate_density()
    );
}
