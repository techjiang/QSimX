use std::env;

use qsimx::prelude::{
    Circuit, Gate, GroverOracle, StateVector, StateVectorSimulator, grover_amplification_gain,
    grover_search, grover_success_curve, qft, qft_gate_density, qft_roundtrip_fidelity,
    qft_sparsity, recommended_grover_iterations,
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
    w.w_state(&targets);
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
        "QSimX CLI\n\nUsage:\n  qsimx-cli bell [shots]\n  qsimx-cli ghz [qubits] [shots]\n  qsimx-cli grover [qubits] [shots]\n  qsimx-cli qft [qubits] [shots]\n  qsimx-cli qasm\n  qsimx-cli standard [qubits]\n  qsimx-cli inspect [qubits]\n  qsimx-cli verify-m2 [qubits]\n  qsimx-cli verify-m2-plus [qubits]\n  qsimx-cli verify-m2-final [qubits]\n  qsimx-cli verify-m2-compact [qubits]\n  qsimx-cli version\n"
    );
}
