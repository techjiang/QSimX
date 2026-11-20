use std::env;

use qsimx::prelude::{Circuit, StateVectorSimulator};

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
            run_grover(qubits, args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024));
        }
        "qft" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_qft(qubits, args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1024));
        }
        "standard" => run_standard(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3)),
        "inspect" => {
            let qubits = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
            run_inspect(qubits);
        }
        "summary9" => run_summary9(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2)),
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
            println!("zero probability q0: {:.6}", result.state().zero_probability(0).unwrap());
            println!("phase |0>: {:.6}", result.state().phase(0));
            println!("state fidelity: {:.6}", result.state_fidelity(&reference).unwrap());
            println!("gate names: {:?}", circuit.gates().iter().map(|gate| gate.name()).collect::<Vec<_>>());
            println!("qasm names: {:?}", circuit.gates().iter().map(|gate| gate.qasm_name()).collect::<Vec<_>>());
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
            println!("nonzero amplitudes: {}", result.state().amplitudes().iter().filter(|a| a.norm_sqr() > 1e-12).count());
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
            println!("one probability q0: {:.6}", result.state().one_probability(0).unwrap());
            println!("most likely: {:?}", result.most_likely());
            println!("qasm names: {:?}", circuit.gates().iter().map(|g| g.qasm_name()).collect::<Vec<_>>());
        }
        Err(err) => eprintln!("error: {err}"),
    }
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
            let max = result.sorted_counts().iter().map(|(_, c)| *c).max().unwrap_or(1);
            println!("histogram:");
            for (state, count) in result.sorted_counts() {
                let filled = (count * 32) / max;
                let bar: String = std::iter::repeat('#').take(filled).collect();
                println!("  {state:0width$b} | {bar:<32} {count}", width = circuit.qubits());
            }
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn print_usage() {
    println!(
        "QSimX CLI\n\nUsage:\n  qsimx-cli bell [shots]\n  qsimx-cli ghz [qubits] [shots]\n  qsimx-cli grover [qubits] [shots]\n  qsimx-cli qft [qubits] [shots]\n  qsimx-cli qasm\n  qsimx-cli standard [qubits]\n  qsimx-cli inspect [qubits]\n  qsimx-cli version\n"
    );
}
