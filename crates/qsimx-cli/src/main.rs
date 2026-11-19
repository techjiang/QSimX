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

fn run(circuit: &Circuit, shots: usize) {
    match StateVectorSimulator::new().run(circuit, shots) {
        Ok(result) => {
            println!("circuit:\n{}\n", circuit.diagram());
            println!("probabilities: {:?}", result.probabilities());
            println!("counts:");
            for (state, count) in result.sorted_counts() {
                println!("  {state:0width$b}: {count}", width = circuit.qubits());
            }
        }
        Err(err) => eprintln!("error: {err}"),
    }
}

fn print_usage() {
    println!(
        "QSimX CLI\n\nUsage:\n  qsimx-cli bell [shots]\n  qsimx-cli ghz [qubits] [shots]\n  qsimx-cli qasm\n"
    );
}
