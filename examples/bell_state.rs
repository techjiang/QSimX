use qsimx::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut circuit = Circuit::new(2);
    circuit.h(0).cnot(0, 1);

    let result = StateVectorSimulator::new().run(&circuit, 1024)?;
    println!(" Bell state");
    println!("{}", circuit.diagram());
    println!("probabilities: {:?}", result.probabilities());
    for (state, count) in result.sorted_counts() {
        println!("|{state:02b}> x {count}");
    }
    Ok(())
}
