# QSimX · Quantum Computing Simulator

在笔记本上运行量子计算机。

QSimX 是一个用 Rust 实现的高性能量子计算模拟器，目标是让量子电路模拟、算法实验和可视化变得简单。

## 当前进度

当前仓库包含里程碑 M1 的可用骨架：

- 状态向量模拟器（最多 30 qubit，受内存限制）
- 标准量子门：H / X / Y / Z / S / T / RX / RY / RZ / CNOT / CZ / SWAP / Toffoli / 自定义 2×2 门
- 电路构建器与 QASM 导出
- 计算基测量、Bloch 向量、泡利期望值
- CLI 示例：Bell 态、GHZ 态
- 单元测试与集成测试

后续里程碑：M2（Grover/QFT + 可视化）、M3（噪声/GPU/MPS/Python 绑定）、M4（Web/WASM/导出/CLI 生态）。

## 快速开始

```bash
cargo test -p qsimx-core
cargo run -p qsimx-cli -- bell 1024
cargo run -p qsimx-cli -- ghz 3 1024
```

## Rust 示例

```rust
use qsimx::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut circuit = Circuit::new(3);
    circuit.h(0).cnot(0, 1).cnot(0, 2);

    let result = StateVectorSimulator::new().run(&circuit, 1024)?;
    println!("{:?}", result.probabilities());
    println!("{:?}", result.bloch_vector(0)?);
    println!("{}", circuit.to_qasm()?);
    Ok(())
}
```

## 许可证

Apache-2.0 © 2026 科技酱
