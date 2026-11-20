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

## 第七轮迭代

1. 新增 `StateVector::amplitude_string()`，统一振幅文本格式。
2. 新增 `StateVector::cumulative_probabilities()` 与 `Gate::name()`。
3. 新增 `Circuit::depth()` 与 `Circuit::gate_counts()`。
4. 新增 `SimulationResult::shannon_entropy()`、`top_counts()`，以及 CLI `inspect` 命令。
5. 补齐回归测试与文档；`git diff --check` 无空白错误。

## 快速开始

```bash
cargo test -p qsimx-core
cargo run -p qsimx-cli -- bell 1024
cargo run -p qsimx-cli -- ghz 3 1024
cargo run -p qsimx-cli -- grover 3 1024
cargo run -p qsimx-cli -- qft 3 1024
cargo run -p qsimx-cli -- inspect 3
cargo run -p qsimx-cli -- summary9 2
cargo run -p qsimx-cli -- version
```

## 迭代记录

本轮连续 5 次增量更新：

1. `SimulationResult::histogram_json()` 导出 JSON 直方图
2. CLI 新增 `grover` / `qft` 示例电路，核心新增受控相位门 `CP`
3. `StateVector::fidelity()` 计算状态保真度及测试
4. CLI 输出终端文本直方图，新增 `version` 命令
5. 补齐 README 与迭代说明

## 第八轮迭代

1. `StateVector::zero_probability()` 计算单量子比特零态概率
2. `StateVector::phase()` 返回基态振幅相位
3. `Gate::qasm_name()` 提供稳定 QASM 门名
4. `Circuit::summary()` 汇总量子比特数、门数与深度
5. `SimulationResult::state_fidelity()` 计算结果态与参考态保真度

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- standard 3
```

## 第九轮迭代

1. `StateVector::one_probability()` 计算单量子比特一态概率
2. `Circuit::bell_pair()` 快速构造 Bell 态电路
3. `SimulationResult::most_likely()` 返回最可能测量结果及其频率
4. CLI 新增 `summary9 [qubits]` 命令，集中展示本轮指标
5. 补齐 3 项回归测试与 README 记录

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- summary9 2
```


## M2：Grover / QFT 完整实现

本轮先把开放 PR 收敛：PR #10 关闭并指向更高版本分支，PR #11 合并到当前基线。随后把 M2 的两个核心算法补成库级 API：

- `qft(targets, qubits)` 与 `inverse_qft(targets, qubits)` 返回真实电路和终态
- `grover_search(qubits, oracle, iterations)` 支持标记态 oracle 与自定义电路 oracle
- `GroverOracle::marked()` 构造单解 oracle，扩散算子由精确 QFT 组合生成
- CLI 增加 `grover-m2` / `qft-m2`，输出门数、标记态概率和电路图
- 集成测试覆盖 QFT 基态映射、QFT 逆变换、Grover 放大与非法输入

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- grover-m2 3
cargo run -p qsimx-cli -- qft-m2 3
```

当前环境若没有 Rust 工具链，上述命令需在 Rust 2024 环境中执行；截图只展示真实执行过的命令输出。

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
