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
cargo run -p qsimx-cli -- analyze 3
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

## 第十轮迭代：M2 验证与稳定性

本轮在已有 M2 算法基础上做 5 次增量更新，重点是把 Grover 的迭代数与验证输出做实：

1. 新增 `grover_iterations()` / `recommended_grover_iterations()`，使用 `floor(pi/4 * sqrt(N/M))` 避免小搜索空间过冲，并约束非法输入返回 1 次。
2. CLI 的 `grover-m2` 改用推荐迭代数，避免 4 qubit 概率被过冲拉低。
3. 新增 `Circuit::w_state()`，快速构造 active target 上的均匀叠加态。
4. 新增 `StateVector::significant_basis_count()` 与 `SimulationResult::effective_outcomes()`，便于统计显著基态与有效测量输出数。
5. CLI 新增 `verify-m2`，集中输出 Grover / QFT / W 态指标，并补齐回归测试。

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- verify-m2 3
```

本机已安装 Rust 1.100 与 C 工具链后真实执行：`cargo test --workspace` 23/23 通过；`verify-m2 3` 输出 Grover 2 次迭代、标记态概率 0.500000、QFT 7 门、2 个显著振幅、QFT 态归一化。

## 第十一轮迭代

1. `StateVector::entanglement_entropy()` 计算单量子比特纠缠熵，`purity()` 计算约化态纯度。
2. `SimulationResult::entropy_stability()` 估计测量分布的熵稳定性。
3. `qft_roundtrip_fidelity()` 验证 QFT 与逆 QFT 的状态保真度。
4. `grover_success_curve()` 输出从 0 次开始的 Grover 成功率曲线。
5. CLI 新增 `verify-m2-plus`，集中展示纠缠熵、纯度、QFT 往返保真度和 Grover 曲线。

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- verify-m2-plus 3
```

## 第十二轮迭代

1. 新增 `StateVector::total_variation_distance()`，计算两个量子态测量分布的总变差距离。
2. 新增 `Circuit::state_purity()`，报告量子态最大振幅概率作为纯度快照。
3. 新增 `qft_sparsity()`，统计 QFT 输出中超过阈值的显著模式数。
4. CLI 新增 `verify-m2-final`，集中输出 TVD、纯度、QFT 模式数、往返保真度与 Grover 成功率。
5. 补充 3 项回归测试，覆盖新 API 的边界与数值行为。

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- verify-m2-final 3
```

真实执行：`cargo test --workspace` 30/30 通过；`verify-m2-final 3` 输出 TVD 0.500000、GHZ 最大振幅概率 0.500000、QFT 显著模式 8、往返保真度 1.000000、Grover 迭代 2 次且标记态概率 0.500000。


## 第十三轮迭代

1. 新增 `StateVector::probability_span()`，报告基态概率极差。
2. 新增 `SimulationResult::histogram_text()`，把测量直方图渲染为可复用文本。
3. 新增 `Circuit::distinct_touching_qubits()` 与 `Circuit::add_layered()`，支持按活跃量子比特补层。
4. 新增 `qft_gate_density()` 与 `grover_amplification_gain()`，补齐 M2 的小型分析指标。
5. CLI 新增 `verify-m2-compact`，集中输出本轮指标与文本直方图，并补充 4 项回归测试。

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- verify-m2-compact 3
```

## 第十四轮迭代

1. 新增 `StateVector::basis_parity()`，返回基态的按位奇偶校验。
2. 新增 `SimulationResult::parity_counts()`，按奇偶校验汇总测量计数。
3. 新增 `Circuit::w_state()` 的边界校验，并保持其均匀叠加语义。
4. 新增 `grover_success_cumulative()` 与 `qft_average_probability()`，补充算法累计指标。
5. CLI 新增 `verify-m2-report`，集中输出本轮指标，并补充 5 项回归测试。

验证命令：

```bash
cargo test --workspace
cargo run -p qsimx-cli -- verify-m2-report 3
```

## 第十五轮迭代

1. 新增 `StateVector::conditional_bit_probability()`，计算固定一个量子比特条件下的目标位概率。
2. 新增 `StateVector::basis_bitwise_and()`，汇总基态中所有置位比特位置的按位与。
3. 新增 `qft_gate_symmetry()`，比较 QFT 与逆 QFT 的门数对称性。
4. CLI 新增 `verify-m2-cond`，集中输出条件概率、基态位运算与 QFT 对称性指标。
5. 补充 3 项回归测试，覆盖新 API 的边界条件与数值结果。

验证命令：

```bash
cargo test --workspace
cargo run -q -p qsimx-cli -- verify-m2-cond 3
```

## 第十六轮迭代

1. 新增 `StateVector::jensen_shannon_divergence()`，衡量两个量子态测量分布的 JS 散度。
2. 新增 `Circuit::max_gates_per_qubit()`，统计单条量子线上的最大门数量。
3. 新增 `qft_density_gap()`，比较 QFT 与逆 QFT 的门密度差异。
4. CLI 新增 `verify-m2-dist`，集中输出 Grover 迭代数、JS 散度、TVD、往返保真度与电路负载。
5. 补充 4 项回归测试，覆盖新 API 的数值与边界条件。

验证命令：

```bash
cargo test --workspace
cargo run -q -p qsimx-cli -- verify-m2-dist 3
```

## 第十七轮迭代

1. 新增 `StateVector::even_parity_probability()`，计算偶校验基态概率质量。
2. 新增 `StateVector::renyi2_entropy()`，计算二阶 Rényi 熵。
3. 新增 `Circuit::min_gates_per_touched_qubit()`，统计活跃量子线上的最小门负载。
4. 新增 `qft_gate_balance()`，衡量 QFT 与逆 QFT 的门数平衡度。
5. CLI 新增 `verify-m2-balance`，并补充 4 项回归测试。

验证命令：

```bash
cargo test --workspace
cargo run -q -p qsimx-cli -- verify-m2-balance 3
```

## 第十八轮迭代

1. 新增 `StateVector::odd_parity_probability()` 与 `max_probability_step()`，补齐概率分布奇偶与跳变分析。
2. 新增 `Circuit::active_qubit_count()`，确认活跃量子线数量。
3. 新增 `grover_recommended_gain()`，量化推荐迭代次数带来的成功率提升。
4. 新增 `qft_density_balance()`，衡量 QFT 与逆 QFT 的平均门密度。
5. CLI 新增 `verify-m2-parity [qubits]`，并补充 4 项回归测试。

验证命令：

```bash
cargo test --workspace
cargo run -q -p qsimx-cli -- verify-m2-parity 3
```

## 第十九轮迭代

1. 新增 `StateVector::collision_probability()`，计算两次独立测量落入同一基态的概率。
2. 新增 `StateVector::participation_ratio()`，给出有效测量结果数量。
3. 新增 `Circuit::parameterized_gate_count()`，统计依赖连续角度参数的门数量。
4. CLI 新增 `verify-m2-diagnostics`，集中展示 Grover 峰值、QFT 对称性与 Bell 态统计。
5. 补充 2 项回归测试，覆盖基态/均匀态与参数门计数。

验证命令：

```bash
cargo test --workspace
cargo run -q -p qsimx-cli -- verify-m2-diagnostics 3
```
