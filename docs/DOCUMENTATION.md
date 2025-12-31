# 项目文档

## 系统架构文档

### 整体架构

本系统采用模块化设计，主要包括以下几个模块：

```
┌─────────────────────────────────────────────────────────┐
│                    多目标跳频信号定位系统                    │
└─────────────────────────────────────────────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        ▼                   ▼                   ▼
   ┌─────────┐         ┌─────────┐         ┌─────────┐
   │信号生成  │         │信号处理  │         │定位分析  │
   └─────────┘         └─────────┘         └─────────┘
        │                   │                   │
        ├─ 跳频图案         ├─ 信道化滤波       ├─ TDOA估计
        ├─ 调制解调         ├─ 能量检测         ├─ 定位算法
        └─ 多目标合成       ├─ 信号分选         └─ 误差分析
                            └─ 特征提取
```

## 关键算法说明

### 1. 信道化滤波器设计

**多相滤波器组结构：**

- 原型滤波器：低通FIR
- 滤波器阶数：自适应计算
- 频率响应：矩形或升余弦
- 实现方式：FFT快速卷积

**核心代码位置：** 在各主程序的信道化部分

### 2. 能量检测算法

**自适应阈值设定：**

$$\text{Threshold} = \mu_{noise} + k \times \sigma_{noise}$$

其中：
- $\mu_{noise}$: 噪声均值估计
- $\sigma_{noise}$: 噪声标准差
- $k$: 阈值因子（通常3-5）

### 3. DBSCAN聚类分选

**参数设置：**
- Epsilon: 邻域半径，根据特征空间自适应
- MinPts: 最小点数，通常设为5-10

**特征空间：**
- 时间特征（归一化时间戳）
- 频率特征（信道索引）
- 幅度特征（归一化能量）
- DOA特征（到达角估计）
- 带宽特征
- 稳定性特征

### 4. TDOA定位算法

**GCC-PHAT时延估计：**

$$\tau = \arg\max_\tau R_{12}(\tau)$$

其中 $R_{12}(\tau)$ 是广义互相关函数。

**最小二乘定位：**

迭代求解目标位置 $(x, y)$ 使得：

$$\min \sum_{i=1}^{N} (d_i - d_1 - c\tau_i)^2$$

## API文档

### generate_multi_target_signals

```matlab
function RxSignals = generate_multi_target_signals(TargetPositions, TargetPowers, HopPatterns, StationPos, Fs, T_hop, NumHops, SymbolRates)
```

**功能：** 生成多目标跳频接收信号

**输入参数：**
- `TargetPositions`: 目标位置矩阵 [N×2]，单位：米
- `TargetPowers`: 目标功率向量 [N×1]
- `HopPatterns`: 跳频图案细胞数组 {N×1}
- `StationPos`: 观测站位置 [M×2]
- `Fs`: 采样率（Hz）
- `T_hop`: 跳频持续时间（秒）
- `NumHops`: 跳频次数
- `SymbolRates`: 符号速率向量 [N×1]

**输出：**
- `RxSignals`: 接收信号矩阵 [M×Samples]

### multi_target_signal_sorting

```matlab
function multi_target_signal_sorting(RxSignals, StationPos, Fs, T_hop)
```

**功能：** 多目标信号分选主函数

**输入参数：**
- `RxSignals`: 多站接收信号 [M×Samples]
- `StationPos`: 观测站位置 [M×2]
- `Fs`: 采样率（Hz）
- `T_hop`: 跳频持续时间（秒）

**输出：**
- 可视化结果（图形）
- 分选结果打印

### extract_pulse_features

```matlab
function features = extract_pulse_features(PulseList, ChannelOutputs, StationPos)
```

**功能：** 提取脉冲多维特征

**输入参数：**
- `PulseList`: 脉冲列表结构体数组
- `ChannelOutputs`: 信道输出细胞数组
- `StationPos`: 观测站位置

**输出：**
- `features`: 特征矩阵 [N×6]
  - 列1: 时间戳
  - 列2: 频率
  - 列3: 幅度
  - 列4: DOA
  - 列5: 带宽
  - 列6: 稳定性

## 使用示例

### 示例1：基础双目标定位

```matlab
% 设置参数
NumTargets = 2;
TargetPositions = [5000, 8000; 12000, 6000];
TargetPowers = [100, 150];

% 运行主程序
run('20251113/hop_freq_dpd_2_40m_5_targets_add_hop_cic.m')
```

### 示例2：五目标复杂场景

```matlab
% 修改参数
NumTargets = 5;
TargetPositions = [5000, 8000; 12000, 6000; 8000, 12000; 6000, 5000; 10000, 10000];

% 调整跳频图案
HopPatterns = { ... };  % 自定义5个跳频图案

% 运行仿真
```

## 性能优化建议

### 计算效率优化

1. **使用向量化操作**
   - 避免for循环处理大数组
   - 利用MATLAB内置向量函数

2. **预分配内存**
   ```matlab
   signals = zeros(NumStations, TotalSamples);
   ```

3. **并行处理**
   - 对独立信道并行处理
   - 使用parfor处理多站数据

### 精度改进

1. **提高采样率** - 更好的时延分辨率
2. **增加信道数** - 更精细的频率分辨率
3. **优化滤波器设计** - 减少旁瓣干扰

## 常见问题

### Q1: 分选效果不佳？

**解决方案：**
- 调整DBSCAN参数（Epsilon、MinPts）
- 增加特征维度权重
- 检查信噪比是否足够

### Q2: 定位误差过大？

**解决方案：**
- 增加观测站数量
- 优化站址几何布局
- 提高时延估计精度

### Q3: 运行速度慢？

**解决方案：**
- 减少数据长度
- 降低采样率（在满足奈奎斯特定理前提下）
- 使用并行计算

## 参考文献

1. Peterson, R. L., Ziemer, R. E., & Borth, D. E. (1995). Introduction to spread spectrum communications.
2. Poisel, R. A. (2011). Electronic warfare target location methods.
3. Torrieri, D. (2015). Principles of spread-spectrum communication systems.

## 版本历史

- **v1.0** (2025-12) - 初始版本
  - 基础多目标信号生成
  - 信道化滤波实现
  - TDOA定位功能

## 后续开发计划

- [ ] 添加实时数据接口
- [ ] 机器学习分类器
- [ ] GPU加速
- [ ] GUI界面开发
