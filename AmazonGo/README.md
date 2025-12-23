# 亚马逊棋价值网络训练系统

基于蒙特卡洛模拟的价值网络实现，用于亚马逊棋游戏的局面评估。

## 文件说明

- `value_network.cpp` - MLP价值网络实现，包含前向传播、反向传播和训练功能
- `data_generator.cpp` - 数据生成pipeline，使用文件IO记录更大depth的蒙特卡洛树
- `integration_test.cpp` - 集成测试程序，展示完整的使用流程
- `Amazon.cpp` - 原有的亚马逊棋GUI程序（来自父目录）

## 编译方法

```bash
# 编译价值网络
g++ -std=c++17 -o value_network value_network.cpp

# 编译数据生成器
g++ -std=c++17 -o data_generator data_generator.cpp

# 编译集成测试
g++ -std=c++17 -o integration_test integration_test.cpp
```

## 使用方法

### 1. 生成训练数据

```bash
# 生成训练数据
# 参数: 游戏数量 最大深度 每个位置的模拟次数
./data_generator 10 50 100

# 示例：生成10个游戏，最大深度50，每个位置模拟100次
```

这会在当前目录生成：
- `training_data.bin` - 训练数据文件
- `mc_tree_*.node` - 蒙特卡洛树节点文件

### 2. 训练价值网络

修改 `integration_test.cpp` 中的参数，或者创建自定义的训练程序：

```cpp
// 创建网络: 64(输入) -> 128 -> 64 -> 1(输出)
vector<int> layer_sizes = {64, 128, 64, 1};
ValueNetwork network(layer_sizes);

// 加载训练数据
DataGenerator generator("mc_tree");
auto training_data = generator.load_training_data("training_data.bin");

// 训练网络
train_network(network, training_data, 100, 0.001);  // 100个epoch，学习率0.001

// 保存模型
network.save_model("trained_model.model");
```

### 3. 使用训练好的模型

```cpp
// 加载模型
ValueNetwork network({64, 128, 64, 1});
network.load_model("trained_model.model");

// 评估局面
int board[8][8] = { /* 棋盘状态 */ };
double value = network.forward(board);  // 获取局面价值估计
```

## 网络架构

- **输入**: 64个节点 (8x8棋盘flatten)
  - 1: 黑子位置
  - -1: 白子位置
  - 0: 空位
  - 0.5: 障碍位置

- **隐藏层**: 可配置的ReLU激活层

- **输出**: 1个节点，使用tanh激活，输出范围[-1, 1]
  - 正值: 对黑方有利
  - 负值: 对白方有利

## 数据生成流程

1. **初始化**: 从标准亚马逊棋初始局面开始

2. **蒙特卡洛模拟**:
   - 对每个可能走法进行多次随机模拟
   - 模拟到最大深度或终局
   - 使用quick_score函数评估叶子节点

3. **反向传播**:
   - 将叶子节点的得分反向传播更新路径上的所有节点
   - 每个节点的mc_score是其所有子节点的平均得分

4. **数据收集**:
   - 收集树中所有节点的局面状态和对应的mc_score
   - 保存为二进制训练数据文件

## 文件IO设计

为了支持大深度蒙特卡洛树，系统使用文件存储：

- **节点文件**: `mc_tree_{node_id}.node`
  - 包含完整的节点状态信息
  - 支持大树的持久化存储

- **训练数据**: `training_data.bin`
  - 二进制格式存储局面和对应的价值
  - 高效的批量加载

- **模型文件**: `*.model`
  - 保存训练好的网络权重和偏置
  - 文本格式，便于调试和跨平台使用

## 性能优化

- **缓存机制**: 最近使用的节点保存在内存中
- **批量训练**: 支持mini-batch梯度下降
- **文件缓冲**: 使用二进制IO提高读写效率
- **内存限制**: 可配置最大缓存大小

## 扩展建议

1. **并行化**: 数据生成和网络训练可以并行化
2. **数据增强**: 通过旋转和镜像扩展训练数据
3. **更深的网络**: 尝试更深的网络架构或残差连接
4. **正则化**: 添加dropout或L2正则化防止过拟合
5. **学习率调度**: 使用学习率衰减或自适应优化器

## 与原有系统的集成

训练好的价值网络可以集成到原有的 `Amazon.cpp` MCTS算法中：

1. 在MCTS的选择阶段使用价值网络进行局面评估
2. 在叶子节点扩展时使用价值网络初始化节点价值
3. 结合传统MCTS和价值网络的混合搜索算法

这样可以在保持MCTS优势的同时，利用神经网络的学习能力提高搜索效率。



