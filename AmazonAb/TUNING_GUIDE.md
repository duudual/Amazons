# Alpha-Beta剪枝参数调优指南

## 快速调整参数

### 1. 搜索深度调整

在文件开头找到：
```cpp
int searchDepth = 4; // Alpha-Beta搜索深度
```

**推荐设置：**
- 快速决策（<1秒）：深度 2-3
- 平衡性能（1-5秒）：深度 4-5
- 强力搜索（5-30秒）：深度 6-7
- 极限搜索（>30秒）：深度 8+

**注意**：每增加1层深度，搜索时间约增加10-50倍（取决于剪枝效率）

### 2. 评估函数权重调整

在`evaluateBoard`函数中找到：
```cpp
int totalScore = territoryScore * 100 +      // 领地控制
                 mobilityScore * 5 +          // 移动能力
                 centerScore * 3 +            // 中心控制
                 positionScore * 2 +          // 棋子位置
                 (myCloseTerritory - oppCloseTerritory) * 20;  // 近距离领地
```

**调优建议：**

#### 开局阶段（前10回合）
- 增加`mobilityScore`权重到 8-10
- 增加`centerScore`权重到 5-8
- 保持领地控制权重最高

#### 中盘阶段（中间阶段）
- 使用默认权重平衡各项指标
- 领地控制 100，移动能力 5，中心控制 3

#### 残局阶段（走法<20）
- 大幅增加`territoryScore`权重到 150-200
- 降低移动能力权重到 2-3
- 重点争夺剩余空位

### 3. 置换表大小调整

```cpp
const int TT_SIZE = 1048576; // 2^20 entries
```

**内存占用估算：**
- 当前设置（2^20）：约 40-50 MB
- 增大到 2^22：约 160-200 MB（更好的缓存）
- 减小到 2^18：约 10-13 MB（节省内存）

**建议：**
- 比赛环境：使用 2^20 或更大
- 内存受限：使用 2^18
- 服务器运行：可以使用 2^24

### 4. 针对不同对手的策略

#### 对抗进攻型对手
- 增加`mobilityScore`权重
- 增加搜索深度
- 重视保持行动自由

#### 对抗防守型对手
- 增加`territoryScore`权重
- 增加`centerScore`权重
- 尽早占据中心区域

#### 对抗随机对手
- 降低搜索深度（节省时间）
- 使用默认权重即可

## 性能监控

### GUI版本
在`alphaBetaSearch`函数中，取消注释调试输出：
```cpp
cout << "[Alpha-Beta] Depth " << depth << ": score = " << score << endl;
```

这会显示每层的搜索分数，帮助你了解：
- 搜索是否在改善
- 评估函数是否合理
- 是否有异常情况

### 置换表命中率
添加统计代码：
```cpp
// 在置换表查询处添加
static int ttHits = 0, ttQueries = 0;
ttQueries++;
if(ttEntry->hash == hash && ttEntry->depth >= depth) {
    ttHits++;
    if(ttQueries % 1000 == 0) {
        cout << "TT Hit Rate: " << (100.0 * ttHits / ttQueries) << "%" << endl;
    }
}
```

**良好的命中率：**
- 30-50%：正常水平
- 50-70%：优秀
- <20%：考虑增大置换表

## 测试建议

### 1. 基准测试
设置固定的测试局面，记录：
- 搜索时间
- 搜索节点数
- 置换表命中率
- 最终胜率

### 2. 自我对战
让不同参数版本互相对战：
```bash
# 深度3 vs 深度4
# 权重A vs 权重B
```

### 3. 局面分析
保存关键局面，分析AI的决策：
- 是否抓住机会
- 是否避免失误
- 评估分数是否合理

## 高级优化

### 1. 动态深度调整
```cpp
int calculateDepth(int board[GRIDSIZE][GRIDSIZE]) {
    int moveCount = getAllMovesOnBoard(board, currBotColor).size();
    if(moveCount < 10) return 6;  // 残局加深
    if(moveCount > 100) return 3; // 开局浅搜
    return 4;  // 中盘正常
}
```

### 2. 时间管理
```cpp
#include <chrono>
auto startTime = chrono::high_resolution_clock::now();
// ... 搜索代码 ...
auto elapsed = chrono::duration_cast<chrono::milliseconds>(
    chrono::high_resolution_clock::now() - startTime
).count();
if(elapsed > timeLimit) break; // 超时退出
```

### 3. 开局优化
```cpp
// 对于前几步，使用预定义的好位置
if(turnID <= 2) {
    // 返回中心区域的移动
}
```

## 调试技巧

### 1. 评估函数验证
```cpp
// 测试已知好/坏局面的评分
int testBoard[GRIDSIZE][GRIDSIZE] = {...};
int score = evaluateBoard(testBoard, grid_black);
cout << "Test position score: " << score << endl;
```

### 2. 搜索树可视化
```cpp
// 在alphaBetaSearch中添加
string indent(depth * 2, ' ');
cout << indent << "Depth " << depth << ", score " << score << endl;
```

### 3. 走法质量检查
```cpp
// 记录每次的最佳走法
ofstream log("moves.log", ios::app);
log << "Turn " << turnID << ": " << bestMove.x0 << "," << bestMove.y0 
    << " -> " << bestMove.x1 << "," << bestMove.y1 << endl;
```

## 常见问题

**Q: 搜索太慢怎么办？**
A: 降低搜索深度，检查移动排序是否有效，增大置换表。

**Q: AI走出明显失误？**
A: 检查评估函数权重，增加搜索深度，分析该局面的评分。

**Q: 内存占用过大？**
A: 减小置换表大小，或使用更小的哈希函数。

**Q: 如何提高胜率？**
A: 增加搜索深度，优化评估函数，改进移动排序，收集对手数据。
