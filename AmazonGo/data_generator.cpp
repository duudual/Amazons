#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <random>
#include <ctime>
#include <cstring>

#define GRIDSIZE 8
#define CELL 80
#define OBSTACLE 2
#define judge_black 0
#define judge_white 1
#define grid_black 1
#define grid_white -1

using namespace std;

// 随机数生成器
random_device rd;
mt19937 gen(rd());

// 走法结构
struct Move {
	int x0, y0; // 起始位置
	int x1, y1; // 目标位置
	int x2, y2; // 障碍位置
	Move(int x0_, int y0_, int x1_, int y1_, int x2_, int y2_)
		: x0(x0_), y0(y0_), x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
};

// 训练数据样本
struct TrainingSample {
    int board[GRIDSIZE][GRIDSIZE];
    double mc_score;  // 蒙特卡洛分数

    TrainingSample() {
        memset(board, 0, sizeof(board));
        mc_score = 0.0;
    }

    void copy_board(const int src[GRIDSIZE][GRIDSIZE]) {
        for (int i = 0; i < GRIDSIZE; ++i) {
            for (int j = 0; j < GRIDSIZE; ++j) {
                board[i][j] = src[i][j];
            }
        }
    }
};

// 蒙特卡洛节点（用于文件存储的简化版本）
struct MCNode {
    string node_id;  // 节点唯一标识符
    int board[GRIDSIZE][GRIDSIZE];
    string parent_id;
    Move move;       // 从父节点到此节点的走法
    int color;       // 当前行动颜色
    double W;        // 累积胜分
    int N;          // 访问次数
    vector<string> children_ids;

    MCNode() : move(-1,-1,-1,-1,-1,-1), W(0.0), N(0), color(0) {
        memset(board, 0, sizeof(board));
    }
};

// 全局变量
int dx[] = {-1,-1,-1,0,0,1,1,1};
int dy[] = {-1,0,1,-1,1,-1,0,1};

// 工具函数
inline bool inMap(int x, int y) {
    if (x < 0 || x >= GRIDSIZE || y < 0 || y >= GRIDSIZE)
        return false;
    return true;
}

void copyBoard(int dst[GRIDSIZE][GRIDSIZE], int src[GRIDSIZE][GRIDSIZE]) {
    for (int i = 0; i < GRIDSIZE; ++i) {
        for (int j = 0; j < GRIDSIZE; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

bool canReachOnBoard(int board[GRIDSIZE][GRIDSIZE], int x1, int y1, int x2, int y2) {
    if (x1 == x2 && y1 == y2) return true;

    int dx_dir = 0, dy_dir = 0;
    if (x2 != x1) dx_dir = (x2 > x1) ? 1 : -1;
    if (y2 != y1) dy_dir = (y2 > y1) ? 1 : -1;

    if (dx_dir != 0 && dy_dir != 0) {
        if (abs(x2 - x1) != abs(y2 - y1)) return false;
    } else if (dx_dir == 0 && dy_dir == 0) {
        return false;
    }

    int steps = max(abs(x2 - x1), abs(y2 - y1));
    for (int i = 1; i < steps; ++i) {
        int x = x1 + dx_dir * i;
        int y = y1 + dy_dir * i;
        if (board[x][y] != 0) return false;
    }
    return true;
}

bool ProcStepOnBoard(int board[GRIDSIZE][GRIDSIZE], int x0, int y0, int x1, int y1, int x2, int y2, int color) {
    if ((!inMap(x0, y0)) || (!inMap(x1, y1)) || (!inMap(x2, y2)))
        return false;
    if (board[x0][y0] != color || board[x1][y1] != 0)
        return false;
    if ((board[x2][y2] != 0) && !(x2 == x0 && y2 == y0))
        return false;
    if (!canReachOnBoard(board, x0, y0, x1, y1))
        return false;
    board[x0][y0] = 0;
    board[x1][y1] = color;
    board[x2][y2] = OBSTACLE;
    return true;
}

vector<Move> getAllMovesOnBoard(int board[GRIDSIZE][GRIDSIZE], int color) {
    vector<Move> moves;
    for (int i = 0; i < GRIDSIZE; ++i) {
        for (int j = 0; j < GRIDSIZE; ++j) {
            if (board[i][j] != color) continue;
            for (int k = 0; k < 8; ++k) {
                for (int delta1 = 1; delta1 < GRIDSIZE; delta1++) {
                    int xx = i + dx[k] * delta1;
                    int yy = j + dy[k] * delta1;
                    if (board[xx][yy] != 0 || !inMap(xx, yy))
                        break;
                    for (int l = 0; l < 8; ++l) {
                        for (int delta2 = 1; delta2 < GRIDSIZE; delta2++) {
                            int xxx = xx + dx[l] * delta2;
                            int yyy = yy + dy[l] * delta2;
                            if (!inMap(xxx, yyy))
                                break;
                            if (board[xxx][yyy] != 0 && !(i == xxx && j == yyy))
                                break;
                            if (board[i][j] == color && board[xx][yy] == 0) {
                                if ((board[xxx][yyy] == 0) || (i == xxx && j == yyy)) {
                                    moves.push_back(Move(i, j, xx, yy, xxx, yyy));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return moves;
}

// BFS计算领土
int blackDist[GRIDSIZE][GRIDSIZE];
int whiteDist[GRIDSIZE][GRIDSIZE];

void bfsTerritory(int board[GRIDSIZE][GRIDSIZE], int color, int distMap[GRIDSIZE][GRIDSIZE]) {
    for(int i=0; i<GRIDSIZE; i++)
        for(int j=0; j<GRIDSIZE; j++)
            distMap[i][j] = 9999;

    queue<pair<int, int>> q;

    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == color) {
                distMap[i][j] = 0;
                q.push({i, j});
            }
        }
    }

    while(!q.empty()) {
        auto [cx, cy] = q.front();
        q.pop();

        for(int k=0; k<8; k++) {
            int nx = cx + dx[k];
            int ny = cy + dy[k];

            if(inMap(nx, ny) && board[nx][ny] == 0) {
                if(distMap[nx][ny] > distMap[cx][cy] + 1) {
                    distMap[nx][ny] = distMap[cx][cy] + 1;
                    q.push({nx, ny});
                }
            }
        }
    }
}

double quickScore(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    bfsTerritory(board, grid_black, blackDist);
    bfsTerritory(board, grid_white, whiteDist);

    double score = 0;
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == 0) {
                int b = blackDist[i][j];
                int w = whiteDist[i][j];
                if(b < w) score -= 1.0;
                else if (w < b) score += 1.0;
                else if (b!=9999) score += 0.0;
            }
        }
    }
    if (myColor == grid_white) return score;
    else return -score;
}

// 文件操作类 - 用于管理蒙特卡洛树节点
class MCTreeStorage {
private:
    string base_filename;
    unordered_map<string, MCNode> node_cache;
    int max_cache_size;

public:
    MCTreeStorage(const string& filename, int cache_size = 1000)
        : base_filename(filename), max_cache_size(cache_size) {}

    // 生成节点ID
    string generate_node_id(const int board[GRIDSIZE][GRIDSIZE], int color) {
        string id = to_string(color) + "_";
        for (int i = 0; i < GRIDSIZE; ++i) {
            for (int j = 0; j < GRIDSIZE; ++j) {
                id += to_string(board[i][j]) + ",";
            }
        }
        return id;
    }

    // 保存节点到文件
    void save_node(const MCNode& node) {
        string filename = base_filename + "_" + node.node_id + ".node";
        ofstream file(filename, ios::binary);

        // 写入基本信息
        size_t id_len = node.node_id.length();
        file.write(reinterpret_cast<char*>(&id_len), sizeof(size_t));
        file.write(node.node_id.c_str(), id_len);

        size_t parent_len = node.parent_id.length();
        file.write(reinterpret_cast<char*>(&parent_len), sizeof(size_t));
        file.write(node.parent_id.c_str(), parent_len);

        // 写入棋盘
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.board[0][0])), sizeof(int) * GRIDSIZE * GRIDSIZE);

        // 写入走法
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.x0)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.y0)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.x1)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.y1)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.x2)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.move.y2)), sizeof(int));

        // 写入统计信息
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.color)), sizeof(int));
        file.write(reinterpret_cast<char*>(const_cast<double*>(&node.W)), sizeof(double));
        file.write(reinterpret_cast<char*>(const_cast<int*>(&node.N)), sizeof(int));

        // 写入子节点数量和ID
        size_t children_count = node.children_ids.size();
        file.write(reinterpret_cast<char*>(&children_count), sizeof(size_t));
        for (const auto& child_id : node.children_ids) {
            size_t child_len = child_id.length();
            file.write(reinterpret_cast<char*>(&child_len), sizeof(size_t));
            file.write(child_id.c_str(), child_len);
        }

        file.close();
    }

    // 从文件加载节点
    bool load_node(const string& node_id, MCNode& node) {
        // 首先检查缓存
        if (node_cache.find(node_id) != node_cache.end()) {
            node = node_cache[node_id];
            return true;
        }

        string filename = base_filename + "_" + node_id + ".node";
        ifstream file(filename, ios::binary);

        if (!file.is_open()) return false;

        // 读取基本信息
        size_t id_len;
        file.read(reinterpret_cast<char*>(&id_len), sizeof(size_t));
        node.node_id.resize(id_len);
        file.read(&node.node_id[0], id_len);

        size_t parent_len;
        file.read(reinterpret_cast<char*>(&parent_len), sizeof(size_t));
        node.parent_id.resize(parent_len);
        file.read(&node.parent_id[0], parent_len);

        // 读取棋盘
        file.read(reinterpret_cast<char*>(&node.board[0][0]), sizeof(int) * GRIDSIZE * GRIDSIZE);

        // 读取走法
        file.read(reinterpret_cast<char*>(&node.move.x0), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.move.y0), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.move.x1), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.move.y1), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.move.x2), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.move.y2), sizeof(int));

        // 读取统计信息
        file.read(reinterpret_cast<char*>(&node.color), sizeof(int));
        file.read(reinterpret_cast<char*>(&node.W), sizeof(double));
        file.read(reinterpret_cast<char*>(&node.N), sizeof(int));

        // 读取子节点
        size_t children_count;
        file.read(reinterpret_cast<char*>(&children_count), sizeof(size_t));
        node.children_ids.resize(children_count);
        for (size_t i = 0; i < children_count; ++i) {
            size_t child_len;
            file.read(reinterpret_cast<char*>(&child_len), sizeof(size_t));
            node.children_ids[i].resize(child_len);
            file.read(&node.children_ids[i][0], child_len);
        }

        file.close();

        // 添加到缓存
        if (node_cache.size() < max_cache_size) {
            node_cache[node_id] = node;
        }

        return true;
    }

    // 更新节点统计信息
    void update_node_stats(const string& node_id, double delta_W, int delta_N) {
        MCNode node;
        if (load_node(node_id, node)) {
            node.W += delta_W;
            node.N += delta_N;

            // 更新缓存
            node_cache[node_id] = node;

            // 保存到文件
            save_node(node);
        }
    }

    // 添加子节点
    void add_child(const string& parent_id, const string& child_id) {
        MCNode parent;
        if (load_node(parent_id, parent)) {
            parent.children_ids.push_back(child_id);
            node_cache[parent_id] = parent;
            save_node(parent);
        }
    }

};

// 数据生成器
class DataGenerator {
private:
    MCTreeStorage storage;
    vector<TrainingSample> training_data;
    int max_depth;
    int simulations_per_position;

public:
    DataGenerator(const string& storage_prefix, int max_d = 100, int sims = 1000)
        : storage(storage_prefix), max_depth(max_d), simulations_per_position(sims) {}

    // 生成初始棋盘
    void initialize_board(int board[GRIDSIZE][GRIDSIZE]) {
        memset(board, 0, sizeof(int) * GRIDSIZE * GRIDSIZE);
        // 标准亚马逊棋初始布局
        board[0][2] = board[2][0] = board[5][0] = board[7][2] = grid_black;
        board[0][5] = board[2][7] = board[5][7] = board[7][5] = grid_white;
    }

    // 模拟单次游戏直到最大深度
    double simulate_game(int board[GRIDSIZE][GRIDSIZE], int current_color, int max_depth) {
        int sim_board[GRIDSIZE][GRIDSIZE];
        copyBoard(sim_board, board);

        for (int d = 0; d < max_depth; ++d) {
            vector<Move> moves = getAllMovesOnBoard(sim_board, current_color);
            if (moves.empty()) {
                // 当前玩家无法移动，游戏结束
                return (current_color == grid_black) ? -1.0 : 1.0;
            }

            // 随机选择一个走法
            uniform_int_distribution<> dist(0, moves.size() - 1);
            Move m = moves[dist(gen)];

            ProcStepOnBoard(sim_board, m.x0, m.y0, m.x1, m.y1, m.x2, m.y2, current_color);
            current_color = -current_color;
        }

        // 达到最大深度，使用快速评分
        return quickScore(sim_board, grid_black);
    }

    // 生成训练数据
    void generate_training_data(int num_games) {
        cout << "Generating training data with " << num_games << " games..." << endl;

        for (int game = 0; game < num_games; ++game) {
            cout << "Game " << (game + 1) << "/" << num_games << endl;

            int root_board[GRIDSIZE][GRIDSIZE];
            initialize_board(root_board);

            // 为每个游戏创建根节点
            MCNode root_node;
            copyBoard(root_node.board, root_board);
            root_node.color = grid_black;
            root_node.node_id = "root_" + to_string(game);
            root_node.parent_id = "";
            root_node.move = Move(-1, -1, -1, -1, -1, -1);

            storage.save_node(root_node);

            // 从根节点开始进行蒙特卡洛模拟
            expand_and_simulate(root_node.node_id, max_depth);

            // 从这个游戏中收集训练数据
            collect_training_data(root_node.node_id);
        }

        cout << "Generated " << training_data.size() << " training samples." << endl;
    }

    // 扩展并模拟节点
    void expand_and_simulate(const string& node_id, int remaining_depth) {
        if (remaining_depth <= 0) return;

        MCNode node;
        if (!storage.load_node(node_id, node)) return;

        // 获取所有可能的走法
        vector<Move> moves = getAllMovesOnBoard(node.board, node.color);
        if (moves.empty()) return;

        // 为每个走法创建子节点并进行模拟
        for (const auto& move : moves) {
            int child_board[GRIDSIZE][GRIDSIZE];
            copyBoard(child_board, node.board);
            ProcStepOnBoard(child_board, move.x0, move.y0, move.x1, move.y1, move.x2, move.y2, node.color);

            MCNode child_node;
            copyBoard(child_node.board, child_board);
            child_node.color = -node.color;
            child_node.node_id = node_id + "_child_" + to_string(node.children_ids.size());
            child_node.parent_id = node_id;
            child_node.move = move;

            // 进行多次模拟
            double total_score = 0.0;
            for (int s = 0; s < simulations_per_position; ++s) {
                double score = simulate_game(child_board, child_node.color, remaining_depth - 1);
                total_score += score;
            }

            child_node.W = total_score;
            child_node.N = simulations_per_position;

            storage.save_node(child_node);
            storage.add_child(node_id, child_node.node_id);

            // 递归扩展
            if (remaining_depth > 1) {
                expand_and_simulate(child_node.node_id, remaining_depth - 1);
            }
        }

        // 反向传播：更新父节点的统计信息
        backpropagate_scores(node_id);
    }

    // 反向传播分数
    void backpropagate_scores(const string& node_id) {
        MCNode node;
        if (!storage.load_node(node_id, node)) return;

        if (node.children_ids.empty()) return;

        // 计算平均分数
        double total_score = 0.0;
        int total_visits = 0;

        for (const auto& child_id : node.children_ids) {
            MCNode child;
            if (storage.load_node(child_id, child)) {
                total_score += child.W;
                total_visits += child.N;
            }
        }

        // 更新当前节点
        node.W = total_score / node.children_ids.size();
        node.N = total_visits;

        storage.save_node(node);

        // 递归更新父节点
        if (!node.parent_id.empty()) {
            backpropagate_scores(node.parent_id);
        }
    }

    // 收集训练数据
    void collect_training_data(const string& node_id) {
        MCNode node;
        if (!storage.load_node(node_id, node)) return;

        // 添加当前节点的数据
        TrainingSample sample;
        sample.copy_board(node.board);
        sample.mc_score = node.W / max(1, node.N);  // 平均分数
        training_data.push_back(sample);

        // 递归收集子节点的数据
        for (const auto& child_id : node.children_ids) {
            collect_training_data(child_id);
        }
    }

    // 保存训练数据到文件
    void save_training_data(const string& filename) {
        ofstream file(filename, ios::binary);

        size_t num_samples = training_data.size();
        file.write(reinterpret_cast<char*>(&num_samples), sizeof(size_t));

        for (const auto& sample : training_data) {
            file.write(reinterpret_cast<char*>(const_cast<int*>(&sample.board[0][0])), sizeof(int) * GRIDSIZE * GRIDSIZE);
            file.write(reinterpret_cast<char*>(const_cast<double*>(&sample.mc_score)), sizeof(double));
        }

        file.close();
        cout << "Training data saved to " << filename << endl;
    }

    // 加载训练数据
    vector<TrainingSample> load_training_data(const string& filename) {
        vector<TrainingSample> data;
        ifstream file(filename, ios::binary);

        if (!file.is_open()) {
            cerr << "Cannot open training data file: " << filename << endl;
            return data;
        }

        size_t num_samples;
        file.read(reinterpret_cast<char*>(&num_samples), sizeof(size_t));

        data.resize(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            file.read(reinterpret_cast<char*>(&data[i].board[0][0]), sizeof(int) * GRIDSIZE * GRIDSIZE);
            file.read(reinterpret_cast<char*>(&data[i].mc_score), sizeof(double));
        }

        file.close();
        cout << "Loaded " << data.size() << " training samples from " << filename << endl;
        return data;
    }

    // 获取训练数据
    const vector<TrainingSample>& get_training_data() const {
        return training_data;
    }
};

int main(int argc, char* argv[]) {
    srand(time(0));

    // 参数解析
    int num_games = 10;  // 默认10个游戏
    int max_depth = 50;  // 默认最大深度50
    int simulations = 100; // 默认每个位置100次模拟

    if (argc > 1) num_games = atoi(argv[1]);
    if (argc > 2) max_depth = atoi(argv[2]);
    if (argc > 3) simulations = atoi(argv[3]);

    cout << "Data Generator Configuration:" << endl;
    cout << "Number of games: " << num_games << endl;
    cout << "Max depth: " << max_depth << endl;
    cout << "Simulations per position: " << simulations << endl;

    // 创建数据生成器
    DataGenerator generator("mc_tree", max_depth, simulations);

    // 生成训练数据
    generator.generate_training_data(num_games);

    // 保存训练数据
    generator.save_training_data("training_data.bin");

    cout << "Data generation completed!" << endl;
    return 0;
}
