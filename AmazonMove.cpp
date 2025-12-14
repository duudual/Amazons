#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include<vector>
#include<queue>
#include<algorithm>
#include<climits>
#include<cmath>


#define GRIDSIZE 8
#define OBSTACLE 2
#define judge_black 0
#define judge_white 1
#define grid_black 1
#define grid_white -1

using namespace std;

int currBotColor; // 我所执子颜色（1为黑，-1为白，棋盘状态亦同）
int gridInfo[GRIDSIZE][GRIDSIZE] = { 0 }; // 先x后y，记录棋盘状态
int dx[] = { -1,-1,-1,0,0,1,1,1 };
int dy[] = { -1,0,1,-1,1,-1,0,1 };

// 判断是否在地图内
inline bool inMap(int x, int y)
{
	if (x < 0 || x >= GRIDSIZE || y < 0 || y >= GRIDSIZE)
		return false;
	return true;
}

struct Node{
	int id;
	int x,y;
	Node(int id_,int x_,int y_):id(id_),x(x_),y(y_){}
};

// 走法结构
struct Move {
	int x0, y0; // 起始位置
	int x1, y1; // 目标位置
	int x2, y2; // 障碍位置
	Move(int x0_, int y0_, int x1_, int y1_, int x2_, int y2_) 
		: x0(x0_), y0(y0_), x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
};

// // 在坐标处落子，检查是否合法或模拟落子
// bool ProcStep(int x0, int y0, int x1, int y1, int x2, int y2, int color, bool check_only)
// {
// 	if ((!inMap(x0, y0)) || (!inMap(x1, y1)) || (!inMap(x2, y2)))
// 		return false;
// 	if (gridInfo[x0][y0] != color || gridInfo[x1][y1] != 0)
// 		return false;
// 	if ((gridInfo[x2][y2] != 0) && !(x2 == x0 && y2 == y0))
// 		return false;
// 	if (!check_only)
// 	{
// 		gridInfo[x0][y0] = 0;
// 		gridInfo[x1][y1] = color;
// 		gridInfo[x2][y2] = OBSTACLE;
// 	}
// 	return true;
// }



// // 检查两个点是否可以直接按照皇后走法到达
// bool canReach(int x1, int y1, int x2, int y2) {
// 	if (x1 == x2 && y1 == y2) return true;
	
// 	int dx_dir = 0, dy_dir = 0;
// 	if (x2 != x1) dx_dir = (x2 > x1) ? 1 : -1;
// 	if (y2 != y1) dy_dir = (y2 > y1) ? 1 : -1;
	
// 	if (dx_dir != 0 && dy_dir != 0) {
// 		// 对角线
// 		if (abs(x2 - x1) != abs(y2 - y1)) return false;
// 	} else if (dx_dir == 0 && dy_dir == 0) {
// 		return false;
// 	}
	
// 	// 检查路径上的所有格子是否为空
// 	int steps = max(abs(x2 - x1), abs(y2 - y1));
// 	for (int i = 1; i < steps; ++i) {
// 		int x = x1 + dx_dir * i;
// 		int y = y1 + dy_dir * i;
// 		if (gridInfo[x][y] != 0) return false;
// 	}
// 	return true;
// }

// // queen为起点构建图，然后运行Dijkstra计算最短距离
// vector<int> dijstra(const vector<Node>& squares, const Node& queen) {
// 	const int INF = 1e5;
// 	int n = GRIDSIZE*GRIDSIZE; // squares + queen起点
// 	vector<int> dist(n, INF);
	
// 	// 构建图：节点0是queen起点，节点1到n-1是squares
// 	// 使用邻接表，边权都是1
// 	vector<vector<int>> g(n);
	
// 	// 构建边：queen-squares
// 	for (int i = 0; i < squares.size(); ++i) {
// 		if (canReach(queen.x, queen.y, squares[i].x, squares[i].y)) {
// 			g[0].push_back(i + 1);
// 			g[i + 1].push_back(0);
// 		}
// 	}
	
// 	// 构建边：squares之间
// 	for (int i = 0; i < squares.size(); ++i) {
// 		for (int j = i + 1; j < squares.size(); ++j) {
// 			if (canReach(squares[i].x, squares[i].y, squares[j].x, squares[j].y)) {
// 				g[i + 1].push_back(j + 1);
// 				g[j + 1].push_back(i + 1);
// 			}
// 		}
// 	}
	
// 	// Dijkstra算法
// 	dist[0] = 0; // queen起点距离为0
// 	using P = pair<int, int>; // dist, index
// 	priority_queue<P, vector<P>, greater<P>> pq;
// 	pq.push({ 0, 0 });
	
// 	while (!pq.empty()) {
// 		auto [d, u] = pq.top();
// 		pq.pop();
// 		if (d != dist[u]) continue;
// 		for (int v : g[u]) {
// 			if (dist[v] > d + 1) { // 边权都是1
// 				dist[v] = d + 1;
// 				pq.push({ dist[v], v });
// 			}
// 		}
// 	}
	
// 	// 返回距离数组，索引0是queen起点，索引1到n-1对应squares
// 	return dist;
// }

// //局面判断函数,可到达的空地数目越多,到达所需要的步数越短，得分越高
// //可直接按照皇后走法走到的两点形成edge,确定起点,可构建无向非负图.用Dijsktra寻找最短路.
// double judegeNow(int color){
// 	vector<Node> squares;
// 	vector<Node> queens;
// 	for (int i = 0; i < GRIDSIZE; ++i) {
// 		for (int j = 0; j < GRIDSIZE; ++j) {
// 			if(gridInfo[i][j]==0){
// 				squares.push_back(Node(i*GRIDSIZE+j,i,j));
// 			}
// 			if(gridInfo[i][j]==color){
// 				queens.push_back(Node(i*GRIDSIZE+j,i,j));
// 			}
// 		}
// 	}

// 	const int INF = 1e5;
// 	vector<int> glodis(squares.size(),INF);

// 	for(auto queen:queens){
// 		// 以当前皇后为起点构建图，计算到所有squares的最短距离
// 		vector<int> dist = dijstra(squares, queen);
// 		for (int i = 1; i < dist.size(); ++i) {
// 			if (dist[i] != INF) {
// 				glodis[i-1] = min(glodis[i-1], dist[i]); //glodis没有起点,i-1
// 				// glodis[i] += dist[i];// 可选：所有的最短距离相加,需要将起始距离改为0
// 			}
// 		}
// 	}

// 	double score = 0;
// 	for(auto dis:glodis){
// 		if(dis!=INF){
// 			score += 1.0/dis;
// 		}
// 	}
// 	// 对局结束
// 	if(score == 0){
// 		score = INF;
// 	}
// 	return score;
// }

// // 复用示例代码
// vector<Move> getAllMoves(int color) {
// 	vector<Move> moves;
// 	for (int i = 0; i < GRIDSIZE; ++i) {
// 		for (int j = 0; j < GRIDSIZE; ++j) {
// 			if (gridInfo[i][j] != color) continue;
// 			for (int k = 0; k < 8; ++k) {
// 				for (int delta1 = 1; delta1 < GRIDSIZE; delta1++) {
// 					int xx = i + dx[k] * delta1;
// 					int yy = j + dy[k] * delta1;
// 					if (gridInfo[xx][yy] != 0 || !inMap(xx, yy))
// 						break;
// 					for (int l = 0; l < 8; ++l) {
// 						for (int delta2 = 1; delta2 < GRIDSIZE; delta2++) {
// 							int xxx = xx + dx[l] * delta2;
// 							int yyy = yy + dy[l] * delta2;
// 							if (!inMap(xxx, yyy))
// 								break;
// 							if (gridInfo[xxx][yyy] != 0 && !(i == xxx && j == yyy))
// 								break;
// 							if (ProcStep(i, j, xx, yy, xxx, yyy, color, true)) {
// 								moves.push_back(Move(i, j, xx, yy, xxx, yyy));
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return moves;
// }

void copyBoard(int dst[GRIDSIZE][GRIDSIZE], int src[GRIDSIZE][GRIDSIZE]) {
    for (int i = 0; i < GRIDSIZE; ++i) {
        for (int j = 0; j < GRIDSIZE; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

// 在指定棋盘上检查两个点是否可以直接按照皇后走法到达
bool canReachOnBoard(int board[GRIDSIZE][GRIDSIZE], int x1, int y1, int x2, int y2) {
	if (x1 == x2 && y1 == y2) return true;
	
	int dx_dir = 0, dy_dir = 0;
	if (x2 != x1) dx_dir = (x2 > x1) ? 1 : -1;
	if (y2 != y1) dy_dir = (y2 > y1) ? 1 : -1;
	
	if (dx_dir != 0 && dy_dir != 0) {
		// 对角线
		if (abs(x2 - x1) != abs(y2 - y1)) return false;
	} else if (dx_dir == 0 && dy_dir == 0) {
		return false;
	}
	
	// 检查路径上的所有格子是否为空
	int steps = max(abs(x2 - x1), abs(y2 - y1));
	for (int i = 1; i < steps; ++i) {
		int x = x1 + dx_dir * i;
		int y = y1 + dy_dir * i;
		if (board[x][y] != 0) return false;
	}
	return true;
}

// 在指定棋盘上执行走法
bool ProcStepOnBoard(int board[GRIDSIZE][GRIDSIZE], int x0, int y0, int x1, int y1, int x2, int y2, int color) {
	if ((!inMap(x0, y0)) || (!inMap(x1, y1)) || (!inMap(x2, y2)))
		return false;
	if (board[x0][y0] != color || board[x1][y1] != 0)
		return false;
	if ((board[x2][y2] != 0) && !(x2 == x0 && y2 == y0))
		return false;
	board[x0][y0] = 0;
	board[x1][y1] = color;
	board[x2][y2] = OBSTACLE;
	return true;
}

// 在指定棋盘上获取所有合法走法
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
							// 检查走法是否合法
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

// 在指定棋盘上运行Dijkstra
vector<int> dijstraOnBoard(int board[GRIDSIZE][GRIDSIZE], const vector<Node>& squares, const Node& queen) {
	const int INF = 1e5;
	int n = GRIDSIZE*GRIDSIZE;
	vector<int> dist(n, INF);
	
	vector<vector<int>> g(n);
	
	// 构建边：queen-squares
	for (int i = 0; i < squares.size(); ++i) {
		if (canReachOnBoard(board, queen.x, queen.y, squares[i].x, squares[i].y)) {
			g[0].push_back(i + 1);
			g[i + 1].push_back(0);
		}
	}
	
	// 构建边：squares之间
	for (int i = 0; i < squares.size(); ++i) {
		for (int j = i + 1; j < squares.size(); ++j) {
			if (canReachOnBoard(board, squares[i].x, squares[i].y, squares[j].x, squares[j].y)) {
				g[i + 1].push_back(j + 1);
				g[j + 1].push_back(i + 1);
			}
		}
	}
	
	// Dijkstra算法
	dist[0] = 0;
	using P = pair<int, int>;
	priority_queue<P, vector<P>, greater<P>> pq;
	pq.push({ 0, 0 });
	
	while (!pq.empty()) {
		auto top = pq.top();
		int d = top.first;
		int u = top.second;
		pq.pop();
		if (d != dist[u]) continue;
		for (int v : g[u]) {
			if (dist[v] > d + 1) {
				dist[v] = d + 1;
				pq.push({ dist[v], v });
			}
		}
	}
	
	return dist;
}

// 在指定棋盘上评估局面
double judegeNowOnBoard(int board[GRIDSIZE][GRIDSIZE], int color) {
	vector<Node> squares;
	vector<Node> queens;
	for (int i = 0; i < GRIDSIZE; ++i) {
		for (int j = 0; j < GRIDSIZE; ++j) {
			if(board[i][j] == 0){
				squares.push_back(Node(i*GRIDSIZE+j,i,j));
			}
			if(board[i][j] == color){
				queens.push_back(Node(i*GRIDSIZE+j,i,j));
			}
		}
	}

	const int INF = 1e5;
	vector<int> glodis(squares.size(), INF);

	for(auto queen:queens){
		vector<int> dist = dijstraOnBoard(board, squares, queen);
		for (int i = 1; i < dist.size(); ++i) {
			if (dist[i] != INF) {
				glodis[i-1] = min(glodis[i-1], dist[i]);
			}
		}
	}

	double score = 0;
	for(auto dis:glodis){
		if(dis != INF){
			score += 1.0/dis;
		}
	}
	if(score == 0){
		score = INF;
	}
	return score;
}

// < -------- MC --------- >//
struct MCTSNode {
    // 当前节点所在局面（下完自己的上一手）
    int board[GRIDSIZE][GRIDSIZE];

    Move move;           // 本节点对应的走法（根节点无效）
    int color;           // 当前节点是哪个颜色行动（根节点为 currBotColor）
    
    double W;            // 胜利得分累积
    int N;               // 访问次数
    bool fullyExpanded;  // 是否所有子节点已扩展
    
    vector<MCTSNode*> children;
    MCTSNode* parent;

    MCTSNode(int b[GRIDSIZE][GRIDSIZE], int c, Move m, MCTSNode* p)
        : move(m), color(c), W(0), N(0), fullyExpanded(false), parent(p)
    {
        for (int i = 0; i < GRIDSIZE; i++)
            for (int j = 0; j < GRIDSIZE; j++)
                board[i][j] = b[i][j];
    }
};
double UCT(MCTSNode* parent, MCTSNode* child, double C = 1.414) {
    if (child->N == 0) return 1e9; // 优先探索未访问节点
    return child->W / child->N + C * sqrt(log(parent->N + 1) / child->N);
}
MCTSNode* selectByUCT(MCTSNode* node) {
    while (!node->children.empty()) {
        MCTSNode* best = nullptr;
        double bestUCT = -1e10;
        for (auto ch : node->children) {
            double u = UCT(node, ch);
            if (u > bestUCT) {
                bestUCT = u;
                best = ch;
            }
        }
        node = best;
    }
    return node;
}
MCTSNode* expand(MCTSNode* node) {
    vector<Move> moves = getAllMovesOnBoard(node->board, node->color);
    if (moves.empty()) return node; // 终局
    
    // 如果节点还没有子节点，全扩展
    if (node->children.empty()) {
        for (auto &mv : moves) {
            int newBoard[GRIDSIZE][GRIDSIZE];
            copyBoard(newBoard, node->board);
            ProcStepOnBoard(newBoard, mv.x0, mv.y0, mv.x1, mv.y1, mv.x2, mv.y2, node->color);

            node->children.push_back(new MCTSNode(newBoard, -node->color, mv, node));
        }
    }

    // 返回一个未访问过的子节点
    for (auto ch : node->children) {
        if (ch->N == 0) return ch;
    }

    // 如果都访问过，则认为已完全展开
    node->fullyExpanded = true;
    return node;
}
// 快速评估一个走法的质量（用于模拟阶段的选择）
double quickEvaluateMove(int board[GRIDSIZE][GRIDSIZE], const Move& m, int color) {
    int testBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(testBoard, board);
    ProcStepOnBoard(testBoard, m.x0, m.y0, m.x1, m.y1, m.x2, m.y2, color);
    
    double myScore = judegeNowOnBoard(testBoard, color);
    double oppScore = judegeNowOnBoard(testBoard, -color);
    
    return myScore - oppScore;
}

// 基于评估值的加权随机选择
Move selectMoveByWeight(const vector<Move>& moves, int board[GRIDSIZE][GRIDSIZE], int color, double temperature = 1.0) {
    if (moves.empty()) return Move(-1, -1, -1, -1, -1, -1);
    
    vector<double> scores;
    double minScore = 1e10, maxScore = -1e10;
    
    // 计算每个走法的评估值
    for (const auto& m : moves) {
        double score = quickEvaluateMove(board, m, color);
        scores.push_back(score);
        minScore = min(minScore, score);
        maxScore = max(maxScore, score);
    }
    
    // 归一化并应用softmax（带温度参数）
    vector<double> weights;
    double sum = 0;
    for (double s : scores) {
        // 归一化到[0, 1]，然后应用exp
        double normalized = (maxScore - minScore > 1e-6) ? 
            (s - minScore) / (maxScore - minScore) : 0.5;
        double weight = exp(normalized / temperature);
        weights.push_back(weight);
        sum += weight;
    }
    
    // 加权随机选择
    double r = (double)rand() / RAND_MAX * sum;
    double cumsum = 0;
    for (int i = 0; i < moves.size(); i++) {
        cumsum += weights[i];
        if (r <= cumsum) {
            return moves[i];
        }
    }
    
    return moves[moves.size() - 1];
}

double simulate(int board[GRIDSIZE][GRIDSIZE], int rootColor, int currentColor, int maxDepth = 10) {
    int simBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(simBoard, board);

    for (int d = 0; d < maxDepth; d++) {
        vector<Move> mv = getAllMovesOnBoard(simBoard, currentColor);
        if (mv.empty()) break;
        // 使用基于评估的加权随机选择，温度参数随深度增加（后期更随机）
        double temp = 0.5 + d * 0.1; // 温度从0.5逐渐增加到2.5
        Move m = selectMoveByWeight(mv, simBoard, currentColor, temp);
        ProcStepOnBoard(simBoard, m.x0, m.y0, m.x1, m.y1, m.x2, m.y2, currentColor);
        currentColor = -currentColor;
    }

    double myScore = judegeNowOnBoard(simBoard, rootColor);
    double oppScore = judegeNowOnBoard(simBoard, -rootColor);

    return myScore - oppScore;
}
void backpropagate(MCTSNode* node, double value) {
    while (node != nullptr) {
        node->N++;
        node->W += value; // value >0 对 rootColor 好
        value = -value;   // 视角反转
        node = node->parent;
    }
}
Move monteCarloSearch(int color, int ITER) {
    int rootBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(rootBoard, gridInfo);

    // 根节点（其 move 无意义）
    MCTSNode* root = new MCTSNode(rootBoard, color, Move(-1, -1, -1, -1, -1, -1), nullptr);

    for (int i = 0; i < ITER; i++) {
        // 1) Selection
        MCTSNode* selected = selectByUCT(root);

        // 2) Expansion
        MCTSNode* leaf = expand(selected);

        // 3) Simulation
        double value = simulate(leaf->board, color, leaf->color);

        // 4) Backpropagation
        backpropagate(leaf, value);
    }

    // 最终选择访问次数最多的子节点
    MCTSNode* best = nullptr;
    int bestN = -1;
    for (auto ch : root->children) {
        if (ch->N > bestN) {
            bestN = ch->N;
            best = ch;
        }
    }

    if (best == nullptr) return Move(-1,-1,-1,-1,-1,-1);
    return best->move;
}
// < -------- MC --------- >//



int main()
{
	int x0, y0, x1, y1, x2, y2;

	// 初始化棋盘
	gridInfo[0][(GRIDSIZE - 1) / 3] = gridInfo[(GRIDSIZE - 1) / 3][0]
		= gridInfo[GRIDSIZE - 1 - ((GRIDSIZE - 1) / 3)][0]
		= gridInfo[GRIDSIZE - 1][(GRIDSIZE - 1) / 3] = grid_black;
	gridInfo[0][GRIDSIZE - 1 - ((GRIDSIZE - 1) / 3)] = gridInfo[(GRIDSIZE - 1) / 3][GRIDSIZE - 1]
		= gridInfo[GRIDSIZE - 1 - ((GRIDSIZE - 1) / 3)][GRIDSIZE - 1]
		= gridInfo[GRIDSIZE - 1][GRIDSIZE - 1 - ((GRIDSIZE - 1) / 3)] = grid_white;


	int turnID;
	cin >> turnID;

	// 读入到当前回合为止，自己和对手的所有行动，从而把局面恢复到当前回合
	currBotColor = grid_white; // 先假设自己是白方
	for (int i = 0; i < turnID; i++)
	{
		// 根据这些输入输出逐渐恢复状态到当前回合

		// 首先是对手行动
		cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2;
		if (x0 == -1)
			currBotColor = grid_black; // 第一回合收到坐标是-1, -1，说明我是黑方
		else
			// ProcStep(x0, y0, x1, y1, x2, y2, -currBotColor, false); // 模拟对方落子
			ProcStepOnBoard(gridInfo,x0, y0, x1, y1, x2, y2, -currBotColor); // 模拟对方落子

																	// 然后是自己当时的行动
																	// 对手行动总比自己行动多一个
		if (i < turnID - 1)
		{
			cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2;
			if (x0 >= 0)
				// ProcStep(x0, y0, x1, y1, x2, y2, currBotColor, false); // 模拟己方落子
				ProcStepOnBoard(gridInfo, x0, y0, x1, y1, x2, y2, currBotColor); // 模拟己方落子
		}
	}

	// 做出决策（你只需修改以下部分）

	// 使用蒙特卡洛模拟选择最优走法
	srand(time(0));
	Move bestMove = monteCarloSearch(currBotColor, 30); // 每个走法模拟30次
	
	int startX, startY, resultX, resultY, obstacleX, obstacleY;
	if (bestMove.x0 >= 0) {
		startX = bestMove.x0;
		startY = bestMove.y0;
		resultX = bestMove.x1;
		resultY = bestMove.y1;
		obstacleX = bestMove.x2;
		obstacleY = bestMove.y2;
	}
	else
	{
		startX = -1;
		startY = -1;
		resultX = -1;
		resultY = -1;
		obstacleX = -1;
		obstacleY = -1;
	}

	// 决策结束，输出结果（你只需修改以上部分）
	cout << startX << ' ' << startY << ' ' << resultX << ' ' << resultY << ' ' << obstacleX << ' ' << obstacleY << endl;
	return 0;
}