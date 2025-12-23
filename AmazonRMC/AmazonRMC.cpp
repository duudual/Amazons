#define UNICODE
#define _UNICODE

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include<vector>
#include<queue>
#include<algorithm>
#include<climits>
#include<cmath>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <limits>

#define GRIDSIZE 8
#define CELL 80
#define OBSTACLE 2
#define judge_black 0
#define judge_white 1
#define grid_black 1
#define grid_white -1

using namespace std;

// 棋盘状态快照（用于存盘 & 复盘）
struct GameSnapshot {
    int board[GRIDSIZE][GRIDSIZE];
    int currentTurnColor;
};

int currBotColor; // 我所执子颜色（1为黑，-1为白，棋盘状态亦同）
int gridInfo[GRIDSIZE][GRIDSIZE] = { 0 }; // 先x后y，记录棋盘状态
int dx[] = { -1,-1,-1,0,0,1,1,1 };
int dy[] = { -1,0,1,-1,1,-1,0,1 };

int maxDepth = 40;//模拟最大深度
int iterations = 80000;//模拟次数
double c = 1.414;//探索参数

// 路径结构：存储从起点到终点的所有点
struct Path {
    vector<pair<int, int>> points;  // 路径上的所有点（包括起点和终点）
    Path() {}
    Path(vector<pair<int, int>> pts) : points(pts) {}
    bool contains(int x, int y) const {
        for(auto& p : points) {
            if(p.first == x && p.second == y) return true;
        }
        return false;
    }
};

// 游戏状态枚举
enum GameStatus {
    GAME_ONGOING = 0,    // 游戏进行中
    BLACK_WIN = 1,       // 黑方胜
    WHITE_WIN = 2,       // 白方胜
    DRAW = 3             // 平局
};


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
    Move() : x0(-1), y0(-1), x1(-1), y1(-1), x2(-1), y2(-1) {}
};


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
	if (!canReachOnBoard(board, x0, y0, x1, y1))
		return false;
	board[x0][y0] = 0;
	board[x1][y1] = color;
	board[x2][y2] = OBSTACLE;
	return true;
}

// 获取某个棋子可以移动到的所有路径（从起点到终点的所有点）
vector<Path> getValidMovePaths(int board[GRIDSIZE][GRIDSIZE], int x0, int y0) {
	vector<Path> paths;
	if (!inMap(x0, y0) || board[x0][y0] == 0) return paths;
	
	for (int k = 0; k < 8; ++k) {
		vector<pair<int, int>> pathPoints;
		pathPoints.push_back({x0, y0});  // 包含起点
		for (int delta = 1; delta < GRIDSIZE; delta++) {
			int x1 = x0 + dx[k] * delta;
			int y1 = y0 + dy[k] * delta;
			if (!inMap(x1, y1) || board[x1][y1] != 0)
				break;
			pathPoints.push_back({x1, y1});
			paths.push_back(Path(pathPoints));  // 每个路径都包含从起点到当前点的所有点
		}
	}
	return paths;
}

// 获取从某个位置可以放置障碍的所有路径（包括原位置）
// board应该是移动后的棋盘状态（原位置已空，新位置有棋子）
vector<Path> getValidObstaclePaths(int board[GRIDSIZE][GRIDSIZE], int x1, int y1, int x0, int y0) {
	vector<Path> paths;
	if (!inMap(x1, y1)) return paths;
	
	for (int k = 0; k < 8; ++k) {
		vector<pair<int, int>> pathPoints;
		pathPoints.push_back({x1, y1});  // 包含起点
		for (int delta = 1; delta < GRIDSIZE; delta++) {
			int x2 = x1 + dx[k] * delta;
			int y2 = y1 + dy[k] * delta;
			if (!inMap(x2, y2))
				break;
			if (board[x2][y2] != 0)
				break;
			pathPoints.push_back({x2, y2});
			paths.push_back(Path(pathPoints));  // 每个路径都包含从起点到当前点的所有点
		}
	}
	// 也可以放回原位置（如果原位置在可到达范围内且为空）
	if (inMap(x0, y0) && board[x0][y0] == 0 && canReachOnBoard(board, x1, y1, x0, y0)) {
		bool hasOriginalPath = false;
		for(auto& path : paths) {
			if(path.contains(x0, y0)) {
				hasOriginalPath = true;
				break;
			}
		}
		if(!hasOriginalPath) {
			vector<pair<int, int>> originalPath;
			// 计算从x1,y1到x0,y0的路径
			int dx_dir = 0, dy_dir = 0;
			if (x0 != x1) dx_dir = (x0 > x1) ? 1 : -1;
			if (y0 != y1) dy_dir = (y0 > y1) ? 1 : -1;
			int steps = max(abs(x0 - x1), abs(y0 - y1));
			for (int i = 0; i <= steps; ++i) {
				int x = x1 + dx_dir * i;
				int y = y1 + dy_dir * i;
				originalPath.push_back({x, y});
			}
			paths.push_back(Path(originalPath));
		}
	}
	return paths;
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

// 占地记录
int blackDist[GRIDSIZE][GRIDSIZE];
int whiteDist[GRIDSIZE][GRIDSIZE];

// bfs计算从特定颜色棋子出发，到达全图的最短步数
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

// 得分
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
    // 如果是白方(grid_white = -1)，正分好；如果是黑方(1)，负分好
    if (myColor == grid_white) return score; 
    else return -score;
}

// 计算黑方得分（从黑方视角）
double getBlackScore(int board[GRIDSIZE][GRIDSIZE]) {
    return quickScore(board, grid_black);
}

// 计算白方得分（从白方视角）
double getWhiteScore(int board[GRIDSIZE][GRIDSIZE]) {
    return quickScore(board, grid_white);
}



// 判断游戏是否结束
GameStatus checkGameStatus(int board[GRIDSIZE][GRIDSIZE]) {
    vector<Move> blackMoves = getAllMovesOnBoard(board, grid_black);
    vector<Move> whiteMoves = getAllMovesOnBoard(board, grid_white);
    
    bool blackCanMove = !blackMoves.empty();
    bool whiteCanMove = !whiteMoves.empty();
    
    if (!blackCanMove && !whiteCanMove) {
        // 双方都无法移动，比较得分
        double blackScore = getBlackScore(board);
        double whiteScore = getWhiteScore(board);
        if (blackScore > whiteScore) return BLACK_WIN;
        else if (whiteScore > blackScore) return WHITE_WIN;
        else return DRAW;
    } else if (!blackCanMove) {
        // 黑方无法移动，白方胜
        return WHITE_WIN;
    } else if (!whiteCanMove) {
        // 白方无法移动，黑方胜
        return BLACK_WIN;
    }
    
    return GAME_ONGOING;
}

// < -------- MCTS --------- >//
struct MCTSNode {
    Move move;                    // 该节点对应的走法
    int wins;                     // 获胜次数
    int visits;                   // 访问次数
    vector<MCTSNode*> children;   // 子节点列表
    MCTSNode* parent;             // 父节点指针
    int board[GRIDSIZE][GRIDSIZE]; // 当前节点的棋盘状态
    int currentPlayer;            // 当前玩家颜色

    MCTSNode(Move m, int boardState[GRIDSIZE][GRIDSIZE], int player, MCTSNode* p = nullptr)
        : move(m), wins(0), visits(0), parent(p), currentPlayer(player) {
        copyBoard(board, boardState);
    }

    MCTSNode(int boardState[GRIDSIZE][GRIDSIZE], int player)
        : move(), wins(0), visits(0), parent(nullptr), currentPlayer(player) {
        copyBoard(board, boardState);
    }

    ~MCTSNode() {
        for (auto child : children) {
            delete child;
        }
    }

    // 检查是否为叶节点（未完全扩展）
    bool isLeaf() const {
        return children.empty();
    }

    // 检查是否为完全展开的节点
    bool isFullyExpanded() {
        if (isLeaf()) return false;
        vector<Move> allMoves = getAllMovesOnBoard(board, currentPlayer);
        return children.size() == allMoves.size();
    }

    // 获取UCT值
    double getUCTValue(double c = 1.414) const {  // c是探索参数，通常用sqrt(2)
        if (visits == 0) return numeric_limits<double>::max();  // 未访问过的节点优先级最高

        double exploitation = static_cast<double>(wins) / visits;
        double exploration = c * sqrt(log(parent->visits) / visits);
        return exploitation + exploration;
    }
};

// 选择阶段：使用UCT算法从根节点选择到叶节点
MCTSNode* selectNode(MCTSNode* root, double c = 1.414) {
    MCTSNode* current = root;

    while (!current->isLeaf()) {
        double bestUCT = -numeric_limits<double>::max();
        MCTSNode* bestChild = nullptr;

        for (MCTSNode* child : current->children) {
            double uctValue = child->getUCTValue(c);
            if (uctValue > bestUCT) {
                bestUCT = uctValue;
                bestChild = child;
            }
        }

        if (bestChild == nullptr) break;  // 不应该发生，但为了安全
        current = bestChild;
    }

    return current;
}

// 扩展阶段：为叶节点添加一个子节点
MCTSNode* expandNode(MCTSNode* node) {
    // 获取当前节点的所有合法走法
    vector<Move> possibleMoves = getAllMovesOnBoard(node->board, node->currentPlayer);

    // 找出尚未扩展的走法
    vector<Move> unexpandedMoves;
    for (const Move& move : possibleMoves) {
        bool alreadyExpanded = false;
        for (MCTSNode* child : node->children) {
            if (child->move.x0 == move.x0 && child->move.y0 == move.y0 &&
                child->move.x1 == move.x1 && child->move.y1 == move.y1 &&
                child->move.x2 == move.x2 && child->move.y2 == move.y2) {
                alreadyExpanded = true;
                break;
            }
        }
        if (!alreadyExpanded) {
            unexpandedMoves.push_back(move);
        }
    }

    if (unexpandedMoves.empty()) {
        return node;  // 已经完全扩展
    }

    // 随机选择一个未扩展的走法进行扩展
    int randomIndex = rand() % unexpandedMoves.size();
    Move selectedMove = unexpandedMoves[randomIndex];

    // 创建新节点
    int newBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(newBoard, node->board);
    ProcStepOnBoard(newBoard, selectedMove.x0, selectedMove.y0,
                   selectedMove.x1, selectedMove.y1,
                   selectedMove.x2, selectedMove.y2, node->currentPlayer);

    MCTSNode* newNode = new MCTSNode(selectedMove, newBoard, -node->currentPlayer, node);
    node->children.push_back(newNode);

    return newNode;
}

// 模拟阶段：从指定节点开始随机模拟直到终局
double simulateGame(MCTSNode* node,int maxDepth) {
    int simBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(simBoard, node->board);
    int currentPlayer = node->currentPlayer;

    // 限制模拟的最大步数，避免无限循环
    int steps = 0;

    while (steps < maxDepth) {
        // 检查游戏是否结束
        GameStatus status = checkGameStatus(simBoard);
        if (status != GAME_ONGOING) {
            // 返回结果：从根玩家的视角来看
            if (status == BLACK_WIN) {
                return (node->parent->currentPlayer == grid_black) ? 1.0 : 0.0;
            } else if (status == WHITE_WIN) {
                return (node->parent->currentPlayer == grid_white) ? 1.0 : 0.0;
            } else { // DRAW
                return 0.5;
            }
        }

        // 获取当前玩家的所有合法走法
        vector<Move> moves = getAllMovesOnBoard(simBoard, currentPlayer);
        if (moves.empty()) {
            // 当前玩家无路可走，游戏结束
            GameStatus status = checkGameStatus(simBoard);
            if (status == BLACK_WIN) {
                return (node->parent->currentPlayer == grid_black) ? 1.0 : 0.0;
            } else if (status == WHITE_WIN) {
                return (node->parent->currentPlayer == grid_white) ? 1.0 : 0.0;
            } else {
                return 0.5;
            }
        }

        // 随机选择一个走法
        int randomIndex = rand() % moves.size();
        Move selectedMove = moves[randomIndex];

        // 执行走法
        ProcStepOnBoard(simBoard, selectedMove.x0, selectedMove.y0,
                       selectedMove.x1, selectedMove.y1,
                       selectedMove.x2, selectedMove.y2, currentPlayer);

        // 切换玩家
        currentPlayer = -currentPlayer;
        steps++;
    }

    // 达到最大步数，使用启发式评分
    double score = quickScore(simBoard, node->parent->currentPlayer);
    return (score > 0) ? 1.0 : (score < 0) ? 0.0 : 0.5;
}

// 回溯阶段：将模拟结果沿着路径回传
void backpropagate(MCTSNode* node, double result) {
    MCTSNode* current = node;
    while (current != nullptr) {
        current->visits++;
        if (result > 0.5) {  // 获胜
            current->wins++;
        } else if (result < 0.5) {  // 失败
            // wins不增加
        } else {  // 平局
            current->wins += 0.5;
        }

        // 切换视角（因为下一层是从对手的视角看）
        result = 1.0 - result;
        current = current->parent;
    }
}

// 标准蒙特卡洛树搜索主函数
Move monteCarloSearch(int color, int iterations, int c) {
    int rootBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(rootBoard, gridInfo);

    // 创建根节点
    MCTSNode* root = new MCTSNode(rootBoard, color);

    // 检查是否有合法走法
    vector<Move> allMoves = getAllMovesOnBoard(rootBoard, color);
    if (allMoves.empty()) {
        delete root;
        return Move(-1, -1, -1, -1, -1, -1);
    }
    if (allMoves.size() == 1) {
        delete root;
        return allMoves[0];
    }

    // 进行指定次数的MCTS迭代
    for (int i = 0; i < iterations; i++) {
        // 1. 选择阶段：从根节点选择到叶节点
        MCTSNode* selectedNode = selectNode(root,c);

        // 2. 扩展阶段：扩展叶节点
        MCTSNode* expandedNode;
        if (selectedNode->isFullyExpanded()) {
            expandedNode = selectedNode;
        } else {
            expandedNode = expandNode(selectedNode);
        }

        // 3. 模拟阶段：从扩展节点进行随机模拟
        double simulationResult = simulateGame(expandedNode, maxDepth);

        // 4. 回溯阶段：将结果回传到根节点
        backpropagate(expandedNode, simulationResult);
    }

    // 选择访问次数最多的走法作为最终决策
    MCTSNode* bestChild = nullptr;
    int maxVisits = -1;

    for (MCTSNode* child : root->children) {
        if (child->visits > maxVisits) {
            maxVisits = child->visits;
            bestChild = child;
        }
    }

    Move bestMove = bestChild->move;

    // 清理内存
    delete root;

    return bestMove;
}
// < -------- MC --------- >//

// ===================== 主函数 =====================
int main(){
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

	// 使用MCTS选择最优走法
	srand(time(0));
	Move bestMove = monteCarloSearch(currBotColor, iterations, c); // 使用全局超参数
	
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
