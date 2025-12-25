#define UNICODE
#define _UNICODE

#include <SFML/Graphics.hpp>
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
#include <chrono>

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

int searchDepth = 3; // Alpha-Beta搜索深度
const int K1_QUEEN_MOVES = 80;   // 第一阶段保留的皇后移动位置数
const int K2_FULL_MOVES = 50;    // 第二阶段保留的完整走法数
const int INF = 1000000000; // 无穷大值

// 性能分析统计
struct PerfStats {
    double totalGenMoveTime = 0;  // move候选两阶段总耗时（毫秒）
    double totalTime = 0;  // 这次操作总耗时（毫秒）
    
    void reset() {
        totalGenMoveTime = 0;
        totalTime = 0;
    }
    
    void printStats() const {
        cout << "\n========== 时间统计 ==========\n";
        cout << "Move候选生成时间: " << totalGenMoveTime << " ms\n";
        cout << "总耗时: " << totalTime << " ms\n";
        cout << "================================\n";
    }
} perfStats;

// 权重参数
const int WEIGHT_CONNECTIVITY = 10;  // 联通点权重
const int WEIGHT_MOBILITY = 3;       // 灵活性权重
const int WEIGHT_CENTER = 1;         // 中心性权重


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

// GUI 相关变量
int currentTurn = grid_black;
bool humanIsBlack = true;
int selX = -1, selY = -1;
int moveX = -1, moveY = -1;
vector<Path> validMovePaths;  // 可移动的路径
vector<Path> validObstaclePaths;  // 可放置障碍的路径
GameStatus gameStatus = GAME_ONGOING;  // 游戏状态
bool needAIMove = false;  // 是否需要AI走棋（用于延迟AI走棋，先显示人类走棋结果）




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
                    if (!inMap(xx, yy) || board[xx][yy] != 0)
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

// < --------- score ------------>

// 计算单个位置的灵活性（周围8个方向可移动的格子数）
int countPositionMobility(int board[GRIDSIZE][GRIDSIZE], int x, int y) {
    int mobility = 0;
    for(int k = 0; k < 8; k++) {
        for(int delta = 1; delta < GRIDSIZE; delta++) {
            int nx = x + dx[k] * delta;
            int ny = y + dy[k] * delta;
            if(!inMap(nx, ny) || board[nx][ny] != 0) break;
            mobility++;
        }
    }
    return mobility;
}

// 计算位置的中心性得分（越靠近中心越高）
int centerScore(int x, int y) {
    int center = GRIDSIZE / 2;
    // 到中心的曼哈顿距离，最大为7+7=14，转换为得分（越小越好）
    int dist = abs(x - center) + abs(y - center);
    return max(0, 8 - dist);  // 中心为8分，边角为0分
}

// 计算某方所有Queen的联通区域大小（共享联通区域的Queen只计算一次）
// 同时累加灵活性和中心性
struct ColorScore {
    int connectivity;  // 联通点总数
    int mobility;      // 灵活性总和
    int centerBonus;   // 中心性加成
};

ColorScore calculateColorScore(int board[GRIDSIZE][GRIDSIZE], int color) {
    ColorScore result = {0, 0, 0};
    
    // 找到所有该颜色的Queen位置
    vector<pair<int,int>> queens;
    for(int i = 0; i < GRIDSIZE; i++) {
        for(int j = 0; j < GRIDSIZE; j++) {
            if(board[i][j] == color) {
                queens.push_back({i, j});
                // 累加Queen位置的中心性
                result.centerBonus += centerScore(i, j);
            }
        }
    }
    
    if(queens.empty()) return result;
    
    // 使用并查集思想：BFS从所有Queen同时出发，计算联通区域
    // 已访问标记
    bool visited[GRIDSIZE][GRIDSIZE] = {false};
    queue<pair<int,int>> q;
    
    // 所有Queen作为起点
    for(auto& qpos : queens) {
        visited[qpos.first][qpos.second] = true;
        q.push(qpos);
    }
    
    while(!q.empty()) {
        auto [cx, cy] = q.front();
        q.pop();
        
        // 计算该点的灵活性（可移动方向数）
        int localMobility = 0;
        
        for(int k = 0; k < 8; k++) {
            // 皇后走法：沿着方向一直走
            for(int delta = 1; delta < GRIDSIZE; delta++) {
                int nx = cx + dx[k] * delta;
                int ny = cy + dy[k] * delta;
                
                if(!inMap(nx, ny)) break;
                
                // 遇到障碍物或对方棋子，停止
                if(board[nx][ny] == OBSTACLE || board[nx][ny] == -color) break;
                
                // 遇到己方Queen，共享联通区域（不重复计算）
                if(board[nx][ny] == color) break;
                
                // 空格，计入联通区域
                if(!visited[nx][ny]) {
                    visited[nx][ny] = true;
                    result.connectivity++;
                    q.push({nx, ny});
                }
                localMobility++;
            }
        }
        
        result.mobility += localMobility;
    }
    
    return result;
}

// 新的局面评估函数
int evaluateBoard(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    // 快速检查终局：是否有一方无法移动
    bool myCanMove = false, oppCanMove = false;
    
    for(int i = 0; i < GRIDSIZE && (!myCanMove || !oppCanMove); i++) {
        for(int j = 0; j < GRIDSIZE && (!myCanMove || !oppCanMove); j++) {
            if(board[i][j] == myColor || board[i][j] == -myColor) {
                int color = board[i][j];
                for(int k = 0; k < 8 && !(color == myColor ? myCanMove : oppCanMove); k++) {
                    int nx = i + dx[k];
                    int ny = j + dy[k];
                    if(inMap(nx, ny) && board[nx][ny] == 0) {
                        if(color == myColor) myCanMove = true;
                        else oppCanMove = true;
                    }
                }
            }
        }
    }
    
    if(!myCanMove && !oppCanMove) {
        // 双方都无法移动，比较联通区域
        ColorScore myScore = calculateColorScore(board, myColor);
        ColorScore oppScore = calculateColorScore(board, -myColor);
        if(myScore.connectivity > oppScore.connectivity) return INF - 1;
        else if(myScore.connectivity < oppScore.connectivity) return -INF + 1;
        else return 0;
    }
    
    if(!myCanMove) return -INF + 1;  // 我方输了
    if(!oppCanMove) return INF - 1;  // 对方输了
    
    // 计算双方得分
    ColorScore myScore = calculateColorScore(board, myColor);
    ColorScore oppScore = calculateColorScore(board, -myColor);
    
    // 综合评分：联通点差 * 权重 + 灵活性差 * 权重 + 中心性差 * 权重
    int totalScore = (myScore.connectivity - oppScore.connectivity) * WEIGHT_CONNECTIVITY
                   + (myScore.mobility - oppScore.mobility) * WEIGHT_MOBILITY
                   + (myScore.centerBonus - oppScore.centerBonus) * WEIGHT_CENTER;
    
    return totalScore;
}

// 计算黑方得分（从黑方视角）
double getBlackScore(int board[GRIDSIZE][GRIDSIZE]) {
    return evaluateBoard(board, grid_black);
}

// 计算白方得分（从白方视角）
double getWhiteScore(int board[GRIDSIZE][GRIDSIZE]) {
    return evaluateBoard(board, grid_white);
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

// < -------- Alpha-Beta剪枝 --------- >

// 置换表结构
struct TTEntry {
    long long hash;
    int depth;
    int score;
    int flag; // 0: EXACT, 1: LOWERBOUND, 2: UPPERBOUND
    Move bestMove;
    
    TTEntry() : hash(0), depth(-1), score(0), flag(0), bestMove() {}
};

const int TT_SIZE = 1048576; // 2^20 entries
TTEntry transpositionTable[TT_SIZE];

// Zobrist哈希值
long long zobristTable[GRIDSIZE][GRIDSIZE][3]; // [x][y][piece: 0=black, 1=white, 2=obstacle]
long long zobristBlackTurn;

// 初始化Zobrist哈希表
void initZobrist() {
    srand(20251225); // 固定种子以保证一致性
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            for(int k=0; k<3; k++) {
                zobristTable[i][j][k] = ((long long)rand() << 32) | rand();
            }
        }
    }
    zobristBlackTurn = ((long long)rand() << 32) | rand();
}

// 计算棋盘的Zobrist哈希值
long long computeHash(int board[GRIDSIZE][GRIDSIZE], int currentPlayer) {
    long long hash = 0;
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == grid_black) {
                hash ^= zobristTable[i][j][0];
            } else if(board[i][j] == grid_white) {
                hash ^= zobristTable[i][j][1];
            } else if(board[i][j] == OBSTACLE) {
                hash ^= zobristTable[i][j][2];
            }
        }
    }
    if(currentPlayer == grid_black) hash ^= zobristBlackTurn;
    return hash;
}

// ========== 两阶段走法生成 ==========

// 皇后移动位置结构（第一阶段）
struct QueenMove {
    int x0, y0;  // 起始位置
    int x1, y1;  // 目标位置
    int score;   // 评分（灵活性 + 中心性）
    
    QueenMove() : x0(-1), y0(-1), x1(-1), y1(-1), score(0) {}
    QueenMove(int x0_, int y0_, int x1_, int y1_, int s = 0) 
        : x0(x0_), y0(y0_), x1(x1_), y1(y1_), score(s) {}
};

// 第一阶段：评估皇后移动位置的得分
// 评分标准：灵活性（周围可移动格子数）+ 中心性 + 3x3空白点（远离障碍和边界）
int scoreQueenMove(int board[GRIDSIZE][GRIDSIZE], int x0, int y0, int x1, int y1, int color) {
    // 临时移动皇后
    int tempBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(tempBoard, board);
    tempBoard[x0][y0] = 0;
    tempBoard[x1][y1] = color;
    
    // 计算灵活性（从新位置可到达的格子数）
    int mobility = countPositionMobility(tempBoard, x1, y1);
    
    // 计算中心性
    int center = centerScore(x1, y1);
    
    // 统计3x3范围内的空白点个数（越多说明离障碍和边界越远）
    int emptyCount = 0;
    for(int di = -1; di <= 1; di++) {
        for(int dj = -1; dj <= 1; dj++) {
            int ni = x1 + di;
            int nj = y1 + dj;
            if(inMap(ni, nj) && tempBoard[ni][nj] == 0) {
                emptyCount++;
            }
        }
    }
    
    return mobility * 3 + center + emptyCount * 2;  // 灵活性权重最高，空白点次之
}

// 第一阶段：生成并筛选前K1个最佳皇后移动位置
vector<QueenMove> generateTopQueenMoves(int board[GRIDSIZE][GRIDSIZE], int color, int k1) {
    vector<QueenMove> allQueenMoves;
    
    // 遍历所有该颜色的皇后
    for(int i = 0; i < GRIDSIZE; i++) {
        for(int j = 0; j < GRIDSIZE; j++) {
            if(board[i][j] != color) continue;
            
            // 生成该皇后的所有可能移动位置
            for(int k = 0; k < 8; k++) {
                for(int delta = 1; delta < GRIDSIZE; delta++) {
                    int x1 = i + dx[k] * delta;
                    int y1 = j + dy[k] * delta;
                    
                    if(!inMap(x1, y1) || board[x1][y1] != 0) break;
                    
                    // 评估这个移动位置
                    int score = scoreQueenMove(board, i, j, x1, y1, color);
                    allQueenMoves.push_back(QueenMove(i, j, x1, y1, score));
                }
            }
        }
    }
    
    // 按得分降序排序
    sort(allQueenMoves.begin(), allQueenMoves.end(), 
         [](const QueenMove& a, const QueenMove& b) { return a.score > b.score; });
    
    // 只保留前K1个
    if((int)allQueenMoves.size() > k1) {
        allQueenMoves.resize(k1);
    }
    
    return allQueenMoves;
}

// 第二阶段：对走法评分（完整走法 = 皇后移动 + 射箭）
// 使用快速启发式评分：隔断性 + 阻止对手灵活性
int scoreFullMove(int board[GRIDSIZE][GRIDSIZE], const Move& move, int color) {
    // 不需要执行完整走法，只需评估障碍位置
    int x2 = move.x2, y2 = move.y2;
    
    // 1. 隔断性：统计障碍位置周围8个方向的障碍物/边界个数（越多越好，说明容易形成闭合）
    int blockCount = 0;
    for(int k = 0; k < 8; k++) {
        int nx = x2 + dx[k];
        int ny = y2 + dy[k];
        // 边界或已有障碍物/棋子
        if(!inMap(nx, ny) || board[nx][ny] == OBSTACLE) {
            blockCount++;
        }
    }
    
    // 2. 阻止对手灵活性：统计障碍位置周围一圈的对手棋子个数（越多越好，说明阻挡了对手的移动）
    int oppCount = 0;
    int oppColor = -color;
    for(int k = 0; k < 8; k++) {
        int nx = x2 + dx[k];
        int ny = y2 + dy[k];
        if(inMap(nx, ny) && board[nx][ny] == oppColor) {
            oppCount++;
        }
    }
    
    // 综合评分：隔断性权重更高
    return blockCount * 10 + oppCount * 5;
}

// 第二阶段：从选定的皇后移动生成完整走法，并筛选前K2个
vector<Move> generateTopFullMoves(int board[GRIDSIZE][GRIDSIZE], 
                                   const vector<QueenMove>& queenMoves, 
                                   int color, int k2) {
    
    vector<pair<int, Move>> scoredMoves;  // (score, move)
    
    for(const auto& qm : queenMoves) {
        // 临时移动皇后
        int tempBoard[GRIDSIZE][GRIDSIZE];
        copyBoard(tempBoard, board);
        tempBoard[qm.x0][qm.y0] = 0;
        tempBoard[qm.x1][qm.y1] = color;
        // 生成所有可能的射箭位置
        for(int k = 0; k < 8; k++) {
            for(int delta = 1; delta < GRIDSIZE; delta++) {
                int x2 = qm.x1 + dx[k] * delta;
                int y2 = qm.y1 + dy[k] * delta;
                
                if(!inMap(x2, y2)) break;
                
                // 可以射到原位置（已经空了）或空格
                if(tempBoard[x2][y2] != 0 && !(x2 == qm.x0 && y2 == qm.y0)) break;
                
                // 创建完整走法
                Move fullMove(qm.x0, qm.y0, qm.x1, qm.y1, x2, y2);
                int score = scoreFullMove(board, fullMove, color);
                
                scoredMoves.push_back({score, fullMove});
            }
        }
        
    }
    // 按得分降序排序
    sort(scoredMoves.begin(), scoredMoves.end(),
         [](const pair<int,Move>& a, const pair<int,Move>& b) { return a.first > b.first; });
    
    // 只保留前K2个
    vector<Move> result;
    int count = min(k2, (int)scoredMoves.size());
    for(int i = 0; i < count; i++) {
        result.push_back(scoredMoves[i].second);
    }
    
    return result;
}

// 两阶段走法生成主函数
vector<Move> generateMovesWithTwoPhase(int board[GRIDSIZE][GRIDSIZE], int color) {
    auto phaseStartTime = chrono::high_resolution_clock::now();
    
    // 第一阶段：生成并筛选皇后移动位置
    vector<QueenMove> topQueenMoves = generateTopQueenMoves(board, color, K1_QUEEN_MOVES);
    if(topQueenMoves.empty()) {
        return vector<Move>();  // 无法移动
    }
    
    // 第二阶段：生成并筛选完整走法
    vector<Move> topFullMoves = generateTopFullMoves(board, topQueenMoves, color, K2_FULL_MOVES);

    auto phaseEndTime = chrono::high_resolution_clock::now();
    auto phaseDuration = chrono::duration_cast<chrono::microseconds>(phaseEndTime - phaseStartTime).count();
    perfStats.totalGenMoveTime += phaseDuration / 1000.0;  // 转换为毫秒
    
    return topFullMoves;
}

// ========== 两阶段走法生成结束 ==========

// Alpha-Beta剪枝搜索（带置换表 + 两阶段走法生成）
int alphaBetaSearch(int board[GRIDSIZE][GRIDSIZE], int depth, int alpha, int beta, 
                    int color, Move& bestMove, long long hash) {
    // 查询置换表
    int ttIndex = (hash & 0x7FFFFFFFFFFFFFFF) % TT_SIZE;
    TTEntry* ttEntry = &transpositionTable[ttIndex];
    
    if(ttEntry->hash == hash && ttEntry->depth >= depth) {
        if(ttEntry->flag == 0) { // EXACT
            bestMove = ttEntry->bestMove;
            return ttEntry->score;
        } else if(ttEntry->flag == 1) { // LOWERBOUND
            alpha = max(alpha, ttEntry->score);
        } else if(ttEntry->flag == 2) { // UPPERBOUND
            beta = min(beta, ttEntry->score);
        }
        if(alpha >= beta) {
            bestMove = ttEntry->bestMove;
            return ttEntry->score;
        }
    }
    
    const int alphaOrig = alpha;
    const int betaOrig = beta;

    // 使用两阶段走法生成
    vector<Move> moves = generateMovesWithTwoPhase(board, color);
    
    // 终止条件
    if(depth == 0 || moves.empty()) {
        int score = evaluateBoard(board, color);
        return score;
    }
    
    int bestScore = -INF;
    Move localBestMove;
    
    for(const Move& move : moves) {
        // 执行走法
        int newBoard[GRIDSIZE][GRIDSIZE];
        copyBoard(newBoard, board);
        if(!ProcStepOnBoard(newBoard, move.x0, move.y0, move.x1, move.y1, move.x2, move.y2, color)) {
            continue; // 防御性检查
        }
        
        // 计算新哈希值
        long long newHash = hash;
        // 移除旧位置的棋子
        if(board[move.x0][move.y0] == grid_black) newHash ^= zobristTable[move.x0][move.y0][0];
        else if(board[move.x0][move.y0] == grid_white) newHash ^= zobristTable[move.x0][move.y0][1];
        // 添加新位置的棋子
        if(color == grid_black) newHash ^= zobristTable[move.x1][move.y1][0];
        else newHash ^= zobristTable[move.x1][move.y1][1];
        // 添加障碍物
        newHash ^= zobristTable[move.x2][move.y2][2];
        // 切换玩家
        newHash ^= zobristBlackTurn;
        
        Move oppBestMove;
        int score = -alphaBetaSearch(newBoard, depth - 1, -beta, -alpha, -color, oppBestMove, newHash);
        
        if(score > bestScore) {
            bestScore = score;
            localBestMove = move;
        }
        
        alpha = max(alpha, score);
        if(alpha >= beta) {
            break; // Beta剪枝
        }
    }
    
    // 存储到置换表
    ttEntry->hash = hash;
    ttEntry->depth = depth;
    ttEntry->score = bestScore;
    ttEntry->bestMove = localBestMove;
    
    if(bestScore <= alphaOrig) {
        ttEntry->flag = 2; // UPPERBOUND
    } else if(bestScore >= betaOrig) {
        ttEntry->flag = 1; // LOWERBOUND
    } else {
        ttEntry->flag = 0; // EXACT
    }
    
    bestMove = localBestMove;
    
    return bestScore;
}

// 迭代加深搜索
Move iterativeDeepeningSearch(int color, int maxDepth) {
    auto totalStartTime = chrono::high_resolution_clock::now();
    
    // 重置性能统计
    perfStats.reset();
    
    // initZobrist();
    Move bestMove;
    long long hash = computeHash(gridInfo, color);
    
    // 逐步增加深度
    for(int depth = 1; depth <= maxDepth; depth++) {
        Move currentBestMove;
        int score = alphaBetaSearch(gridInfo, depth, -INF, INF, color, currentBestMove, hash);
        
        if(currentBestMove.x0 >= 0) {
            bestMove = currentBestMove;
        }
        
        // 输出调试信息
        cout << "[Alpha-Beta] Depth " << depth << ": score = " << score << endl;
    }
    
    auto totalEndTime = chrono::high_resolution_clock::now();
    auto totalDuration = chrono::duration_cast<chrono::microseconds>(totalEndTime - totalStartTime).count();
    perfStats.totalTime = totalDuration / 1000.0;  // 转换为毫秒
    
    // 输出详细的性能统计
    perfStats.printStats();
    
    return bestMove;
}
// < -------- Alpha-Beta剪枝结束 --------- >

// ===================== 保存 / 读取 =====================
void saveGame(){
    ofstream out("save.txt");
    out<<currentTurn<<"\n";
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++) out<<gridInfo[i][j]<<" ";
        out<<"\n";
    }
    out.close();
    cout<<"[Saved]\n";
}

void loadGame(){
    ifstream in("save.txt");
    if(!in.is_open()) {
        cout<<"[Load failed: file not found]\n";
        return;
    }
    in>>currentTurn;
    for(int i=0;i<8;i++)
        for(int j=0;j<8;j++) in>>gridInfo[i][j];
    in.close();
    selX=selY=moveX=moveY=-1;
    validMovePaths.clear();
    validObstaclePaths.clear();
    gameStatus = checkGameStatus(gridInfo);  // 检查游戏状态
    cout<<"[Loaded]\n";
}

// ===================== AI 移动 =====================
void aiMove(){
    if(gameStatus != GAME_ONGOING) return;  // 游戏已结束，不再移动
    // 清除用户的选择状态（防止在AI走棋时用户有未完成的选择）
    selX=selY=moveX=moveY=-1;
    validMovePaths.clear();
    validObstaclePaths.clear();

    // 根据人类颜色设置AI颜色
    int aiColor = humanIsBlack ? grid_white : grid_black;

    // 使用 Alpha-Beta 剪枝算法进行决策
    Move bestMove = iterativeDeepeningSearch(aiColor, searchDepth);
    cout<<"AI selected move: ("<<bestMove.x0<<","<<bestMove.y0<<") -> ("
        <<bestMove.x1<<","<<bestMove.y1<<") with obstacle at ("
        <<bestMove.x2<<","<<bestMove.y2<<")\n";
    if(bestMove.x0 >= 0) {
        ProcStepOnBoard(gridInfo, bestMove.x0, bestMove.y0, bestMove.x1, bestMove.y1,
                       bestMove.x2, bestMove.y2, currentTurn);
        currentTurn = -currentTurn;
        gameStatus = checkGameStatus(gridInfo);  // 检查游戏状态
    }
    cout<<"AI move done\n";
}

// ===================== 主函数 =====================
int main(){
    // 初始化Zobrist哈希表
    initZobrist();

    srand(time(0));

	// 初始化棋盘
    gridInfo[0][2]=gridInfo[2][0]=gridInfo[5][0]=gridInfo[7][2]=grid_black;
    gridInfo[0][5]=gridInfo[2][7]=gridInfo[5][7]=gridInfo[7][5]=grid_white;
    gameStatus = GAME_ONGOING;
    needAIMove = false;

    sf::RenderWindow win(sf::VideoMode(640,820),"Amazon Chess");
    sf::Font font;
    bool fontLoaded = false;

    // 尝试多个可能的字体路径
    const char* fontPaths[] = {
        "arial.ttf",                                    // Current directory
        "C:/Windows/Fonts/arial.ttf",                   // Windows system fonts directory (lowercase)
        "C:/Windows/Fonts/Arial.ttf",                  // Windows system fonts directory (uppercase)
        "C:/Windows/Fonts/msyh.ttc",                   // Microsoft YaHei (fallback if Arial unavailable)
        "C:/Windows/Fonts/simsun.ttc",                 // SimSun (fallback if Arial unavailable)
        "C:/Windows/Fonts/calibri.ttf"                 // Calibri (fallback if Arial unavailable)
    };

    for(int i = 0; i < sizeof(fontPaths)/sizeof(fontPaths[0]); i++) {
        fontLoaded = font.loadFromFile(fontPaths[i]);
        if(fontLoaded) {
            cout << "[Font loaded: " << fontPaths[i] << "]\n";
            break;
        }
    }

    if(!fontLoaded) {
        cout << "[Warning: Could not load any font file, text may not display]\n";
        cout << "[Please place arial.ttf in the program directory or ensure system fonts are accessible]\n";
    }

    // 创建说明文本和得分显示
    vector<sf::Text> helpTexts;
    sf::Text blackScoreText, whiteScoreText, gameStatusText;
    if(fontLoaded) {
        vector<string> helpStrings = {
            "Controls:",
            "S - Save game",
            "L - Load game",
            "A - Manual AI move",
            "C - Switch human color",
            "R - Reset game"
        };

        for(size_t i = 0; i < helpStrings.size(); i++) {
            sf::Text text;
            text.setFont(font);
            text.setString(helpStrings[i]);
            text.setCharacterSize(14);
            text.setFillColor(sf::Color::Black);
            text.setPosition(10, 640 + i * 18);
            helpTexts.push_back(text);
        }

        // 得分显示
        blackScoreText.setFont(font);
        blackScoreText.setCharacterSize(16);
        blackScoreText.setFillColor(sf::Color::Black);
        blackScoreText.setPosition(350, 640);

        whiteScoreText.setFont(font);
        whiteScoreText.setCharacterSize(16);
        whiteScoreText.setFillColor(sf::Color::Black);
        whiteScoreText.setPosition(350, 660);

        gameStatusText.setFont(font);
        gameStatusText.setCharacterSize(18);
        gameStatusText.setFillColor(sf::Color::Red);
        gameStatusText.setPosition(350, 680);
    }

    while(win.isOpen()){
        sf::Event e;
        while(win.pollEvent(e)){
            if(e.type==sf::Event::Closed) win.close();

            if(e.type==sf::Event::KeyPressed){
                if(e.key.code==sf::Keyboard::S) saveGame();
                if(e.key.code==sf::Keyboard::L) loadGame();
                if(e.key.code==sf::Keyboard::A) {
                    // 清除用户的选择状态
                    selX=selY=moveX=moveY=-1;
                    validMovePaths.clear();
                    validObstaclePaths.clear();
                    // 检查是否是AI的回合
                    bool isAITurn = (humanIsBlack && currentTurn == grid_white) ||
                                   (!humanIsBlack && currentTurn == grid_black);
                    if(isAITurn && gameStatus == GAME_ONGOING) {
                        aiMove();

                    }
                }
                if(e.key.code==sf::Keyboard::C) {
                    humanIsBlack=!humanIsBlack;
                    selX=selY=moveX=moveY=-1;
                    validMovePaths.clear();
                    validObstaclePaths.clear();
                }
                if(e.key.code==sf::Keyboard::R){
                    memset(gridInfo,0,sizeof(gridInfo));
                    gridInfo[0][2]=gridInfo[2][0]=gridInfo[5][0]=gridInfo[7][2]=grid_black;
                    gridInfo[0][5]=gridInfo[2][7]=gridInfo[5][7]=gridInfo[7][5]=grid_white;
                    currentTurn = grid_black;
                    gameStatus = GAME_ONGOING;
                    needAIMove = false;
                    selX=selY=moveX=moveY=-1;
                    validMovePaths.clear();
                    validObstaclePaths.clear();
                }
            }

            if(e.type==sf::Event::MouseButtonPressed){
                int x=e.mouseButton.x/CELL;
                int y=e.mouseButton.y/CELL;
                if(!inMap(x,y)) continue;

                // 检查是否是人类的回合
                bool isHumanTurn = (humanIsBlack && currentTurn == grid_black) ||
                                  (!humanIsBlack && currentTurn == grid_white);

                if(!isHumanTurn) continue;
                if(gameStatus != GAME_ONGOING) continue;  // 游戏已结束，不允许操作

                // 第一步：选择棋子
                if(selX==-1 && gridInfo[x][y]==currentTurn){
                    selX=x; selY=y;
                    moveX=moveY=-1;
                    validMovePaths = getValidMovePaths(gridInfo, selX, selY);
                    validObstaclePaths.clear();
                }
                // 第二步：点击路径上的点确认移动目标
                else if(selX!=-1 && moveX==-1){
                    // 检查是否点击了某个可移动路径上的点
                    Path* clickedPath = nullptr;
                    for(auto& path : validMovePaths) {
                        if(path.contains(x, y)) {
                            clickedPath = &path;
                            break;
                        }
                    }
                    if(clickedPath && clickedPath->points.size() > 1) {
                        // 找到路径的终点（最后一个点）
                        auto& lastPoint = clickedPath->points.back();
                        moveX = lastPoint.first;
                        moveY = lastPoint.second;
                        // 创建临时棋盘，模拟移动后的状态（原位置变空，新位置有棋子）
                        int tempBoard[GRIDSIZE][GRIDSIZE];
                        copyBoard(tempBoard, gridInfo);
                        tempBoard[selX][selY] = 0;  // 原位置变空
                        tempBoard[moveX][moveY] = currentTurn;  // 新位置有棋子
                        // 计算可以放置障碍的路径（从移动后的位置发射，也可以放回原位置）
                        validObstaclePaths = getValidObstaclePaths(tempBoard, moveX, moveY, selX, selY);
                    } else if(gridInfo[x][y]==currentTurn) {
                        // 如果点击了另一个自己的棋子，重新选择
                        selX=x; selY=y;
                        moveX=moveY=-1;
                        validMovePaths = getValidMovePaths(gridInfo, selX, selY);
                        validObstaclePaths.clear();
                    } else {
                        // 点击了无效位置，取消选择
                        selX=selY=moveX=moveY=-1;
                        validMovePaths.clear();
                        validObstaclePaths.clear();
                    }
                }
                // 第三步：点击障碍路径上的点完成走棋
                else if(selX!=-1 && moveX!=-1){
                    // 检查是否点击了某个可放置障碍路径上的点
                    Path* clickedPath = nullptr;
                    for(auto& path : validObstaclePaths) {
                        if(path.contains(x, y)) {
                            clickedPath = &path;
                            break;
                        }
                    }
                    if(clickedPath && clickedPath->points.size() > 1) {
                        // 找到路径的终点（最后一个点）
                        auto& lastPoint = clickedPath->points.back();
                        int obstacleX = lastPoint.first;
                        int obstacleY = lastPoint.second;
                        // 验证完整走法是否合法
                        if(ProcStepOnBoard(gridInfo, selX, selY, moveX, moveY, obstacleX, obstacleY, currentTurn)) {
                            selX=selY=moveX=moveY=-1;
                            validMovePaths.clear();
                            validObstaclePaths.clear();
                            currentTurn=-currentTurn;
                            gameStatus = checkGameStatus(gridInfo);  // 检查游戏状态
                            // 人类走完后，如果是AI的回合，设置标志延迟AI走棋（先显示人类走棋的结果）
                            bool isAITurn = (humanIsBlack && currentTurn == grid_white) ||
                                           (!humanIsBlack && currentTurn == grid_black);
                            if(isAITurn && gameStatus == GAME_ONGOING) {
                                needAIMove = true;  // 设置标志，在主循环中延迟执行
                            }
                        } else {
                            // 如果走法不合法，重置选择
                            selX=selY=moveX=moveY=-1;
                            validMovePaths.clear();
                            validObstaclePaths.clear();
                        }
                    } else {
                        // 点击了无效位置，取消选择
                        selX=selY=moveX=moveY=-1;
                        validMovePaths.clear();
                        validObstaclePaths.clear();
                    }
                }
            }
        }

        // 如果需要AI走棋，先绘制当前状态（显示人类走棋的结果），然后再让AI走棋
        if(needAIMove) {
            // 先绘制一次，显示人类走棋的结果（障碍）
            win.clear(sf::Color::White);
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    sf::RectangleShape cell(sf::Vector2f(CELL-1,CELL-1));
                    cell.setPosition(i*CELL,j*CELL);
                    // 国际象棋棋盘配色：浅色格子（米色）和深色格子（深棕色）
                cell.setFillColor((i+j)%2?sf::Color(181,136,99):sf::Color(240,217,181));
                    win.draw(cell);

                    if(gridInfo[i][j]!=0){
                        sf::CircleShape c(CELL/2-5);
                        c.setPosition(i*CELL+5,j*CELL+5);
                    if(gridInfo[i][j]==grid_black) c.setFillColor(sf::Color::Black);
                    else if(gridInfo[i][j]==grid_white) {
                        // 白色棋子使用深灰色边框，内部浅灰色，以便在浅色格子上也能看清
                        c.setFillColor(sf::Color(200,200,200));
                        c.setOutlineThickness(2);
                        c.setOutlineColor(sf::Color::Black);
                    }
                    else c.setFillColor(sf::Color::Red);
                    win.draw(c);
                    }
                }
            }

            // 绘制说明文本和得分
            if(fontLoaded) {
                for(auto& text : helpTexts) {
                    win.draw(text);
                }

                double blackScore = getBlackScore(gridInfo);
                double whiteScore = getWhiteScore(gridInfo);
                // 格式化得分，保留1位小数
                char scoreBuf1[50], scoreBuf2[50];
                sprintf(scoreBuf1, "Black Score: %.1f", blackScore);
                sprintf(scoreBuf2, "White Score: %.1f", whiteScore);
                blackScoreText.setString(scoreBuf1);
                whiteScoreText.setString(scoreBuf2);
                win.draw(blackScoreText);
                win.draw(whiteScoreText);

                string statusStr = "Game Ongoing - " + string(currentTurn == grid_black ? "Black's Turn" : "White's Turn");
                gameStatusText.setString(statusStr);
                gameStatusText.setFillColor(sf::Color::Black);
                win.draw(gameStatusText);
            }

            win.display();  // 先显示人类走棋的结果

            // 然后让AI走棋
            needAIMove = false;
            aiMove();
            // 注意：这里不使用continue，让AI移动的结果在正常的绘制循环中显示
        }

        // 绘制
        win.clear(sf::Color::White);
        for(int i=0;i<8;i++){
            for(int j=0;j<8;j++){
                sf::RectangleShape cell(sf::Vector2f(CELL-1,CELL-1));
                cell.setPosition(i*CELL,j*CELL);
                // 国际象棋棋盘配色：浅色格子（米色）和深色格子（深棕色）
                cell.setFillColor((i+j)%2?sf::Color(181,136,99):sf::Color(240,217,181));
                win.draw(cell);

                if(gridInfo[i][j]!=0){
                    sf::CircleShape c(CELL/2-5);
                    c.setPosition(i*CELL+5,j*CELL+5);
                    if(gridInfo[i][j]==grid_black) c.setFillColor(sf::Color::Black);
                    else if(gridInfo[i][j]==grid_white) {
                        // 白色棋子使用深灰色边框，内部浅灰色，以便在浅色格子上也能看清
                        c.setFillColor(sf::Color(200,200,200));
                        c.setOutlineThickness(2);
                        c.setOutlineColor(sf::Color::Black);
                    }
                    else c.setFillColor(sf::Color::Red);
                    win.draw(c);
                }
            }
        }

        // 绘制选中的棋子（绿色高亮）
        if(selX != -1 && selY != -1) {
            sf::RectangleShape highlight(sf::Vector2f(CELL-2,CELL-2));
            highlight.setPosition(selX*CELL+1, selY*CELL+1);
            highlight.setFillColor(sf::Color(0,255,0,150));
            win.draw(highlight);
        }

        // 绘制可移动的路径（蓝色高亮，显示路径上的所有点）
        if(selX != -1 && selY != -1 && moveX == -1) {
            for(auto& path : validMovePaths) {
                for(auto& pos : path.points) {
                    // 跳过起点（已经用绿色高亮显示了）
                    if(pos.first == selX && pos.second == selY) continue;
                    sf::RectangleShape highlight(sf::Vector2f(CELL-2,CELL-2));
                    highlight.setPosition(pos.first*CELL+1, pos.second*CELL+1);
                    highlight.setFillColor(sf::Color(0,0,255,120));
                    win.draw(highlight);
                }
            }
        }

        // 绘制已选择的移动目标路径（深蓝色高亮）
        if(moveX != -1 && moveY != -1) {
            // 找到对应的路径并高亮显示
            for(auto& path : validMovePaths) {
                if(path.contains(moveX, moveY)) {
                    for(auto& pos : path.points) {
                        sf::RectangleShape highlight(sf::Vector2f(CELL-2,CELL-2));
                        highlight.setPosition(pos.first*CELL+1, pos.second*CELL+1);
                        highlight.setFillColor(sf::Color(0,0,200,150));
                        win.draw(highlight);
                    }
                    break;
                }
            }
        }

        // 绘制可放置障碍的路径（黄色高亮，显示路径上的所有点）
        if(moveX != -1 && moveY != -1) {
            for(auto& path : validObstaclePaths) {
                for(auto& pos : path.points) {
                    // 跳过起点（移动后的位置）
                    if(pos.first == moveX && pos.second == moveY) continue;
                    sf::RectangleShape highlight(sf::Vector2f(CELL-2,CELL-2));
                    highlight.setPosition(pos.first*CELL+1, pos.second*CELL+1);
                    highlight.setFillColor(sf::Color(255,255,0,100));
                    win.draw(highlight);
                }
            }
        }

        // 绘制说明文本（仅在字体加载成功时）
        if(fontLoaded) {
            for(auto& text : helpTexts) {
                win.draw(text);
            }

            // 更新并绘制得分
            double blackScore = getBlackScore(gridInfo);
            double whiteScore = getWhiteScore(gridInfo);
            // 格式化得分，保留1位小数
            char scoreBuf1[50], scoreBuf2[50];
            sprintf(scoreBuf1, "black score: %.1f", blackScore);
            sprintf(scoreBuf2, "white score: %.1f", whiteScore);
            blackScoreText.setString(scoreBuf1);
            whiteScoreText.setString(scoreBuf2);
            win.draw(blackScoreText);
            win.draw(whiteScoreText);

            // 更新并绘制游戏状态
            string statusStr = "";
            switch(gameStatus) {
                case GAME_ONGOING:
                    statusStr = "Game Ongoing - " + string(currentTurn == grid_black ? "Black's Turn" : "White's Turn");
                    gameStatusText.setFillColor(sf::Color::Black);
                    break;
                case BLACK_WIN:
                    statusStr = "Game Over - Black Wins!";
                    gameStatusText.setFillColor(sf::Color::Red);
                    break;
                case WHITE_WIN:
                    statusStr = "Game Over - White Wins!";
                    gameStatusText.setFillColor(sf::Color::Red);
                    break;
                case DRAW:
                    statusStr = "Game Over - Draw!";
                    gameStatusText.setFillColor(sf::Color::Blue);
                    break;
            }
            gameStatusText.setString(statusStr);
            win.draw(gameStatusText);
        }

        win.display();
    }
	return 0;
}
