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

int searchDepth = 2; // Alpha-Beta搜索深度
const int K1_QUEEN_MOVES = 80;   // 第一阶段保留的皇后移动位置数
const int K2_FULL_MOVES = 30;    // 第二阶段保留的完整走法数
const int INF = 1000000000; // 无穷大值

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
// 评分标准：灵活性（周围可移动格子数）+ 中心性
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
    
    return mobility * 3 + center;  // 灵活性权重更高
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
int scoreFullMove(int board[GRIDSIZE][GRIDSIZE], const Move& move, int color) {
    // 执行走法
    int tempBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(tempBoard, board);
    if(!ProcStepOnBoard(tempBoard, move.x0, move.y0, move.x1, move.y1, move.x2, move.y2, color)) {
        return -INF;  // 无效走法
    }
    // 计算局面差（走法得分 = 新局面评分）
    return evaluateBoard(tempBoard, color);
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
                
                if(score > -INF) {
                    scoredMoves.push_back({score, fullMove});
                }
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
    // 第一阶段：生成并筛选皇后移动位置
    vector<QueenMove> topQueenMoves = generateTopQueenMoves(board, color, K1_QUEEN_MOVES);
    if(topQueenMoves.empty()) {
        return vector<Move>();  // 无法移动
    }
    
    // 第二阶段：生成并筛选完整走法
    vector<Move> topFullMoves = generateTopFullMoves(board, topQueenMoves, color, K2_FULL_MOVES);
    
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
    Move bestMove;
    long long hash = computeHash(gridInfo, color);
    
    // 逐步增加深度
    for(int depth = 1; depth <= maxDepth; depth++) {
        Move currentBestMove;
        int score = alphaBetaSearch(gridInfo, depth, -INF, INF, color, currentBestMove, hash);
        
        if(currentBestMove.x0 >= 0) {
            bestMove = currentBestMove;
        }
    }
    
    return bestMove;
}
// < -------- Alpha-Beta剪枝结束 --------- >

// ===================== 主函数 =====================
int main(){
    // 初始化Zobrist哈希表
    initZobrist();
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
			ProcStepOnBoard(gridInfo,x0, y0, x1, y1, x2, y2, -currBotColor); // 模拟对方落子

		// 然后是自己当时的行动
		// 对手行动总比自己行动多一个
		if (i < turnID - 1)
		{
			cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2;
			if (x0 >= 0)
				ProcStepOnBoard(gridInfo, x0, y0, x1, y1, x2, y2, currBotColor); // 模拟己方落子
		}
	}

	// 做出决策（你只需修改以下部分）

	// 使用Alpha-Beta剪枝选择最优走法
	srand(time(0));
	Move bestMove = iterativeDeepeningSearch(currBotColor, searchDepth); // 使用迭代加深的Alpha-Beta搜索
	
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
