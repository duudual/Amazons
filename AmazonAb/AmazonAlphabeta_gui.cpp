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
int MAX_MOVES_PER_NODE = 20;
const int INF = 1000000000; // 无穷大值


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

//服务与随机走法
inline bool pathClear(
    int board[GRIDSIZE][GRIDSIZE],
    int x0, int y0,
    int x1, int y1
) {
    int dx = (x1 > x0) - (x1 < x0);
    int dy = (y1 > y0) - (y1 < y0);

    int x = x0 + dx;
    int y = y0 + dy;

    while (x != x1 || y != y1) {
        if (!inMap(x, y)) return false;
        if (board[x][y] != 0) return false;
        x += dx;
        y += dy;
    }

    // 终点本身必须为空
    return board[x1][y1] == 0;
}

inline bool pathClear(
    int board[GRIDSIZE][GRIDSIZE],
    int x0, int y0,
    int x1, int y1,
    int ex, int ey   // 允许经过/落在的例外格
) {
    int dx = (x1 > x0) - (x1 < x0);
    int dy = (y1 > y0) - (y1 < y0);

    int x = x0 + dx;
    int y = y0 + dy;

    while (x != x1 || y != y1) {
        if (!inMap(x, y)) return false;
        if (board[x][y] != 0 && !(x == ex && y == ey))
            return false;
        x += dx;
        y += dy;
    }

    // 终点检查
    if (board[x1][y1] != 0 && !(x1 == ex && y1 == ey))
        return false;

    return true;
}

// 用于simulate过程中的快速随机
bool getRandomMoveOnBoard(
    int board[GRIDSIZE][GRIDSIZE],
    int color,
    Move& outMove
) {
    vector<pair<int,int>> queens;
    for (int i = 0; i < GRIDSIZE; ++i)
        for (int j = 0; j < GRIDSIZE; ++j)
            if (board[i][j] == color)
                queens.emplace_back(i, j);

    if (queens.empty()) return false;

    for (int attempt = 0; attempt < 30; ++attempt) {
        auto [x0, y0] = queens[rand() % queens.size()];
        int dir1 = rand() % 8;

        int step1 = rand() % GRIDSIZE + 1;
        int x1 = x0 + dx[dir1] * step1;
        int y1 = y0 + dy[dir1] * step1;

        if (!inMap(x1, y1)) continue;
        if (!pathClear(board, x0, y0, x1, y1)) continue;

        int dir2 = rand() % 8;
        int step2 = rand() % GRIDSIZE + 1;
        int x2 = x1 + dx[dir2] * step2;
        int y2 = y1 + dy[dir2] * step2;

        if (!inMap(x2, y2)) continue;
        if (!pathClear(board, x1, y1, x2, y2, x0, y0)) continue;

        outMove = Move(x0, y0, x1, y1, x2, y2);
        return true;
    }
    return false;
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

// 得分 - 快速评估（已弃用，仅保留以供参考）
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
    //归一化
    score /= GRIDSIZE*GRIDSIZE;
    // 如果是白方(grid_white = -1)，正分好；如果是黑方(1)，负分好
    if (myColor == grid_white) return score; 
    else return -score;
}

// 计算棋子的移动能力（可达位置数量）
int countMobility(int board[GRIDSIZE][GRIDSIZE], int color) {
    int mobility = 0;
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == color) {
                for(int k=0; k<8; k++) {
                    for(int delta=1; delta<GRIDSIZE; delta++) {
                        int nx = i + dx[k] * delta;
                        int ny = j + dy[k] * delta;
                        if(!inMap(nx, ny) || board[nx][ny] != 0) break;
                        mobility++;
                    }
                }
            }
        }
    }
    return mobility;
}

// 改进的评估函数
int evaluateBoard(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    // 检查终局状态
    vector<Move> myMoves = getAllMovesOnBoard(board, myColor);
    vector<Move> oppMoves = getAllMovesOnBoard(board, -myColor);
    
    if(myMoves.empty() && oppMoves.empty()) {
        // 双方都无法移动，比较领地
        bfsTerritory(board, grid_black, blackDist);
        bfsTerritory(board, grid_white, whiteDist);
        int myTerritory = 0, oppTerritory = 0;
        for(int i=0; i<GRIDSIZE; i++) {
            for(int j=0; j<GRIDSIZE; j++) {
                if(board[i][j] == 0) {
                    int myDist = (myColor == grid_black) ? blackDist[i][j] : whiteDist[i][j];
                    int oppDist = (myColor == grid_black) ? whiteDist[i][j] : blackDist[i][j];
                    if(myDist < oppDist) myTerritory++;
                    else if(oppDist < myDist) oppTerritory++;
                }
            }
        }
        if(myTerritory > oppTerritory) return INF - 1;
        else if(myTerritory < oppTerritory) return -INF + 1;
        else return 0;
    }
    
    if(myMoves.empty()) return -INF + 1;  // 我方输了
    if(oppMoves.empty()) return INF - 1;  // 对方输了
    
    // 计算各项指标
    bfsTerritory(board, grid_black, blackDist);
    bfsTerritory(board, grid_white, whiteDist);
    
    int myDist[GRIDSIZE][GRIDSIZE], oppDist[GRIDSIZE][GRIDSIZE];
    if(myColor == grid_black) {
        memcpy(myDist, blackDist, sizeof(blackDist));
        memcpy(oppDist, whiteDist, sizeof(whiteDist));
    } else {
        memcpy(myDist, whiteDist, sizeof(whiteDist));
        memcpy(oppDist, blackDist, sizeof(blackDist));
    }
    
    // 1. 领地控制评分（最重要）
    int territoryScore = 0;
    int myCloseTerritory = 0;  // 我方近距离可控区域
    int oppCloseTerritory = 0; // 对方近距离可控区域
    
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == 0) {
                int md = myDist[i][j];
                int od = oppDist[i][j];
                
                if(md < od) {
                    territoryScore += (od - md);  // 差距越大，分数越高
                    if(md <= 2) myCloseTerritory++;
                } else if(od < md) {
                    territoryScore -= (md - od);
                    if(od <= 2) oppCloseTerritory++;
                }
            }
        }
    }
    
    // 2. 移动能力评分（行动力）
    int myMobility = countMobility(board, myColor);
    int oppMobility = countMobility(board, -myColor);
    int mobilityScore = (myMobility - oppMobility);
    
    // 3. 中心控制评分
    int centerScore = 0;
    int center = GRIDSIZE / 2;
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == 0) {
                int distToCenter = abs(i - center) + abs(i - center + 1) + 
                                   abs(j - center) + abs(j - center + 1);
                int md = myDist[i][j];
                int od = oppDist[i][j];
                if(md < od) {
                    centerScore += max(0, 8 - distToCenter);
                } else if(od < md) {
                    centerScore -= max(0, 8 - distToCenter);
                }
            }
        }
    }
    
    // 4. 棋子位置评分（更中心的位置更好）
    int positionScore = 0;
    for(int i=0; i<GRIDSIZE; i++) {
        for(int j=0; j<GRIDSIZE; j++) {
            if(board[i][j] == myColor) {
                int distToCenter = abs(i - center) + abs(i - center + 1) + 
                                   abs(j - center) + abs(j - center + 1);
                positionScore += (8 - distToCenter);
            } else if(board[i][j] == -myColor) {
                int distToCenter = abs(i - center) + abs(i - center + 1) + 
                                   abs(j - center) + abs(j - center + 1);
                positionScore -= (8 - distToCenter);
            }
        }
    }
    
    // 综合评分（权重调整）
    int totalScore = territoryScore * 100 + 
                     mobilityScore * 5 + 
                     centerScore * 3 +
                     positionScore * 2 +
                     (myCloseTerritory - oppCloseTerritory) * 20;
    
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
initZobrist(); // 全局初始化

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

// 对走法进行评分（用于移动排序）
int scoreMoveHeuristic(int board[GRIDSIZE][GRIDSIZE], const Move& move, int color) {
    int score = 0;
    
    // 临时执行走法
    int tempBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(tempBoard, board);
    ProcStepOnBoard(tempBoard, move.x0, move.y0, move.x1, move.y1, move.x2, move.y2, color);
    
    // 评估移动后的局面
    score = evaluateBoard(tempBoard, color);
    
    // 中心位置加分
    int center = GRIDSIZE / 2;
    int distToCenter = abs(move.x1 - center) + abs(move.x1 - center + 1) + 
                       abs(move.y1 - center) + abs(move.y1 - center + 1);
    score += (8 - distToCenter) * 10;
    
    return score;
}

// 移动排序
void orderMoves(int board[GRIDSIZE][GRIDSIZE], vector<Move>& moves, int color) {
    vector<pair<int, int>> scoredMoves;
    
    for(int i=0; i<(int)moves.size(); i++) {

        int score = scoreMoveHeuristic(board, moves[i], color);
        scoredMoves.push_back({score, i});
    }
    
    // 按分数降序排序
    sort(scoredMoves.begin(), scoredMoves.end(), [](const pair<int,int>& a, const pair<int,int>& b) {
        return a.first > b.first;
    });
    
    // 重新排列moves
    vector<Move> sortedMoves;
    for(auto& p : scoredMoves) {
        sortedMoves.push_back(moves[p.second]);
    }
    moves = sortedMoves;
}

// Alpha-Beta剪枝搜索（带置换表）
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

    // 获取所有合法走法
    vector<Move> moves = getAllMovesOnBoard(board, color);
    
    // 终止条件
    if(depth == 0 || moves.empty()) {
        int score = evaluateBoard(board, color);
        return score;
    }
    
    // 移动排序
    orderMoves(board, moves, color);
    
    int bestScore = -INF;
    Move localBestMove;
    int count = 0;
    
    for(const Move& move : moves) {
        if(count++ >= MAX_MOVES_PER_NODE) break; // 限制每个节点的最大走法数
        // 执行走法
        int newBoard[GRIDSIZE][GRIDSIZE];
        copyBoard(newBoard, board);
        if(!ProcStepOnBoard(newBoard, move.x0, move.y0, move.x1, move.y1, move.x2, move.y2, color)) {
            continue; // 防御性检查，理论上不应发生
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
