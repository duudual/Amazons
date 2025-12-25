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

int searchDepth = 4; // Alpha-Beta搜索深度
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
void orderMoves(int board[GRIDSIZE][GRIDSIZE], vector<Move>& moves, int color, const TTEntry* ttEntry) {
    vector<pair<int, int>> scoredMoves;
    
    for(int i=0; i<(int)moves.size(); i++) {
        int score = 0;
        
        // 如果是置换表中的最佳走法，给予最高优先级
        if(ttEntry && ttEntry->depth >= 0) {
            if(moves[i].x0 == ttEntry->bestMove.x0 && moves[i].y0 == ttEntry->bestMove.y0 &&
               moves[i].x1 == ttEntry->bestMove.x1 && moves[i].y1 == ttEntry->bestMove.y1 &&
               moves[i].x2 == ttEntry->bestMove.x2 && moves[i].y2 == ttEntry->bestMove.y2) {
                score = 1000000; // 最高优先级
            }
        }
        
        if(score == 0) {
            score = scoreMoveHeuristic(board, moves[i], color);
        }
        
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
    orderMoves(board, moves, color, ttEntry);
    
    int bestScore = -INF;
    Move localBestMove;
    
    for(const Move& move : moves) {
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
    initZobrist();
    
    Move bestMove;
    long long hash = computeHash(gridInfo, color);
    
    // 逐步增加深度
    for(int depth = 1; depth <= maxDepth; depth++) {
        Move currentBestMove;
        int score = alphaBetaSearch(gridInfo, depth, -INF, INF, color, currentBestMove, hash);
        
        if(currentBestMove.x0 >= 0) {
            bestMove = currentBestMove;
        }
        
        // 输出调试信息（可选）
        // cerr << "Depth " << depth << ": score = " << score << endl;
    }
    
    return bestMove;
}
// < -------- Alpha-Beta剪枝结束 --------- >

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
