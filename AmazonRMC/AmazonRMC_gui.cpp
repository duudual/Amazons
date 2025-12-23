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
double simulateGame(MCTSNode* node,int maxDepth = maxDepth) {
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
Move monteCarloSearch(int color, int iterations = iterations, int maxDepth = maxDepth, int c = c) {
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
        double simulationResult = simulateGame(expandedNode,maxDepth);

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
// < -------- MCTS --------- >//

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

    // 使用 MCTS 算法进行决策
    Move bestMove = monteCarloSearch(aiColor, iterations,maxDepth,c);
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
