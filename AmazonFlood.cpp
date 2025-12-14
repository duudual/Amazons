#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <cmath>
#include <random>
#include <chrono>
#include <cstring>


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

struct Move {
    int x0,y0,x1,y1,x2,y2;
    Move(int a=-1,int b=-1,int c=-1,int d=-1,int e=-1,int f=-1)
        : x0(a),y0(b),x1(c),y1(d),x2(e),y2(f) {}
};

// safe copy for board
inline void copyBoard(int dst[GRIDSIZE][GRIDSIZE], int src[GRIDSIZE][GRIDSIZE]) {
    memcpy(dst, src, sizeof(int) * GRIDSIZE * GRIDSIZE);
}

// safe ProcStepOnBoard (assumes board valid)
bool ProcStepOnBoard(int board[GRIDSIZE][GRIDSIZE], int x0, int y0, int x1, int y1, int x2, int y2, int color) {
    if (!inMap(x0,y0) || !inMap(x1,y1) || !inMap(x2,y2)) return false;
    if (board[x0][y0] != color) return false;
    if (board[x1][y1] != 0) return false;
    if ((board[x2][y2] != 0) && !(x2==x0 && y2==y0)) return false;
    board[x0][y0] = 0;
    board[x1][y1] = color;
    board[x2][y2] = OBSTACLE;
    return true;
}

// safe canReach on a given board (queen move)
bool canReachOnBoard(int board[GRIDSIZE][GRIDSIZE], int x1, int y1, int x2, int y2) {
    if (!inMap(x1,y1) || !inMap(x2,y2)) return false;
    if (x1==x2 && y1==y2) return true;
    int dx_dir = (x2==x1) ? 0 : (x2 > x1 ? 1 : -1);
    int dy_dir = (y2==y1) ? 0 : (y2 > y1 ? 1 : -1);
    if (dx_dir != 0 && dy_dir != 0) {
        if (abs(x2-x1) != abs(y2-y1)) return false;
    }
    if (dx_dir==0 && dy_dir==0) return false;
    int steps = max(abs(x2-x1), abs(y2-y1));
    for (int i=1;i<steps;i++){
        int x = x1 + dx_dir*i;
        int y = y1 + dy_dir*i;
        if (!inMap(x,y)) return false;
        if (board[x][y] != 0) return false;
    }
    return true;
}

// generate all legal moves on given board for color
vector<Move> getAllMovesOnBoard(int board[GRIDSIZE][GRIDSIZE], int color) {
    vector<Move> moves;
    for (int i=0;i<GRIDSIZE;i++){
        for (int j=0;j<GRIDSIZE;j++){
            if (board[i][j] != color) continue;
            // first ray (move queen)
            for (int k=0;k<8;k++){
                for (int delta1=1; delta1<GRIDSIZE; delta1++){
                    int xx = i + dx[k]*delta1;
                    int yy = j + dy[k]*delta1;
                    if (!inMap(xx,yy)) break;
                    if (board[xx][yy] != 0) break;
                    // second ray (shoot arrow from (xx,yy))
                    for (int l=0;l<8;l++){
                        for (int delta2=1; delta2<GRIDSIZE; delta2++){
                            int xxx = xx + dx[l]*delta2;
                            int yyy = yy + dy[l]*delta2;
                            if (!inMap(xxx,yyy)) break;
                            if (board[xxx][yyy] != 0 && !(i==xxx && j==yyy)) break;
                            // valid move
                            moves.emplace_back(i,j,xx,yy,xxx,yyy);
                        }
                    }
                }
            }
        }
    }
    return moves;
}

// fast reachable-count evaluation using BFS flood-fill per queen
// returns total reachable squares for all queens of 'color'
int reachableCount(int board[GRIDSIZE][GRIDSIZE], int color) {
    bool vis[GRIDSIZE][GRIDSIZE] = {0};
    int cnt = 0;
    queue<pair<int,int>> q;
    // push all queens of color as seeds
    for (int i=0;i<GRIDSIZE;i++){
        for (int j=0;j<GRIDSIZE;j++){
            if (board[i][j] == color) {
                if (!vis[i][j]) {
                    vis[i][j] = true;
                    q.push({i,j});
                }
            }
        }
    }
    // We treat queens as sources and BFS by queen-move adjacency (queen can jump along rays to all empty squares)
    while (!q.empty()) {
        auto [x,y] = q.front(); q.pop();
        // explore rays
        for (int dir=0;dir<8;dir++){
            for (int step=1; step<GRIDSIZE; step++){
                int nx = x + dx[dir]*step;
                int ny = y + dy[dir]*step;
                if (!inMap(nx,ny)) break;
                if (board[nx][ny] != 0) break;
                if (!vis[nx][ny]) {
                    vis[nx][ny] = true;
                    cnt++;
                    q.push({nx,ny});
                }
            }
        }
    }
    return cnt;
}

// score heuristic: difference of reachable squares, returns double
double quickReachScore(int board[GRIDSIZE][GRIDSIZE], int color) {
    int my = reachableCount(board, color);
    int opp = reachableCount(board, -color);
    // return difference normalized by board size
    return double(my - opp) / double(GRIDSIZE * GRIDSIZE);
}

// squash score into [-1,1] using tanh to stabilize
inline double squashScore(double raw) {
    return tanh(raw);
}

// ---------------- MCTS ----------------
struct MCTSNode {
    int board[GRIDSIZE][GRIDSIZE];
    Move move; // move that led to this node from parent
    int color; // player to move at this node
    double W;
    int N;
    bool fullyExpanded;
    vector<MCTSNode*> children;
    vector<Move> untriedMoves;
    MCTSNode* parent;
    MCTSNode(int b[GRIDSIZE][GRIDSIZE], int c, Move m, MCTSNode* p) : move(m), color(c), W(0), N(0), fullyExpanded(false), parent(p) {
        copyBoard(board, b);
    }
};

std::mt19937 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());

double UCTValue(MCTSNode* parent, MCTSNode* child, double C = 1.41421356237) {
    if (child->N == 0) {
        // large + tiny random to break ties
        double tie = std::uniform_real_distribution<double>(0.0, 1e-6)(rng);
        return 1e9 + tie;
    }
    double exploitation = child->W / child->N;
    double exploration = C * sqrt(log(double(parent->N + 1)) / double(child->N));
    return exploitation + exploration;
}

MCTSNode* selectByUCT(MCTSNode* node) {
    while (!node->children.empty()) {
        MCTSNode* best = nullptr;
        double bestUCT = -1e100;
        for (auto ch : node->children) {
            double u = UCTValue(node, ch);
            if (u > bestUCT) {
                bestUCT = u;
                best = ch;
            }
        }
        if (best == nullptr) break;
        node = best;
    }
    return node;
}

// expand: expand one untried move (standard MCTS)
MCTSNode* expand(MCTSNode* node) {
    if (node->untriedMoves.empty()) {
        // generate moves and set untriedMoves
        node->untriedMoves = getAllMovesOnBoard(node->board, node->color);
        if (node->untriedMoves.empty()) {
            node->fullyExpanded = true;
            return node; // terminal
        }
    }
    // pick one move at random from untriedMoves to expand
    int idx = std::uniform_int_distribution<int>(0, (int)node->untriedMoves.size()-1)(rng);
    Move mv = node->untriedMoves[idx];
    // remove it from untriedMoves (swap-pop)
    node->untriedMoves[idx] = node->untriedMoves.back();
    node->untriedMoves.pop_back();

    int newBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(newBoard, node->board);
    ProcStepOnBoard(newBoard, mv.x0, mv.y0, mv.x1, mv.y1, mv.x2, mv.y2, node->color);
    MCTSNode* child = new MCTSNode(newBoard, -node->color, mv, node);
    node->children.push_back(child);
    if (node->untriedMoves.empty()) node->fullyExpanded = true;
    return child;
}

// weighted selection among moves using quickReachScore; temperature controls randomness
Move selectMoveByWeight(const vector<Move>& moves, int board[GRIDSIZE][GRIDSIZE], int color, double temperature = 1.0) {
    if (moves.empty()) return Move();
    vector<double> scores;
    scores.reserve(moves.size());
    double minS = 1e100, maxS = -1e100;
    for (auto &m : moves) {
        int tmp[GRIDSIZE][GRIDSIZE];
        copyBoard(tmp, board);
        if (!ProcStepOnBoard(tmp, m.x0,m.y0,m.x1,m.y1,m.x2,m.y2, color)) {
            scores.push_back(-1e6);
            continue;
        }
        double s = quickReachScore(tmp, color);
        scores.push_back(s);
        minS = min(minS, s);
        maxS = max(maxS, s);
    }
    // softmax with temperature (normalize)
    vector<double> weights(moves.size());
    double sum = 0;
    for (size_t i=0;i<moves.size();i++) {
        double norm = (maxS - minS > 1e-9) ? (scores[i]-minS)/(maxS-minS) : 0.5;
        double w = exp(norm / max(1e-9, temperature));
        weights[i] = w;
        sum += w;
    }
    double r = std::uniform_real_distribution<double>(0.0, sum)(rng);
    double cumsum = 0;
    for (size_t i=0;i<moves.size();i++) {
        cumsum += weights[i];
        if (r <= cumsum) return moves[i];
    }
    return moves.back();
}

// simulate (rollout) using weighted random moves guided by quickReachScore
double simulate(int board[GRIDSIZE][GRIDSIZE], int rootColor, int currentColor, int maxDepth = 50) {
    int sim[GRIDSIZE][GRIDSIZE];
    copyBoard(sim, board);

    for (int d=0; d<maxDepth; d++){
        vector<Move> mv = getAllMovesOnBoard(sim, currentColor);
        if (mv.empty()) break;
        double temp = 0.4 + double(d) * 0.08; // gradually more random
        Move m = selectMoveByWeight(mv, sim, currentColor, temp);
        if (m.x0 < 0) break;
        ProcStepOnBoard(sim, m.x0,m.y0,m.x1,m.y1,m.x2,m.y2, currentColor);
        currentColor = -currentColor;
    }
    double raw = quickReachScore(sim, rootColor);
    // we want a reward perspective where positive is good for rootColor
    return squashScore(raw);
}

void backpropagate(MCTSNode* node, double value) {
    // value already from rootColor perspective
    while (node != nullptr) {
        node->N += 1;
        node->W += value;
        value = -value; // alternate perspective
        node = node->parent;
    }
}

void deleteTree(MCTSNode* node) {
    if (!node) return;
    for (auto ch : node->children) deleteTree(ch);
    delete node;
}

// Monte Carlo Search by iterations
Move monteCarloSearch(int color, int ITER, int rootBoardIn[GRIDSIZE][GRIDSIZE]) {
    int rootBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(rootBoard, rootBoardIn);

    MCTSNode* root = new MCTSNode(rootBoard, color, Move(), nullptr);
    // lazily fill root->untriedMoves on first expand/select
    for (int it=0; it<ITER; it++){
        MCTSNode* node = selectByUCT(root);
        // if node is not terminal, expand
        vector<Move> moves = getAllMovesOnBoard(node->board, node->color);
        if (!moves.empty() && !node->fullyExpanded) {
            node = expand(node);
        }
        // simulate
        double value = simulate(node->board, color, node->color, 20);
        // backprop
        backpropagate(node, value);
    }
    // pick best child by visit count
    MCTSNode* best = nullptr;
    int bestN = -1;
    for (auto ch : root->children) {
        if (ch->N > bestN) {
            bestN = ch->N;
            best = ch;
        }
    }
    Move result = Move();
    if (best) result = best->move;
    deleteTree(root);
    return result;
}

// Monte Carlo Search by time limit (ms)
Move monteCarloSearchTime(int color, int timeLimitMs, int rootBoardIn[GRIDSIZE][GRIDSIZE]) {
    int rootBoard[GRIDSIZE][GRIDSIZE];
    copyBoard(rootBoard, rootBoardIn);

    MCTSNode* root = new MCTSNode(rootBoard, color, Move(), nullptr);
    auto t0 = chrono::high_resolution_clock::now();
    int iterations = 0;
    while (true) {
        auto now = chrono::high_resolution_clock::now();
        int elapsed = (int)chrono::duration_cast<chrono::milliseconds>(now - t0).count();
        if (elapsed >= timeLimitMs) break;
        MCTSNode* node = selectByUCT(root);
        vector<Move> moves = getAllMovesOnBoard(node->board, node->color);
        if (!moves.empty() && !node->fullyExpanded) {
            node = expand(node);
        }
        double value = simulate(node->board, color, node->color, 30);
        backpropagate(node, value);
        iterations++;
    }
    // choose best child by visits
    MCTSNode* best = nullptr;
    int bestN = -1;
    for (auto ch : root->children) {
        if (ch->N > bestN) {
            bestN = ch->N;
            best = ch;
        }
    }
    Move result = Move();
    if (best) result = best->move;
    deleteTree(root);
    return result;
}



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
	Move bestMove = monteCarloSearch(currBotColor, 30, gridInfo); // 每个走法模拟30次
	
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