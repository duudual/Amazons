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
const int K1_QUEEN_MOVES = 500;   // 第一阶段保留的皇后移动位置数
const int K2_FULL_MOVES = 80;    // 第二阶段保留的完整走法数
const int INF = 1000000000; // 无穷大值

// 中盘阈值参数
const int MIDGAME_THRESHOLD = 20;  // 手数阈值，20手以后进入中盘
const int ENDGAME_DEPTH = 4;       // 中盘及以后的搜索深度
const int OPENING_DEPTH = 3;       // 开局阶段的搜索深度

// 性能分析统计
struct PerfStats {
    double totalphase1GenMoveTime = 0;
    double totalphase2GenMoveTime = 0;
    int phase1Count = 0;  // 第一阶段调用次数
    int phase2Count = 0;  // 第二阶段调用次数

    double totalTime = 0;  // 这次操作总耗时（毫秒）
    double evaluateBoardTime = 0; // 评估函数总耗时（毫秒）
    double evaluateMoveTime = 0;  // 评估走法总耗时（毫秒）
    int evaluateBoardCount = 0;   // 评估函数调用次数
    int evaluateMoveCount = 0;    // 评估走法调用次数
    double phase2SortMoveTime = 0; // 第二阶段走法排序耗时（毫秒）
    double phase2MoveinGenerateTime = 0; // 第二阶段走法生成耗时（毫秒）
    double phase2PrepMoveTime = 0; // 第二阶段预评估耗时（毫秒）
    double phase2BeforeMoveScore = 0;
    double phase2MoveScoreInF = 0;
    
    void reset() {
        totalphase1GenMoveTime = 0;
        totalphase2GenMoveTime = 0;
        phase1Count = 0;
        phase2Count = 0;
        totalTime = 0;
        evaluateBoardTime = 0;
        evaluateMoveTime = 0;
        evaluateBoardCount = 0;
        evaluateMoveCount = 0;
        phase2SortMoveTime =0;
        phase2MoveinGenerateTime = 0;
        phase2PrepMoveTime = 0;
        phase2BeforeMoveScore =0;
        phase2MoveScoreInF =0;
    }
    
    void printStats() const {
        cout << "\n========== time cost ==========\n";
        cout << "Phase 1 total: " << totalphase1GenMoveTime << " ms";
        if(phase1Count > 0) cout << " | avg: " << (totalphase1GenMoveTime / phase1Count) << " ms/call";
        cout << "\n";
        
        cout << "Phase 2 total: " << totalphase2GenMoveTime << " ms";
        if(phase2Count > 0) cout << " | avg: " << (totalphase2GenMoveTime / phase2Count) << " ms/call";
        cout << "\n";
        
        cout << "Evaluate Board total: " << evaluateBoardTime << " ms";
        if(evaluateBoardCount > 0) cout << " | avg: " << (evaluateBoardTime / evaluateBoardCount) << " ms/call";
        cout << "\n";
        
        cout << "Evaluate Move total: " << evaluateMoveTime << " ms";
        if(evaluateMoveCount > 0) cout << " | avg: " << (evaluateMoveTime / evaluateMoveCount) << " ms/call";
        cout << "\n";
        cout << "Phase 2 Move Preparation total: " << phase2PrepMoveTime << " ms\n";
        cout << "Phase 2 Move Generation total: " << phase2MoveinGenerateTime << " ms\n";
        cout << "Phase 2 Before Move Score total: " << phase2BeforeMoveScore << " ms\n";
        cout<<"phase2 Move Score in FullMove total: "<< phase2MoveScoreInF << " ms\n";
        cout << "Phase 2 Move Sorting total: " << phase2SortMoveTime << " ms\n";
        
        cout << "Total time: " << totalTime << " ms\n";
        cout << "================================\n";
    }
} perfStats;

// 权重参数
const int WEIGHT_CONNECTIVITY = 10;  // 联通点权重
const int WEIGHT_MOBILITY = 3;       // 灵活性权重
const int WEIGHT_CENTER = 1;         // 中心性权重

// 评估模式选择：0 = 原始评估，1 = 高级位棋盘评估
int EVAL_MODE = 1;

// ========== 高级位棋盘评估系统 ==========
typedef unsigned long long u64;

// 权重数组（根据空格数量索引）
double evalW1[56] = { 0,0,0.074938141,2.113308430,1.829644322,1.657095432,1.844633341,1.709531188,1.947853684,1.597983241,1.290230632,1.167080998,0.580433190,0.053694796,-0.011891359,0.527658403,0.975054741,1.056464434,0.872889936,0.678786278,0.323303133,0.346939683,0.260094881,0.264616042,0.246926725,0.189162090,0.141057417,0.105003797,0.100303024,0.094548225,0.075991668,0.056678012,0.052921258,0.046266042,0.047984783,0.029863806,0.046501521,0.035372451,0.036930695,0.026028842,0.023727236,0.008815981,0.005980607,0.009275969,-0.003742765,-0.009044983,-0.010549444,-0.032312561,-0.012371800,-0.037411232,-0.032170087,-0.012665736,-0.025661280,0.030489521,0.038260937,-0.033962481 };
double evalW2[56] = { 0,0,-0.071123712,1.850196242,1.909757733,1.786974669,1.441033125,1.831614137,1.863318324,1.549122095,0.929738343,0.985552371,0.519130468,0.652193844,0.677611768,0.715636075,-0.022837769,0.100610048,-0.172568828,-0.087590173,0.244123191,0.132074401,0.141869083,0.097190522,0.033922136,0.084543742,0.094968185,0.085143305,0.066933967,0.070909552,0.083291963,0.085179411,0.069383688,0.058226194,0.060640641,0.065368392,0.054329056,0.054176800,0.042966656,0.044218585,0.038709428,0.031496748,0.023774918,0.015536518,0.010674691,0.005355561,0.008374230,0.010649841,0.012298613,0.012731116,0.006807683,0.003942049,0.000137917,-0.014996326,-0.007956794,-0.035135169 };
double evalW3[56] = { 0,0,0.285573512,1.067316651,-0.211924806,0.680028081,0.517808735,0.304819107,-0.360255659,0.283211440,-0.118919544,-0.320547521,-0.164136052,-0.199558660,-0.158322647,-0.134399980,-0.166122779,-0.147537082,-0.119008608,-0.115313880,-0.088800281,-0.144199789,-0.093452334,-0.098844223,-0.061869688,-0.042582039,-0.026260521,-0.000397748,-0.004830216,-0.015821418,0.015736405,0.018657653,0.027041025,0.018133903,0.027038572,0.051898617,0.024094300,0.032939505,0.027083304,0.032979585,0.024715593,0.048389383,0.051429220,0.049864184,0.065294459,0.079495110,0.070997894,0.100243673,0.052135076,0.094104849,0.087983906,0.061890475,0.079814047,-0.026632605,-0.043000557,0.082739472 };
double evalW4[56] = { 0,0,0.061070204,2.073594332,1.803056121,1.969837666,1.993263960,1.986328840,1.706978559,1.858220458,1.816977739,1.320441842,1.554476857,1.503799796,1.165330052,0.306896836,0.428119272,0.127804548,0.320249408,0.250915736,0.208399042,0.209280223,0.189255610,0.176422641,0.158465251,0.123073861,0.114771605,0.107638910,0.104587719,0.094095841,0.077729441,0.085181594,0.088369705,0.096417032,0.086386517,0.084225371,0.094143234,0.101544440,0.116312958,0.118276447,0.126188114,0.131197393,0.142377511,0.138548672,0.131359026,0.133364707,0.122829571,0.130566165,0.130473971,0.133546650,0.111887053,0.097306497,0.095673442,0.109085456,0.139188647,0.08377216 };
double evalW5[56] = { 0,0,0.107754663,-0.528144240,-1.379346251,-1.006276369,-0.825565398,-0.430609971,-0.473824382,-0.060512420,-0.942406714,-0.967272699,-0.992075503,-1.073423266,-0.858994186,-0.838604629,-0.738691866,-0.771107376,-0.667914927,-0.620147049,-0.567148566,-0.658371270,-0.578625798,-0.549303949,-0.488953888,-0.356694996,-0.272865117,-0.213500351,-0.175519645,-0.133494422,-0.070355214,-0.038640279,-0.043525659,-0.037141897,-0.031321511,-0.023366986,-0.021908367,-0.013590249,-0.019244174,-0.014974004,-0.009036653,-0.018890575,-0.004873865,-0.017690081,-0.012989196,-0.020389291,-0.031792857,-0.031873308,-0.028690575,-0.024877874,-0.033742540,-0.026430298,-0.027800240,-0.030750951,-0.021546466,-0.037860770 };
double evalW6[56] = { 0,0,-4.567995548,-3.716784000,-3.381752491,-2.938109636,-3.290686846,-2.702058077,-2.855453491,-2.595350504,-2.005551338,-1.637469649,-1.272953510,-1.047218084,-0.879562616,-0.746166408,-0.717715025,-0.592138231,-0.504874051,-0.366238832,-0.404548317,-0.305050164,-0.272765279,-0.249806121,-0.137687027,-0.254965365,-0.089192919,-0.171341106,-0.076427318,-0.220489189,-0.091823421,-0.180632398,-0.106192604,-0.145687878,-0.073370934,-0.149500415,-0.069663063,-0.189439490,-0.027017180,-0.123841494,-0.011969413,-0.017853998,0.068564072,0.002616666,0.172109649,0.157967970,0.158981338,0.351862490,0.168605193,0.400329232,0.361347109,0.220473915,0.434309691,-0.030947627,-0.074839503,0.44456172 };

// 方向掩码（用于位棋盘操作）
const u64 DIR_NORTH = 0x00FFFFFFFFFFFFFFULL;  // 不含最上行
const u64 DIR_SOUTH = 0xFFFFFFFFFFFFFF00ULL;  // 不含最下行
const u64 DIR_WEST  = 0x7F7F7F7F7F7F7F7FULL;  // 不含最左列
const u64 DIR_EAST  = 0xFEFEFEFEFEFEFEFEULL;  // 不含最右列
const u64 DIR_NW    = DIR_NORTH & DIR_WEST;
const u64 DIR_NE    = DIR_NORTH & DIR_EAST;
const u64 DIR_SW    = DIR_SOUTH & DIR_WEST;
const u64 DIR_SE    = DIR_SOUTH & DIR_EAST;

// 位操作工具函数
inline int popcnt(u64 x) {
    int count = 0;
    while(x) { count++; x &= x - 1; }
    return count;
}

inline int lowestOneBit(u64 x) {
    if(x == 0) return -1;
    int pos = 0;
    while(!(x & 1)) { x >>= 1; pos++; }
    return pos;
}

// 高级评估器结构
struct AdvancedEvaluator {
    int mobValues[GRIDSIZE][GRIDSIZE];
    int qcnt1, qcnt2, kcnt1, kcnt2, maxq, maxk;
    u64 umap, p1, p2, blank;
    u64 uq[40], uk[40], q1[40], q2[40], k1[40], k2[40];
    double depthParameter[7] = { 0.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125 };
    double cntParameter[40] = { 0.0, 1.0, 0.95, 0.9, 0.8, 0.5, 0.5 };
    int player;  // 当前评估的玩家视角
    
    void init(int board[GRIDSIZE][GRIDSIZE], int myColor) {
        memset(mobValues, 0, sizeof(mobValues));
        memset(uq, 0, sizeof(uq));
        memset(uk, 0, sizeof(uk));
        memset(q1, 0, sizeof(q1));
        memset(q2, 0, sizeof(q2));
        memset(k1, 0, sizeof(k1));
        memset(k2, 0, sizeof(k2));
        qcnt1 = qcnt2 = kcnt1 = kcnt2 = maxq = maxk = 0;
        umap = p1 = p2 = blank = 0;
        player = myColor;
        
        // 将棋盘转换为位棋盘
        for(int i = 0; i < GRIDSIZE; i++) {
            for(int j = 0; j < GRIDSIZE; j++) {
                int bit = i * 8 + j;
                if(board[i][j] == 0) {
                    umap |= (1ULL << bit);  // 空格
                } else if(board[i][j] == player) {
                    p1 |= (1ULL << bit);    // 己方棋子
                } else if(board[i][j] == -player) {
                    p2 |= (1ULL << bit);    // 对方棋子
                }
                // 障碍物不设置任何位
            }
        }
    }
    
    // 沿方向扩展掩码
    u64 shiftMask(u64 a, int shift, u64 num, u64 direction, u64 can) {
        u64 result = a;
        for(int i = 0; i < 7; i++) {
            result |= (shift > 0 ? result << num : result >> num) & direction & can;
        }
        return result;
    }
    
    // 应用所有8个方向的移位
    u64 applyShifts(u64 temp, u64 mask) {
        u64 result = temp;
        // 8个方向：北、南、西、东、西北、东北、西南、东南
        result |= shiftMask(temp, -1, 1, DIR_WEST, mask);   // 西
        result |= shiftMask(temp, 1, 1, DIR_EAST, mask);    // 东
        result |= shiftMask(temp, -1, 8, DIR_NORTH, mask);  // 北
        result |= shiftMask(temp, 1, 8, DIR_SOUTH, mask);   // 南
        result |= shiftMask(temp, -1, 9, DIR_NW, mask);     // 西北
        result |= shiftMask(temp, 1, 9, DIR_SE, mask);      // 东南
        result |= shiftMask(temp, -1, 7, DIR_NE, mask);     // 东北
        result |= shiftMask(temp, 1, 7, DIR_SW, mask);      // 西南
        return result;
    }
    
    // 皇后BFS（沿皇后移动方向扩展）
    void queenBFS() {
        q1[0] |= p1;
        q2[0] |= p2;
        
        do {
            ++qcnt1;
            u64 temp = q1[qcnt1 - 1];
            q1[qcnt1] = applyShifts(temp, umap | p1);
        } while(q1[qcnt1] != q1[qcnt1 - 1] && qcnt1 < 38);
        
        for(int x = qcnt1 - 1; x >= 1; --x) q1[x] ^= q1[x - 1];
        
        do {
            ++qcnt2;
            u64 temp = q2[qcnt2 - 1];
            q2[qcnt2] = applyShifts(temp, umap | p2);
        } while(q2[qcnt2] != q2[qcnt2 - 1] && qcnt2 < 38);
        
        for(int x = qcnt2 - 1; x >= 1; --x) q2[x] ^= q2[x - 1];
    }
    
    // 国王BFS（一步移动扩展）
    void kingBFS() {
        k1[0] |= p1;
        k2[0] |= p2;
        
        do {
            ++kcnt1;
            u64 a = k1[kcnt1 - 1];
            a |= (a >> 1 & DIR_WEST) | (a << 1 & DIR_EAST) | 
                 (a >> 8 & DIR_NORTH) | (a << 8 & DIR_SOUTH) |
                 (a >> 9 & DIR_NW) | (a << 9 & DIR_SE) | 
                 (a >> 7 & DIR_NE) | (a << 7 & DIR_SW);
            k1[kcnt1] = a & (umap | p1);
        } while(k1[kcnt1] != k1[kcnt1 - 1] && kcnt1 < 38);
        
        for(int x = kcnt1 - 1; x >= 1; --x) k1[x] ^= k1[x - 1];
        
        do {
            ++kcnt2;
            u64 a = k2[kcnt2 - 1];
            a |= (a >> 1 & DIR_WEST) | (a << 1 & DIR_EAST) | 
                 (a >> 8 & DIR_NORTH) | (a << 8 & DIR_SOUTH) |
                 (a >> 9 & DIR_NW) | (a << 9 & DIR_SE) | 
                 (a >> 7 & DIR_NE) | (a << 7 & DIR_SW);
            k2[kcnt2] = a & (umap | p2);
        } while(k2[kcnt2] != k2[kcnt2 - 1] && kcnt2 < 38);
        
        for(int x = kcnt2 - 1; x >= 1; --x) k2[x] ^= k2[x - 1];
    }
    
    // 预处理
    void pretreatment() {
        queenBFS();
        kingBFS();
        maxq = max(qcnt1, qcnt2);
        maxk = max(kcnt1, kcnt2);
        
        for(int i = 1; i < min(qcnt1, qcnt2); ++i) uq[i] |= q1[i] | q2[i];
        for(int i = min(qcnt1, qcnt2); i < maxq; ++i) {
            if(qcnt1 < qcnt2) uq[i] |= q2[i];
            else uq[i] |= q1[i];
        }
        for(int i = 2; i <= maxq; ++i) uq[i] |= uq[i - 1];
        
        for(int i = 1; i < min(kcnt1, kcnt2); ++i) uk[i] |= k1[i] | k2[i];
        for(int i = min(kcnt1, kcnt2); i < maxk; ++i) {
            if(kcnt1 < kcnt2) uk[i] |= k2[i];
            else uk[i] |= k1[i];
        }
        for(int i = 2; i <= maxk; ++i) uk[i] |= uk[i - 1];
    }
    
    // 计算参数
    void computeParameter(double& w1Out, double& w2Out, double& w3Out, double& w4Out) {
        w1Out = w2Out = w3Out = w4Out = 0.0;
        
        // 皇后移动参数
        for(int i = 2; i < min(qcnt1, qcnt2); ++i) {
            w1Out += cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q1[i]);
            w1Out -= cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q2[i]);
        }
        for(int i = min(qcnt1, qcnt2); (qcnt1 != qcnt2) && (i < maxq); ++i) {
            if(qcnt1 < qcnt2) w1Out -= cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q2[i]);
            else w1Out += cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q1[i]);
        }
        w1Out += popcnt(uq[maxq] & (umap ^ q1[qcnt1])) - popcnt(uq[maxq] & (umap ^ q2[qcnt2])) + 0.3 * popcnt(q1[1] & q2[1]);
        
        // 国王移动参数
        for(int i = 2; i < min(kcnt1, kcnt2); ++i) {
            w2Out += cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k1[i]);
            w2Out -= cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k2[i]);
        }
        for(int i = min(kcnt1, kcnt2); (kcnt1 != kcnt2) && (i < maxk); ++i) {
            if(kcnt1 < kcnt2) w2Out -= cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k2[i]);
            else w2Out += cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k1[i]);
        }
        w2Out += popcnt(uk[maxk] & (umap ^ k1[kcnt1])) - popcnt(uk[maxk] & (umap ^ k2[kcnt2])) + 0.3 * popcnt(k1[1] & k2[1]);
        
        // 深度参数
        for(int i = 1; i <= 6 && i < qcnt1; ++i) w3Out -= popcnt(q1[i]) * depthParameter[i];
        for(int i = 1; i <= 6 && i < qcnt2; ++i) w3Out += popcnt(q2[i]) * depthParameter[i];
        
        // 领地参数
        u64 initial1 = k1[kcnt1], initial2 = k2[kcnt2];
        u64 intersect = initial1 & initial2;
        w4Out += popcnt(initial2 & ~initial1) - popcnt(initial1 & ~initial2);
        for(int i = 1; i < kcnt1; ++i) w4Out += popcnt(intersect & k1[i]) * i / 6.0;
        for(int i = 1; i < kcnt2; ++i) w4Out -= popcnt(intersect & k2[i]) * i / 6.0;
    }
    
    // 计算空白值
    double computeBlankValue() {
        double wv1 = 0, wv2 = 0;
        
        // 计算每个空格周围的空格数
        for(blank = umap; blank; blank &= blank - 1) {
            int bit = lowestOneBit(blank);
            u64 ubit = 1ULL << bit;
            ubit |= (ubit >> 1 & DIR_WEST) | (ubit << 1 & DIR_EAST) | 
                    (ubit >> 8 & DIR_NORTH) | (ubit << 8 & DIR_SOUTH) |
                    (ubit >> 9 & DIR_NW) | (ubit << 9 & DIR_SE) | 
                    (ubit >> 7 & DIR_NE) | (ubit << 7 & DIR_SW);
            int cnt = popcnt(ubit & umap);
            int row = bit / 8, column = bit % 8;
            mobValues[row][column] = cnt - 1;
        }
        
        // 对己方棋子计算空白值
        for(u64 piece = p1; piece; piece &= piece - 1) {
            int bit = lowestOneBit(piece);
            int cnt = 0;
            u64 upiece = 1ULL << bit;
            u64 a[40] = { 0ULL };
            a[0] |= upiece;
            do {
                ++cnt;
                u64 result = a[cnt - 1];
                result |= (result >> 1 & DIR_WEST) | (result << 1 & DIR_EAST) | 
                         (result >> 8 & DIR_NORTH) | (result << 8 & DIR_SOUTH) |
                         (result >> 9 & DIR_NW) | (result << 9 & DIR_SE) | 
                         (result >> 7 & DIR_NE) | (result << 7 & DIR_SW);
                a[cnt] = (result & umap) ^ upiece;
            } while(a[cnt] != a[cnt - 1] && cnt < 38);
            
            for(int x = cnt - 1; x >= 1; --x) a[x] ^= a[x - 1];
            for(int i = 1; i < cnt; ++i) {
                for(u64 temp = a[i]; temp; temp &= temp - 1) {
                    int b = lowestOneBit(temp);
                    int row = b / 8, column = b % 8;
                    wv1 += mobValues[row][column] * cntParameter[min(i, 6)] / i;
                }
            }
        }
        
        // 对对方棋子计算空白值
        for(u64 piece = p2; piece; piece &= piece - 1) {
            int bit = lowestOneBit(piece);
            int cnt = 0;
            u64 upiece = 1ULL << bit;
            u64 a[40] = { 0ULL };
            a[0] |= upiece;
            do {
                ++cnt;
                u64 result = a[cnt - 1];
                result |= (result >> 1 & DIR_WEST) | (result << 1 & DIR_EAST) | 
                         (result >> 8 & DIR_NORTH) | (result << 8 & DIR_SOUTH) |
                         (result >> 9 & DIR_NW) | (result << 9 & DIR_SE) | 
                         (result >> 7 & DIR_NE) | (result << 7 & DIR_SW);
                a[cnt] = (result & umap) ^ upiece;
            } while(a[cnt] != a[cnt - 1] && cnt < 38);
            
            for(int x = cnt - 1; x >= 1; --x) a[x] ^= a[x - 1];
            for(int i = 1; i < cnt; ++i) {
                for(u64 temp = a[i]; temp; temp &= temp - 1) {
                    int b = lowestOneBit(temp);
                    int row = b / 8, column = b % 8;
                    wv2 += mobValues[row][column] * cntParameter[min(i, 6)] / i;
                }
            }
        }
        
        return wv1 - wv2;
    }
    
    // 获取评估值
    double getEvaluateValue(int board[GRIDSIZE][GRIDSIZE], int myColor) {
        init(board, myColor);
        
        int blankCnts = popcnt(umap);
        if(blankCnts < 2) blankCnts = 2;
        if(blankCnts > 55) blankCnts = 55;
        
        double p1v = evalW1[blankCnts];
        double p2v = evalW2[blankCnts];
        double p3v = evalW3[blankCnts];
        double p4v = evalW4[blankCnts];
        double p5v = evalW5[blankCnts];
        double p6v = evalW6[blankCnts];
        
        pretreatment();
        
        double a, b, c, d;
        computeParameter(a, b, c, d);
        double e = 0.1 * computeBlankValue();
        
        double v = a * p1v + b * p2v + c * p3v + d * p4v + e * p5v + p6v;
        
        // 归一化到 [0, 1]，然后转换为整数分数
        double normalized = 1.0 - (1.0 / (1.0 + exp(-v)));
        
        // 将 [0, 1] 映射到 [-10000, 10000] 的整数分数
        return (normalized - 0.5) * 20000;
    }
};

// 高级评估函数（供alpha-beta调用）
int evaluateBoardAdvanced(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    auto startTime = chrono::high_resolution_clock::now();
    
    AdvancedEvaluator evaluator;
    int score = (int)evaluator.getEvaluateValue(board, myColor);
    
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> evalDuration = endTime - startTime;
    perfStats.evaluateBoardTime += evalDuration.count();
    perfStats.evaluateBoardCount++;
    
    return score;
}

// ========== 高级位棋盘评估系统结束 ==========


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

// 计算当前棋盘的走法手数（已执行的手数）
int countMovesMade(int board[GRIDSIZE][GRIDSIZE]) {
    int obstacleCount = 0;
    // 统计棋盘上的障碍物个数
    for(int i = 0; i < GRIDSIZE; i++) {
        for(int j = 0; j < GRIDSIZE; j++) {
            if(board[i][j] == OBSTACLE) {
                obstacleCount++;
            }
        }
    }
    // 每一手走法会放置一个障碍物
    return obstacleCount;
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
    auto startTime = chrono::high_resolution_clock::now();

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
    
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> evalDuration = endTime - startTime;
    perfStats.evaluateBoardTime += evalDuration.count();
    perfStats.evaluateBoardCount++;
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
    return blockCount * 2 + oppCount * 5;
}

// 第二阶段：从选定的皇后移动生成完整走法，并筛选前K2个
vector<Move> generateTopFullMoves(int board[GRIDSIZE][GRIDSIZE], 
                                   const vector<QueenMove>& queenMoves, 
                                   int color, int k2) {
    vector<pair<int, Move>> scoredMoves;  // (score, move)
    scoredMoves.reserve(queenMoves.size() * 8 * (GRIDSIZE - 1));
    
    for(const auto& qm : queenMoves) {
        // 逻辑占用查询：避免为每个皇后复制棋盘
        auto cellAt = [&](int x, int y) {
            if(x == qm.x0 && y == qm.y0) return 0;      // 原位置视为空
            if(x == qm.x1 && y == qm.y1) return color;  // 新位置视为己方
            return board[x][y];
            
        };

        // 生成所有可能的射箭位置
        for(int k = 0; k < 8; k++) {
            for(int delta = 1; delta < GRIDSIZE; delta++) {
                const int x2 = qm.x1 + dx[k] * delta;
                const int y2 = qm.y1 + dy[k] * delta;
                if(!inMap(x2, y2)) break;

                // 射线遇到非空则终止；原位置已被视为0，无需特判
                if(cellAt(x2, y2) != 0) break;

                Move fullMove(qm.x0, qm.y0, qm.x1, qm.y1, x2, y2);
                
                int score = scoreFullMove(board, fullMove, color);
                scoredMoves.push_back({score, fullMove});
            }
        }
    }
    // 仅截取前K2：nth_element + 局部排序，避免全量排序
    vector<Move> result;
    if(scoredMoves.empty()) return result;

    if((int)scoredMoves.size() > k2) {
        nth_element(scoredMoves.begin(), scoredMoves.begin() + k2, scoredMoves.end(),
                    [](const pair<int,Move>& a, const pair<int,Move>& b) { return a.first > b.first; });
        scoredMoves.resize(k2);
    }
    sort(scoredMoves.begin(), scoredMoves.end(),
         [](const pair<int,Move>& a, const pair<int,Move>& b) { return a.first > b.first; });

    result.reserve(scoredMoves.size());
    for(auto& item : scoredMoves) {
        result.push_back(item.second);
    }
    return result;
}

// 第二阶段：直接返回
vector<Move> generateTopFullMovesSim(int board[GRIDSIZE][GRIDSIZE], 
                                   const vector<QueenMove>& queenMoves, 
                                   int color, int k2) {
    vector<Move> scoredMoves;  // (score, move)
    int sum = 0;
    for(const auto& qm : queenMoves) {
        // 逻辑占用查询：避免为每个皇后复制棋盘
        auto cellAt = [&](int x, int y) {
            if(x == qm.x0 && y == qm.y0) return 0;      // 原位置视为空
            if(x == qm.x1 && y == qm.y1) return color;  // 新位置视为己方
            return board[x][y];
            
        };

        // 生成所有可能的射箭位置
        for(int k = 0; k < 8; k++) {
            for(int delta = 1; delta < GRIDSIZE; delta++) {
                const int x2 = qm.x1 + dx[k] * delta;
                const int y2 = qm.y1 + dy[k] * delta;
                if(!inMap(x2, y2)) break;
                // 射线遇到非空则终止；原位置已被视为0，无需特判
                if(cellAt(x2, y2) != 0) break;
                Move fullMove(qm.x0, qm.y0, qm.x1, qm.y1, x2, y2);
                scoredMoves.push_back(fullMove);
                sum++;
                if(sum >= k2) {
                    return scoredMoves;
                }
            }
        }
    }
    
    return scoredMoves;
}

// 两阶段走法生成主函数
vector<Move> generateMovesWithTwoPhase(int board[GRIDSIZE][GRIDSIZE], int color) {
    auto phase1StartTime = chrono::high_resolution_clock::now();
    
    // 第一阶段：生成并筛选皇后移动位置
    vector<QueenMove> topQueenMoves = generateTopQueenMoves(board, color, K1_QUEEN_MOVES);
    if(topQueenMoves.empty()) {
        return vector<Move>();  // 无法移动
    }
    auto phase1EndTime = chrono::high_resolution_clock::now();
    auto phase1Duration = chrono::duration_cast<chrono::microseconds>(phase1EndTime - phase1StartTime).count();
    perfStats.totalphase1GenMoveTime += phase1Duration / 1000.0;  // 转换为毫秒
    perfStats.phase1Count++;
    
    // 第二阶段：生成并筛选完整走法
    vector<Move> topFullMoves = generateTopFullMoves(board, topQueenMoves, color, K2_FULL_MOVES);

    auto phaseEndTime = chrono::high_resolution_clock::now();
    auto phaseDuration = chrono::duration_cast<chrono::microseconds>(phaseEndTime - phase1EndTime).count();
    perfStats.totalphase2GenMoveTime += phaseDuration / 1000.0;  // 转换为毫秒
    perfStats.phase2Count++;
    
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
        // 根据评估模式选择评估函数
        int score = (EVAL_MODE == 1) ? evaluateBoardAdvanced(board, color) : evaluateBoard(board, color);
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

    // 根据当前手数动态调整搜索深度
    int moveCount = countMovesMade(gridInfo);
    int currentDepth = (moveCount >= MIDGAME_THRESHOLD) ? ENDGAME_DEPTH : OPENING_DEPTH;
    cout << "[Move count: " << moveCount << ", Search depth: " << currentDepth << "]\n";

    // 使用 Alpha-Beta 剪枝算法进行决策
    Move bestMove = iterativeDeepeningSearch(aiColor, currentDepth);
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
