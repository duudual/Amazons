#define UNICODE
#define _UNICODE

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <limits>
#include <chrono>

#define GRIDSIZE 8
#define OBSTACLE 2
#define judge_black 0
#define judge_white 1
#define grid_black 1
#define grid_white -1

using namespace std;

int currBotColor;
int gridInfo[GRIDSIZE][GRIDSIZE] = { 0 };
int dx[] = { -1,-1,-1,0,0,1,1,1 };
int dy[] = { -1,0,1,-1,1,-1,0,1 };

const int K1_QUEEN_MOVES = 500;
const int K2_FULL_MOVES = 50;    // 第二阶段保留的完整走法数（halfnewScore版本）
const int INF = 1000000000;

const int MIDGAME_THRESHOLD = 28;
const int ENDGAME_DEPTH = 4;
const int OPENING_DEPTH = 3;

struct PerfStats {
    double totalphase1GenMoveTime = 0;
    double totalphase2GenMoveTime = 0;
    int phase1Count = 0;
    int phase2Count = 0;
    double totalTime = 0;
    double evaluateBoardTime = 0;
    int evaluateBoardCount = 0;
    
    void reset() {
        totalphase1GenMoveTime = 0;
        totalphase2GenMoveTime = 0;
        phase1Count = 0;
        phase2Count = 0;
        totalTime = 0;
        evaluateBoardTime = 0;
        evaluateBoardCount = 0;
    }
    
    void printStats() const {
        cerr << "\n========== time cost ==========\n";
        cerr << "Phase 1 total: " << totalphase1GenMoveTime << " ms";
        if(phase1Count > 0) cerr << " | avg: " << (totalphase1GenMoveTime / phase1Count) << " ms/call";
        cerr << "\n";
        cerr << "Phase 2 total: " << totalphase2GenMoveTime << " ms";
        if(phase2Count > 0) cerr << " | avg: " << (totalphase2GenMoveTime / phase2Count) << " ms/call";
        cerr << "\n";
        cerr << "Evaluate Board total: " << evaluateBoardTime << " ms";
        if(evaluateBoardCount > 0) cerr << " | avg: " << (evaluateBoardTime / evaluateBoardCount) << " ms/call";
        cerr << "\n";
        cerr << "Total time: " << totalTime << " ms\n";
        cerr << "================================\n";
    }
} perfStats;

const int WEIGHT_CONNECTIVITY = 10;
const int WEIGHT_MOBILITY = 3;
const int WEIGHT_CENTER = 1;

int EVAL_MODE = 1;

typedef unsigned long long u64;

double evalW1[56] = { 0,0,0.074938141,2.113308430,1.829644322,1.657095432,1.844633341,1.709531188,1.947853684,1.597983241,1.290230632,1.167080998,0.580433190,0.053694796,-0.011891359,0.527658403,0.975054741,1.056464434,0.872889936,0.678786278,0.323303133,0.346939683,0.260094881,0.264616042,0.246926725,0.189162090,0.141057417,0.105003797,0.100303024,0.094548225,0.075991668,0.056678012,0.052921258,0.046266042,0.047984783,0.029863806,0.046501521,0.035372451,0.036930695,0.026028842,0.023727236,0.008815981,0.005980607,0.009275969,-0.003742765,-0.009044983,-0.010549444,-0.032312561,-0.012371800,-0.037411232,-0.032170087,-0.012665736,-0.025661280,0.030489521,0.038260937,-0.033962481 };
double evalW2[56] = { 0,0,-0.071123712,1.850196242,1.909757733,1.786974669,1.441033125,1.831614137,1.863318324,1.549122095,0.929738343,0.985552371,0.519130468,0.652193844,0.677611768,0.715636075,-0.022837769,0.100610048,-0.172568828,-0.087590173,0.244123191,0.132074401,0.141869083,0.097190522,0.033922136,0.084543742,0.094968185,0.085143305,0.066933967,0.070909552,0.083291963,0.085179411,0.069383688,0.058226194,0.060640641,0.065368392,0.054329056,0.054176800,0.042966656,0.044218585,0.038709428,0.031496748,0.023774918,0.015536518,0.010674691,0.005355561,0.008374230,0.010649841,0.012298613,0.012731116,0.006807683,0.003942049,0.000137917,-0.014996326,-0.007956794,-0.035135169 };
double evalW3[56] = { 0,0,0.285573512,1.067316651,-0.211924806,0.680028081,0.517808735,0.304819107,-0.360255659,0.283211440,-0.118919544,-0.320547521,-0.164136052,-0.199558660,-0.158322647,-0.134399980,-0.166122779,-0.147537082,-0.119008608,-0.115313880,-0.088800281,-0.144199789,-0.093452334,-0.098844223,-0.061869688,-0.042582039,-0.026260521,-0.000397748,-0.004830216,-0.015821418,0.015736405,0.018657653,0.027041025,0.018133903,0.027038572,0.051898617,0.024094300,0.032939505,0.027083304,0.032979585,0.024715593,0.048389383,0.051429220,0.049864184,0.065294459,0.079495110,0.070997894,0.100243673,0.052135076,0.094104849,0.087983906,0.061890475,0.079814047,-0.026632605,-0.043000557,0.082739472 };
double evalW4[56] = { 0,0,0.061070204,2.073594332,1.803056121,1.969837666,1.993263960,1.986328840,1.706978559,1.858220458,1.816977739,1.320441842,1.554476857,1.503799796,1.165330052,0.306896836,0.428119272,0.127804548,0.320249408,0.250915736,0.208399042,0.209280223,0.189255610,0.176422641,0.158465251,0.123073861,0.114771605,0.107638910,0.104587719,0.094095841,0.077729441,0.085181594,0.088369705,0.096417032,0.086386517,0.084225371,0.094143234,0.101544440,0.116312958,0.118276447,0.126188114,0.131197393,0.142377511,0.138548672,0.131359026,0.133364707,0.122829571,0.130566165,0.130473971,0.133546650,0.111887053,0.097306497,0.095673442,0.109085456,0.139188647,0.08377216 };
double evalW5[56] = { 0,0,0.107754663,-0.528144240,-1.379346251,-1.006276369,-0.825565398,-0.430609971,-0.473824382,-0.060512420,-0.942406714,-0.967272699,-0.992075503,-1.073423266,-0.858994186,-0.838604629,-0.738691866,-0.771107376,-0.667914927,-0.620147049,-0.567148566,-0.658371270,-0.578625798,-0.549303949,-0.488953888,-0.356694996,-0.272865117,-0.213500351,-0.175519645,-0.133494422,-0.070355214,-0.038640279,-0.043525659,-0.037141897,-0.031321511,-0.023366986,-0.021908367,-0.013590249,-0.019244174,-0.014974004,-0.009036653,-0.018890575,-0.004873865,-0.017690081,-0.012989196,-0.020389291,-0.031792857,-0.031873308,-0.028690575,-0.024877874,-0.033742540,-0.026430298,-0.027800240,-0.030750951,-0.021546466,-0.037860770 };
double evalW6[56] = { 0,0,-4.567995548,-3.716784000,-3.381752491,-2.938109636,-3.290686846,-2.702058077,-2.855453491,-2.595350504,-2.005551338,-1.637469649,-1.272953510,-1.047218084,-0.879562616,-0.746166408,-0.717715025,-0.592138231,-0.504874051,-0.366238832,-0.404548317,-0.305050164,-0.272765279,-0.249806121,-0.137687027,-0.254965365,-0.089192919,-0.171341106,-0.076427318,-0.220489189,-0.091823421,-0.180632398,-0.106192604,-0.145687878,-0.073370934,-0.149500415,-0.069663063,-0.189439490,-0.027017180,-0.123841494,-0.011969413,-0.017853998,0.068564072,0.002616666,0.172109649,0.157967970,0.158981338,0.351862490,0.168605193,0.400329232,0.361347109,0.220473915,0.434309691,-0.030947627,-0.074839503,0.44456172 };

const u64 DIR_NORTH = 0x00FFFFFFFFFFFFFFULL;
const u64 DIR_SOUTH = 0xFFFFFFFFFFFFFF00ULL;
const u64 DIR_WEST  = 0x7F7F7F7F7F7F7F7FULL;
const u64 DIR_EAST  = 0xFEFEFEFEFEFEFEFEULL;
const u64 DIR_NW    = DIR_NORTH & DIR_WEST;
const u64 DIR_NE    = DIR_NORTH & DIR_EAST;
const u64 DIR_SW    = DIR_SOUTH & DIR_WEST;
const u64 DIR_SE    = DIR_SOUTH & DIR_EAST;

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

struct AdvancedEvaluator {
    int mobValues[GRIDSIZE][GRIDSIZE];
    int qcnt1, qcnt2, kcnt1, kcnt2, maxq, maxk;
    u64 umap, p1, p2, blank;
    u64 uq[40], uk[40], q1[40], q2[40], k1[40], k2[40];
    double depthParameter[7] = { 0.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125 };
    double cntParameter[40] = { 0.0, 1.0, 0.95, 0.9, 0.8, 0.5, 0.5 };
    int player;
    
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
        for(int i = 0; i < GRIDSIZE; i++) {
            for(int j = 0; j < GRIDSIZE; j++) {
                int bit = i * 8 + j;
                if(board[i][j] == 0) umap |= (1ULL << bit);
                else if(board[i][j] == player) p1 |= (1ULL << bit);
                else if(board[i][j] == -player) p2 |= (1ULL << bit);
            }
        }
    }
    
    u64 shiftMask(u64 a, int shift, u64 num, u64 direction, u64 can) {
        u64 result = a;
        for(int i = 0; i < 7; i++) {
            result |= (shift > 0 ? result << num : result >> num) & direction & can;
        }
        return result;
    }
    
    u64 applyShifts(u64 temp, u64 mask) {
        u64 result = temp;
        result |= shiftMask(temp, -1, 1, DIR_WEST, mask);
        result |= shiftMask(temp, 1, 1, DIR_EAST, mask);
        result |= shiftMask(temp, -1, 8, DIR_NORTH, mask);
        result |= shiftMask(temp, 1, 8, DIR_SOUTH, mask);
        result |= shiftMask(temp, -1, 9, DIR_NW, mask);
        result |= shiftMask(temp, 1, 9, DIR_SE, mask);
        result |= shiftMask(temp, -1, 7, DIR_NE, mask);
        result |= shiftMask(temp, 1, 7, DIR_SW, mask);
        return result;
    }
    
    void queenBFS() {
        q1[0] |= p1; q2[0] |= p2;
        do {
            ++qcnt1;
            q1[qcnt1] = applyShifts(q1[qcnt1 - 1], umap | p1);
        } while(q1[qcnt1] != q1[qcnt1 - 1] && qcnt1 < 38);
        for(int x = qcnt1 - 1; x >= 1; --x) q1[x] ^= q1[x - 1];
        do {
            ++qcnt2;
            q2[qcnt2] = applyShifts(q2[qcnt2 - 1], umap | p2);
        } while(q2[qcnt2] != q2[qcnt2 - 1] && qcnt2 < 38);
        for(int x = qcnt2 - 1; x >= 1; --x) q2[x] ^= q2[x - 1];
    }
    
    void kingBFS() {
        k1[0] |= p1; k2[0] |= p2;
        do {
            ++kcnt1;
            u64 a = k1[kcnt1 - 1];
            a |= (a >> 1 & DIR_WEST) | (a << 1 & DIR_EAST) | (a >> 8 & DIR_NORTH) | (a << 8 & DIR_SOUTH) | (a >> 9 & DIR_NW) | (a << 9 & DIR_SE) | (a >> 7 & DIR_NE) | (a << 7 & DIR_SW);
            k1[kcnt1] = a & (umap | p1);
        } while(k1[kcnt1] != k1[kcnt1 - 1] && kcnt1 < 38);
        for(int x = kcnt1 - 1; x >= 1; --x) k1[x] ^= k1[x - 1];
        do {
            ++kcnt2;
            u64 a = k2[kcnt2 - 1];
            a |= (a >> 1 & DIR_WEST) | (a << 1 & DIR_EAST) | (a >> 8 & DIR_NORTH) | (a << 8 & DIR_SOUTH) | (a >> 9 & DIR_NW) | (a << 9 & DIR_SE) | (a >> 7 & DIR_NE) | (a << 7 & DIR_SW);
            k2[kcnt2] = a & (umap | p2);
        } while(k2[kcnt2] != k2[kcnt2 - 1] && kcnt2 < 38);
        for(int x = kcnt2 - 1; x >= 1; --x) k2[x] ^= k2[x - 1];
    }
    
    void pretreatment() {
        queenBFS(); kingBFS();
        maxq = max(qcnt1, qcnt2); maxk = max(kcnt1, kcnt2);
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
    
    void computeParameter(double& w1Out, double& w2Out, double& w3Out, double& w4Out) {
        w1Out = w2Out = w3Out = w4Out = 0.0;
        for(int i = 2; i < min(qcnt1, qcnt2); ++i) {
            w1Out += cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q1[i]);
            w1Out -= cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q2[i]);
        }
        for(int i = min(qcnt1, qcnt2); (qcnt1 != qcnt2) && (i < maxq); ++i) {
            if(qcnt1 < qcnt2) w1Out -= cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q2[i]);
            else w1Out += cntParameter[min(i-1, 6)] * popcnt(uq[i - 1] & q1[i]);
        }
        w1Out += popcnt(uq[maxq] & (umap ^ q1[qcnt1])) - popcnt(uq[maxq] & (umap ^ q2[qcnt2])) + 0.3 * popcnt(q1[1] & q2[1]);
        for(int i = 2; i < min(kcnt1, kcnt2); ++i) {
            w2Out += cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k1[i]);
            w2Out -= cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k2[i]);
        }
        for(int i = min(kcnt1, kcnt2); (kcnt1 != kcnt2) && (i < maxk); ++i) {
            if(kcnt1 < kcnt2) w2Out -= cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k2[i]);
            else w2Out += cntParameter[min(i-1, 6)] * popcnt(uk[i - 1] & k1[i]);
        }
        w2Out += popcnt(uk[maxk] & (umap ^ k1[kcnt1])) - popcnt(uk[maxk] & (umap ^ k2[kcnt2])) + 0.3 * popcnt(k1[1] & k2[1]);
        for(int i = 1; i <= 6 && i < qcnt1; ++i) w3Out -= popcnt(q1[i]) * depthParameter[i];
        for(int i = 1; i <= 6 && i < qcnt2; ++i) w3Out += popcnt(q2[i]) * depthParameter[i];
        u64 initial1 = k1[kcnt1], initial2 = k2[kcnt2], intersect = initial1 & initial2;
        w4Out += popcnt(initial2 & ~initial1) - popcnt(initial1 & ~initial2);
        for(int i = 1; i < kcnt1; ++i) w4Out += popcnt(intersect & k1[i]) * i / 6.0;
        for(int i = 1; i < kcnt2; ++i) w4Out -= popcnt(intersect & k2[i]) * i / 6.0;
    }
    
    double computeBlankValue() {
        double wv1 = 0, wv2 = 0;
        for(blank = umap; blank; blank &= blank - 1) {
            int bit = lowestOneBit(blank);
            u64 ubit = 1ULL << bit;
            ubit |= (ubit >> 1 & DIR_WEST) | (ubit << 1 & DIR_EAST) | (ubit >> 8 & DIR_NORTH) | (ubit << 8 & DIR_SOUTH) | (ubit >> 9 & DIR_NW) | (ubit << 9 & DIR_SE) | (ubit >> 7 & DIR_NE) | (ubit << 7 & DIR_SW);
            int cnt = popcnt(ubit & umap);
            mobValues[bit / 8][bit % 8] = cnt - 1;
        }
        auto calcWv = [&](u64 pieces) {
            double val = 0;
            for(u64 piece = pieces; piece; piece &= piece - 1) {
                int bit = lowestOneBit(piece), cnt = 0;
                u64 upiece = 1ULL << bit, a[40] = { 0ULL };
                a[0] |= upiece;
                do {
                    ++cnt;
                    u64 res = a[cnt - 1];
                    res |= (res >> 1 & DIR_WEST) | (res << 1 & DIR_EAST) | (res >> 8 & DIR_NORTH) | (res << 8 & DIR_SOUTH) | (res >> 9 & DIR_NW) | (res << 9 & DIR_SE) | (res >> 7 & DIR_NE) | (res << 7 & DIR_SW);
                    a[cnt] = (res & umap) ^ upiece;
                } while(a[cnt] != a[cnt - 1] && cnt < 38);
                for(int x = cnt - 1; x >= 1; --x) a[x] ^= a[x - 1];
                for(int i = 1; i < cnt; ++i) {
                    for(u64 temp = a[i]; temp; temp &= temp - 1) {
                        int b = lowestOneBit(temp);
                        val += mobValues[b / 8][b % 8] * cntParameter[min(i, 6)] / i;
                    }
                }
            }
            return val;
        };
        wv1 = calcWv(p1); wv2 = calcWv(p2);
        return wv1 - wv2;
    }
    
    double getEvaluateValue(int board[GRIDSIZE][GRIDSIZE], int myColor) {
        init(board, myColor);
        int blankCnts = max(2, min(55, popcnt(umap)));
        pretreatment();
        double a, b, c, d; computeParameter(a, b, c, d);
        double e = 0.1 * computeBlankValue();
        double v = a * evalW1[blankCnts] + b * evalW2[blankCnts] + c * evalW3[blankCnts] + d * evalW4[blankCnts] + e * evalW5[blankCnts] + evalW6[blankCnts];
        return (1.0 - (1.0 / (1.0 + exp(-v))) - 0.5) * 20000;
    }
};

int evaluateBoardAdvanced(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    auto startTime = chrono::high_resolution_clock::now();
    AdvancedEvaluator evaluator;
    int score = (int)evaluator.getEvaluateValue(board, myColor);
    auto endTime = chrono::high_resolution_clock::now();
    perfStats.evaluateBoardTime += chrono::duration<double, std::milli>(endTime - startTime).count();
    perfStats.evaluateBoardCount++;
    return score;
}

struct FastMoveScorer {
    u64 umap, myPieces, oppPieces, myReach, oppReach, contested;
    inline u64 expandOneStep(u64 pieces, u64 canGo) {
        u64 res = pieces;
        res |= (pieces >> 1 & DIR_WEST) | (pieces << 1 & DIR_EAST) | (pieces >> 8 & DIR_NORTH) | (pieces << 8 & DIR_SOUTH) | (pieces >> 9 & DIR_NW) | (pieces << 9 & DIR_SE) | (pieces >> 7 & DIR_NE) | (pieces << 7 & DIR_SW);
        return res & canGo;
    }
    inline u64 expandNSteps(u64 pieces, u64 canGo, int steps) {
        u64 res = pieces;
        for(int i = 0; i < steps; i++) res = expandOneStep(res, canGo);
        return res;
    }
    inline u64 slideOne(u64 a, int shift, u64 num, u64 dir, u64 can) {
        u64 res = a;
        for(int i = 0; i < 7; i++) {
            u64 next = (shift > 0 ? res << num : res >> num) & dir & can;
            if(next == 0) break;
            res |= next;
        }
        return res;
    }
    inline u64 queenReach(u64 pieces, u64 canGo) {
        u64 res = pieces, mask = canGo | pieces;
        res |= slideOne(pieces, -1, 1, DIR_WEST, mask); res |= slideOne(pieces, 1, 1, DIR_EAST, mask);
        res |= slideOne(pieces, -1, 8, DIR_NORTH, mask); res |= slideOne(pieces, 1, 8, DIR_SOUTH, mask);
        res |= slideOne(pieces, -1, 9, DIR_NW, mask); res |= slideOne(pieces, 1, 9, DIR_SE, mask);
        res |= slideOne(pieces, -1, 7, DIR_NE, mask); res |= slideOne(pieces, 1, 7, DIR_SW, mask);
        return res;
    }
    void init(int board[GRIDSIZE][GRIDSIZE], int myColor) {
        umap = myPieces = oppPieces = 0;
        for(int i = 0; i < GRIDSIZE; i++) {
            for(int j = 0; j < GRIDSIZE; j++) {
                int bit = i * 8 + j;
                if(board[i][j] == 0) umap |= (1ULL << bit);
                else if(board[i][j] == myColor) myPieces |= (1ULL << bit);
                else if(board[i][j] == -myColor) oppPieces |= (1ULL << bit);
            }
        }
        myReach = expandNSteps(myPieces, umap | myPieces, 2);
        oppReach = expandNSteps(oppPieces, umap | oppPieces, 2);
        contested = myReach & oppReach;
    }
    int scoreMove(int x0, int y0, int x1, int y1, int x2, int y2, int myColor) {
        int score = 0, bit0 = x0 * 8 + y0, bit1 = x1 * 8 + y1, bit2 = x2 * 8 + y2;
        u64 mask0 = 1ULL << bit0, mask1 = 1ULL << bit1, mask2 = 1ULL << bit2;
        if(contested & mask2) score += 15;
        u64 arrowNeighbor = expandOneStep(mask2, 0xFFFFFFFFFFFFFFFFULL);
        score += popcnt(arrowNeighbor & oppPieces) * 8;
        if((oppReach & mask2) && !(myReach & mask2)) score += 12;
        u64 newUmap = (umap & ~mask2) | mask0; newUmap &= ~mask1;
        score += popcnt(queenReach(mask1, newUmap));
        u64 newMyReach = expandNSteps((myPieces & ~mask0) | mask1, newUmap | ((myPieces & ~mask0) | mask1), 2);
        u64 newOppReach = expandNSteps(oppPieces, newUmap | oppPieces, 2);
        score += (popcnt(newMyReach & ~newOppReach) - popcnt(newOppReach & ~newMyReach)) * 2;
        score += (4 - abs(x1 - 3) - abs(y1 - 3));
        return score;
    }
    int scoreQueenMoveOnly(int x0, int y0, int x1, int y1) {
        int score = 0, bit0 = x0 * 8 + y0, bit1 = x1 * 8 + y1;
        u64 mask0 = 1ULL << bit0, mask1 = 1ULL << bit1;
        u64 newUmap = (umap | mask0) & ~mask1;
        score += popcnt(queenReach(mask1, newUmap)) * 3;
        u64 neighbor = expandOneStep(mask1, 0xFFFFFFFFFFFFFFFFULL) & ~mask1;
        score += popcnt(neighbor & newUmap) * 2;
        score += (4 - abs(x1 - 3) - abs(y1 - 3));
        if(contested & mask1) score += 5;
        if(expandNSteps(oppPieces, umap | oppPieces, 1) & mask1) score += 3;
        return score;
    }
} globalFastScorer;

struct Move {
    int x0, y0, x1, y1, x2, y2;
    Move(int x0_, int y0_, int x1_, int y1_, int x2_, int y2_) : x0(x0_), y0(y0_), x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
    Move() : x0(-1), y0(-1), x1(-1), y1(-1), x2(-1), y2(-1) {}
};

inline bool inMap(int x, int y) { return x >= 0 && x < GRIDSIZE && y >= 0 && y < GRIDSIZE; }

void copyBoard(int dst[GRIDSIZE][GRIDSIZE], int src[GRIDSIZE][GRIDSIZE]) { memcpy(dst, src, sizeof(int)*GRIDSIZE*GRIDSIZE); }

bool canReachOnBoard(int board[GRIDSIZE][GRIDSIZE], int x1, int y1, int x2, int y2) {
    if (x1 == x2 && y1 == y2) return true;
    int dx_dir = (x2 == x1) ? 0 : (x2 > x1 ? 1 : -1), dy_dir = (y2 == y1) ? 0 : (y2 > y1 ? 1 : -1);
    if ((dx_dir != 0 && dy_dir != 0 && abs(x2 - x1) != abs(y2 - y1)) || (dx_dir == 0 && dy_dir == 0)) return false;
    int steps = max(abs(x2 - x1), abs(y2 - y1));
    for (int i = 1; i < steps; ++i) if (board[x1 + dx_dir * i][y1 + dy_dir * i] != 0) return false;
    return true;
}

bool ProcStepOnBoard(int board[GRIDSIZE][GRIDSIZE], int x0, int y0, int x1, int y1, int x2, int y2, int color) {
    if (!inMap(x0, y0) || !inMap(x1, y1) || !inMap(x2, y2) || board[x0][y0] != color || board[x1][y1] != 0 || (board[x2][y2] != 0 && !(x2 == x0 && y2 == y0)) || !canReachOnBoard(board, x0, y0, x1, y1)) return false;
    board[x0][y0] = 0; board[x1][y1] = color; board[x2][y2] = OBSTACLE;
    return true;
}

int countMovesMade(int board[GRIDSIZE][GRIDSIZE]) {
    int cnt = 0;
    for(int i = 0; i < GRIDSIZE; i++) for(int j = 0; j < GRIDSIZE; j++) if(board[i][j] == OBSTACLE) cnt++;
    return cnt;
}

struct ColorScore { int connectivity, mobility, centerBonus; };
ColorScore calculateColorScore(int board[GRIDSIZE][GRIDSIZE], int color) {
    ColorScore res = {0, 0, 0};
    vector<pair<int,int>> queens;
    for(int i = 0; i < GRIDSIZE; i++) for(int j = 0; j < GRIDSIZE; j++) if(board[i][j] == color) { queens.push_back({i, j}); res.centerBonus += (8 - (abs(i - 4) + abs(j - 4))); }
    if(queens.empty()) return res;
    bool visited[GRIDSIZE][GRIDSIZE] = {false}; queue<pair<int,int>> q;
    for(auto& qp : queens) { visited[qp.first][qp.second] = true; q.push(qp); }
    while(!q.empty()) {
        auto [cx, cy] = q.front(); q.pop();
        for(int k = 0; k < 8; k++) {
            for(int d = 1; d < GRIDSIZE; d++) {
                int nx = cx + dx[k] * d, ny = cy + dy[k] * d;
                if(!inMap(nx, ny) || board[nx][ny] == OBSTACLE || board[nx][ny] == -color || board[nx][ny] == color) break;
                if(!visited[nx][ny]) { visited[nx][ny] = true; res.connectivity++; q.push({nx, ny}); }
                res.mobility++;
            }
        }
    }
    return res;
}

int evaluateBoard(int board[GRIDSIZE][GRIDSIZE], int myColor) {
    auto startTime = chrono::high_resolution_clock::now();
    bool mCM = false, oCM = false;
    for(int i = 0; i < GRIDSIZE && (!mCM || !oCM); i++) for(int j = 0; j < GRIDSIZE && (!mCM || !oCM); j++) if(board[i][j] == myColor || board[i][j] == -myColor) {
        for(int k = 0; k < 8; k++) { int nx = i + dx[k], ny = j + dy[k]; if(inMap(nx, ny) && board[nx][ny] == 0) { if(board[i][j] == myColor) mCM = true; else oCM = true; break; } }
    }
    if(!mCM && !oCM) { ColorScore ms = calculateColorScore(board, myColor), os = calculateColorScore(board, -myColor); return ms.connectivity > os.connectivity ? INF - 1 : (ms.connectivity < os.connectivity ? -INF + 1 : 0); }
    if(!mCM) return -INF + 1; if(!oCM) return INF - 1;
    ColorScore ms = calculateColorScore(board, myColor), os = calculateColorScore(board, -myColor);
    int score = (ms.connectivity - os.connectivity) * WEIGHT_CONNECTIVITY + (ms.mobility - os.mobility) * WEIGHT_MOBILITY + (ms.centerBonus - os.centerBonus) * WEIGHT_CENTER;
    perfStats.evaluateBoardTime += chrono::duration<double, std::milli>(chrono::high_resolution_clock::now() - startTime).count();
    perfStats.evaluateBoardCount++;
    return score;
}

struct TTEntry { long long hash; int depth, score, flag; Move bestMove; TTEntry() : hash(0), depth(-1), score(0), flag(0), bestMove() {} };
const int TT_SIZE = 1048576;
TTEntry transpositionTable[TT_SIZE];
long long zobristTable[GRIDSIZE][GRIDSIZE][3], zobristBlackTurn;
void initZobrist() { srand(20251225); for(int i=0; i<GRIDSIZE; i++) for(int j=0; j<GRIDSIZE; j++) for(int k=0; k<3; k++) zobristTable[i][j][k] = ((long long)rand() << 32) | rand(); zobristBlackTurn = ((long long)rand() << 32) | rand(); }
long long computeHash(int board[GRIDSIZE][GRIDSIZE], int currentPlayer) {
    long long h = 0;
    for(int i=0; i<GRIDSIZE; i++) for(int j=0; j<GRIDSIZE; j++) { if(board[i][j] == grid_black) h ^= zobristTable[i][j][0]; else if(board[i][j] == grid_white) h ^= zobristTable[i][j][1]; else if(board[i][j] == OBSTACLE) h ^= zobristTable[i][j][2]; }
    if(currentPlayer == grid_black) h ^= zobristBlackTurn;
    return h;
}

    
struct QueenMove { int x0, y0, x1, y1, score; QueenMove(int x0_, int y0_, int x1_, int y1_, int s = 0) : x0(x0_), y0(y0_), x1(x1_), y1(y1_), score(s) {};QueenMove() : x0(-1), y0(-1), x1(-1), y1(-1), score(0) {} };
vector<QueenMove> generateTopQueenMoves(int board[GRIDSIZE][GRIDSIZE], int color, int k1) {
    vector<QueenMove> qms; globalFastScorer.init(board, color);
    for(int i = 0; i < GRIDSIZE; i++) for(int j = 0; j < GRIDSIZE; j++) if(board[i][j] == color) {
        for(int k = 0; k < 8; k++) for(int d = 1; d < GRIDSIZE; d++) { int x1 = i + dx[k] * d, y1 = j + dy[k] * d; if(!inMap(x1, y1) || board[x1][y1] != 0) break; qms.push_back(QueenMove(i, j, x1, y1, globalFastScorer.scoreQueenMoveOnly(i, j, x1, y1))); }
    }
    if((int)qms.size() > k1) { nth_element(qms.begin(), qms.begin() + k1, qms.end(), [](const QueenMove& a, const QueenMove& b) { return a.score > b.score; }); qms.resize(k1); }
    sort(qms.begin(), qms.end(), [](const QueenMove& a, const QueenMove& b) { return a.score > b.score; });
    return qms;
}

vector<Move> generateTopFullMoves(int board[GRIDSIZE][GRIDSIZE], const vector<QueenMove>& qms, int color, int k2) {
    vector<pair<int, Move>> sms; sms.reserve(qms.size() * 20); globalFastScorer.init(board, color);
    for(const auto& qm : qms) {
        auto cellAt = [&](int x, int y) { if(x == qm.x0 && y == qm.y0) return 0; if(x == qm.x1 && y == qm.y1) return color; return board[x][y]; };
        for(int k = 0; k < 8; k++) for(int d = 1; d < GRIDSIZE; d++) { int x2 = qm.x1 + dx[k] * d, y2 = qm.y1 + dy[k] * d; if(!inMap(x2, y2) || cellAt(x2, y2) != 0) break; sms.push_back({globalFastScorer.scoreMove(qm.x0, qm.y0, qm.x1, qm.y1, x2, y2, color), Move(qm.x0, qm.y0, qm.x1, qm.y1, x2, y2)}); }
    }
    if((int)sms.size() > k2) { nth_element(sms.begin(), sms.begin() + k2, sms.end(), [](const pair<int,Move>& a, const pair<int,Move>& b) { return a.first > b.first; }); sms.resize(k2); }
    sort(sms.begin(), sms.end(), [](const pair<int,Move>& a, const pair<int,Move>& b) { return a.first > b.first; });
    vector<Move> res; for(auto& it : sms) res.push_back(it.second);
    return res;
}

vector<Move> generateMovesWithTwoPhase(int board[GRIDSIZE][GRIDSIZE], int color) {
    auto s1 = chrono::high_resolution_clock::now();
    vector<QueenMove> qms = generateTopQueenMoves(board, color, K1_QUEEN_MOVES);
    if(qms.empty()) return vector<Move>();
    auto s2 = chrono::high_resolution_clock::now();
    perfStats.totalphase1GenMoveTime += chrono::duration<double, std::milli>(s2 - s1).count(); perfStats.phase1Count++;
    vector<Move> fms = generateTopFullMoves(board, qms, color, K2_FULL_MOVES);
    perfStats.totalphase2GenMoveTime += chrono::duration<double, std::milli>(chrono::high_resolution_clock::now() - s2).count(); perfStats.phase2Count++;
    return fms;
}

int alphaBetaSearch(int board[GRIDSIZE][GRIDSIZE], int depth, int alpha, int beta, int color, Move& bestMove, long long hash) {
    int ttIdx = (hash & 0x7FFFFFFFFFFFFFFF) % TT_SIZE; TTEntry* tte = &transpositionTable[ttIdx];
    if(tte->hash == hash && tte->depth >= depth) {
        if(tte->flag == 0) { bestMove = tte->bestMove; return tte->score; }
        else if(tte->flag == 1) alpha = max(alpha, tte->score);
        else if(tte->flag == 2) beta = min(beta, tte->score);
        if(alpha >= beta) { bestMove = tte->bestMove; return tte->score; }
    }
    const int aO = alpha, bO = beta;
    vector<Move> moves = generateMovesWithTwoPhase(board, color);
    if(depth == 0 || moves.empty()) return (EVAL_MODE == 1) ? evaluateBoardAdvanced(board, color) : evaluateBoard(board, color);
    int bS = -INF; Move lBM;
    for(const Move& m : moves) {
        int nB[GRIDSIZE][GRIDSIZE]; copyBoard(nB, board);
        if(!ProcStepOnBoard(nB, m.x0, m.y0, m.x1, m.y1, m.x2, m.y2, color)) continue;
        long long nH = hash;
        if(board[m.x0][m.y0] == grid_black) nH ^= zobristTable[m.x0][m.y0][0]; else nH ^= zobristTable[m.x0][m.y0][1];
        if(color == grid_black) nH ^= zobristTable[m.x1][m.y1][0]; else nH ^= zobristTable[m.x1][m.y1][1];
        nH ^= zobristTable[m.x2][m.y2][2] ^ zobristBlackTurn;
        Move oBM; int s = -alphaBetaSearch(nB, depth - 1, -beta, -alpha, -color, oBM, nH);
        if(s > bS) { bS = s; lBM = m; }
        alpha = max(alpha, s); if(alpha >= beta) break;
    }
    tte->hash = hash; tte->depth = depth; tte->score = bS; tte->bestMove = lBM;
    tte->flag = (bS <= aO) ? 2 : (bS >= bO ? 1 : 0);
    bestMove = lBM; return bS;
}

Move iterativeDeepeningSearch(int color, int maxDepth) {
    auto s = chrono::high_resolution_clock::now(); perfStats.reset();
    Move bM; long long h = computeHash(gridInfo, color);
    for(int d = 1; d <= maxDepth; d++) { Move cBM; alphaBetaSearch(gridInfo, d, -INF, INF, color, cBM, h); if(cBM.x0 >= 0) bM = cBM; }
    perfStats.totalTime = chrono::duration<double, std::milli>(chrono::high_resolution_clock::now() - s).count();
    perfStats.printStats(); return bM;
}

int main() {
    initZobrist();
    int x0, y0, x1, y1, x2, y2, turnID;
    gridInfo[0][2] = gridInfo[2][0] = gridInfo[5][0] = gridInfo[7][2] = grid_black;
    gridInfo[0][5] = gridInfo[2][7] = gridInfo[5][7] = gridInfo[7][5] = grid_white;
    if(!(cin >> turnID)) return 0;
    currBotColor = grid_white;
    for (int i = 0; i < turnID; i++) {
        cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2;
        if (x0 == -1) currBotColor = grid_black;
        else ProcStepOnBoard(gridInfo, x0, y0, x1, y1, x2, y2, -currBotColor);
        if (i < turnID - 1) { cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2; if (x0 >= 0) ProcStepOnBoard(gridInfo, x0, y0, x1, y1, x2, y2, currBotColor); }
    }
    int mC = countMovesMade(gridInfo);
    int cD = (mC >= MIDGAME_THRESHOLD) ? ENDGAME_DEPTH : OPENING_DEPTH;
    srand(time(0));
    Move bM = iterativeDeepeningSearch(currBotColor, cD);
    if (bM.x0 >= 0) cout << bM.x0 << ' ' << bM.y0 << ' ' << bM.x1 << ' ' << bM.y1 << ' ' << bM.x2 << ' ' << bM.y2 << endl;
    else cout << "-1 -1 -1 -1 -1 -1" << endl;
    return 0;
}
