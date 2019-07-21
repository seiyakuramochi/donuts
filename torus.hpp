
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifndef ROUTING_TORUS_H
#define ROUTING_TORUS_H

typedef struct NODE
{
    int index;
    int *value; // n-dimension vector

    double **P;
    /*
     * P: matrix of probabilities
     * P contain probability of node depends on h and d
     * where 1 <= h <= n, 1 <= d <= diameter
    */

    int *neighbors;  // indices
} Node;


using namespace boost::numeric::ublas;
class Torus {

public:

    // ↓ ====== called in constructor
    void allocate();
    void setValues();
    void createConnection();
    void saveNeighbor();
    // ↑ ====== called in constructor

    /*
     * isNeighborとhasLinkは全く同じ挙動になるが、
     * hasLinkはisNeighborの結果を利用する(=adjacencyを参照する)ので速い
     *
     * Function isNeighbor and hasLink has same behavior
     * hasLink is faster than isNeighbor (hasLink uses result of isNeighbor)
     */
    bool isNeighbor(Node *a, Node *b);

    // ↓ ====== just to read/write data
    bool hasLink(Node *a, Node *b);
    void setFaultyLink(Node *a, Node *b);
    bool hasFaultyLink(Node *a, Node *b);
    void setProbability(Node *a, int h, int d, double p);
    double getProbability(Node *a, int h, int d);
    // ↑ ====== just to read/write data

    // 各種距離の計算
    // Functions to calculate distance
    int distance(Node* a, Node* b);
    int hammingDistance(Node* a, Node* b);

    // 近傍ノードのクラス分類
    // Functions to classify neighbor nodes
    bool inPre(Node* neighbor, Node* a, Node* b);
    bool inSpr(Node* neighbor, Node* a, Node* b);

    int k;
    int n;
    Node *nodes; // nodes
    int V; // number of nodes
    int E; // number of edges
    int diameter;
    //int **adjacency; // graph adjacency matrix
    //int **F; // graph faulty links matrix

    // スパース行列
    mapped_matrix<int> adjacency;
    mapped_matrix<int> F;

    // n dimensional vector (string) -> index of a node (in nodes)
    std::unordered_map<std::string, int> node_value_index;

    Torus(int n, int k); // constructor
    ~Torus(); // destructor

    // for debugging
    void printAdjMatrix();
    void printFaultyLinks();
    void printProbabilities();

    // ↓割りと重要な関数たち Very important algorithms
    void setRandomFaultyLinks(double p_faulty);

    std::string get_e_plus(Node *a, int i);
    std::string get_e_minus(Node *a, int i);

    // algorithm P
    void calcRoutingProbabilities();
    double calc_p_i(Node *a, Node *a_plus_e, int h, int d);

    void calcPQ(Node *a, int h, int d, double *p, double *q);
    void XxYxMainProcess(int h, int d, double* p, double* q,
        int out_xp[], int out_yp[]);
    void XpYp(Node *a, int h, int d, double in_p[], double in_q[],
         double out_p[], int out_xp[], int out_yp[]);
    void XqYq(Node *a, int h, int d, double in_p[], double in_q[],
         double out_p[], int out_xp[], int out_yp[]);

    // routing algorithms
    int route(Node *prev, Node* c, Node* t, std::unordered_map<int, bool> visited, int d);
    int brute(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d);

    // breath-first search for optimal routing
    int bfs(Node *c, Node *t, std::unordered_map<int, bool> visited);
};

void sortCouple(double *a, double *b, int len);

#endif //ROUTING_TORUS_H
