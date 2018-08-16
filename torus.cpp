#include <iostream>
#include <cmath>
#include <cassert>
#include <unordered_map>
#include <queue>
#include <random>
#include <string.h>
#include <stdlib.h>
 
#include "bbst.hpp"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "torus.h"

#define DELIVERY_FAIL -1;

using namespace boost::numeric::ublas;

Torus::Torus(int _n, int _k) {
    assert(_k > 2);
    assert(_n >= 0);

    // 各種パラメータを計算
    n = _n;
    k = _k;
    V = static_cast<int>(pow(k, n));
    E = n*V;
    diameter = static_cast<int>(floor(k/2)*n);

    // 隣接行列と故障リンク行列をV*Vのスパース行列として初期化
    adjacency = mapped_matrix<int>(V, V);
    F = mapped_matrix<int>(V, V);

    // メモリを確保する
    allocate();
    // nodes->valueをセットする
    setValues();
    // トーラスの接続行列を計算
    createConnection();
}


void Torus::allocate() {
    // O(n k^n)

    nodes = (Node*)malloc(V*sizeof(Node));

    for(int i=0; i<V; i++){
        nodes[i].index = i;
        nodes[i].value = (int*)malloc(n*sizeof(int));

        nodes[i].P = (float**)malloc(n*sizeof(float*));
        for(int j=0; j<n; j++) {
            nodes[i].P[j] = (float *) malloc(diameter * sizeof(float));
            std::fill_n(nodes[i].P[j], diameter, -1);
        }

        nodes[i].neighbors = (int*)malloc(2*n*sizeof(int));
    }
}


void Torus::setValues() {
    // O(n k^n)

    int m, po;
    using namespace std;
    for(int i=0; i<V; i++) {

        m = i;
        // for all dimension
        for(int j=n-1; j>0; j--) {
            po = int(pow(k, j));
            (nodes)[i].value[n-j-1] = int(m / po);
            m %= po;
        }

        (nodes)[i].value[n-1] = m;

        std::string v = "";
        for(int j=0; j<n; j++)
            v += std::to_string(nodes[i].value[j]) + ",";
        // cout << v << "=" << i << endl;
        node_value_index[v] = i;
    }
}


void Torus::createConnection() {
    // O(n^2 k^n)
    using namespace std;
    for(int i=0; i<V; i++) {
        int i_neighbor = 0;
        for (int j = 0; j < n; j++) {
            std::string v_plus = "", v_minus = "";
            for (int l = 0; l < n; l++) {
                if (l == j) {
                    v_plus += to_string((nodes[i].value[l] + 1) % k) + ",";
                    if(nodes[i].value[l] == 0)
                        v_minus += to_string(k-1) + ",";
                    else
                        v_minus += to_string(nodes[i].value[l] - 1) + ",";
                }else {
                    v_plus += to_string(nodes[i].value[l]) + ",";
                    v_minus += to_string(nodes[i].value[l]) + ",";
                }
            }
            //cout << v_plus;
            //cout << " || ";

            //cout << v_minus;
            //cout << endl;
            //cout << i << ", " << j << ":" << node_value_index.at(v_plus) << endl;
            //cout << nodes[i].value[0] << "," << nodes[i].value[1] <<"→" <<  v_plus << "=" << i <<","<<node_value_index.at(v_plus) << endl;
            adjacency(node_value_index[v_plus], i) = 1;
            nodes[i].neighbors[i_neighbor++] = node_value_index[v_plus];
            nodes[i].neighbors[i_neighbor++] = node_value_index[v_minus];
        }
    }
}


bool Torus::isNeighbor(Node *a, Node *b) {
    // 異なり次元数, 異なっているインデックス
    int hamming_dist = 0, i_different = 0;
    // 異なっているインデックスの値
    int a_differ, b_differ;

    // for all dimension
    for (int i = 0; i < n; i++)
        if(a->value[i] != b->value[i]) {
            hamming_dist++;
            i_different = i;
        }

    if(hamming_dist != 1)
        return false;

    a_differ = a->value[i_different];
    b_differ = b->value[i_different];

    // 1違うパターン
    if(abs(a_differ - b_differ) == 1)
        return true;

    // 0とk-1のパターン
    if((a_differ == k-1 and b_differ == 0) or (a_differ == 0 and b_differ == k-1))
        return true;

    return false;
}


void Torus::printFaultyLinks() {
    for(int i=0; i<V; i++) {
        for (int j=0; j<V; j++)
            if(F(i, j) == 1){
                std::cout << "((";
                for(int l=0; l<n; l++) {
                    std::cout << nodes[i].value[l];
                    if(l < n-1)
                        std::cout << ", ";
                }
                //std::cout << std::endl;
                std::cout << "), (";
                for(int l=0; l<n; l++){
                    std::cout << nodes[j].value[l];
                    if(l < n-1)
                        std::cout << ", ";
                }
                std::cout << "))" << std::endl;

            };
    }
}


void Torus::setRandomFaultyLinks(float p_faulty) {
    // O(k^n)

    int n_faulty = static_cast<int>(p_faulty * E);
    int count_faulty=0;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dice_V(0, V-1);
    std::uniform_int_distribution<int> dice_neighbor(0, 2*n-1);

    assert(0.0 <= p_faulty and p_faulty <= 1.0);
    assert(0 <= n_faulty and n_faulty <= E);

    while(count_faulty < n_faulty) {
        Node *a = &(nodes[dice_V(mt)]);
        Node *b = &(nodes[a->neighbors[dice_neighbor(mt)]]);
        assert(hasLink(a, b));

        if (not hasFaultyLink(a, b)) {
            setFaultyLink(a, b);
            count_faulty++;
        }
    }
}


float Torus::calc_p_i(Node *a, Node *a_plus_e, int h, int d) {
    float p = 0.0;
    if (h == 1){
        if (not hasFaultyLink(a, a_plus_e))
            p += getProbability(a_plus_e, h, d - 1);
    }else if (h == d){
        if (not hasFaultyLink(a, a_plus_e))
            p += getProbability(a_plus_e, h - 1, d - 1);
    }else{
        if (not hasFaultyLink(a, a_plus_e)){
            p += (h - 1) * getProbability(a_plus_e, h - 1, d - 1) +
                 (d - h) * getProbability(a_plus_e, h, d - 1);
            p /= (d - 1);
        }
    }
    return p;
}


void sortCouple(float *a, float *b, int len){
    // aの基準で、aとbをソートする
    using namespace std;

    float *tmp_a = (float*)malloc(len*sizeof(float));
    float *tmp_b = (float*)malloc(len*sizeof(float));

    // まずインデックスの列をソートする
    std::vector<int> argsort;
    argsort.resize(len);
    iota(argsort.begin(), argsort.end(), 0);
    std::sort(
        std::begin(argsort), 
        std::end(argsort),
        [a](int i, int j) { return a[i] < a[j]; }
    );

    // インデックスのソートに基づいてtmpの値をソートする
    for(int i=0; i<len; i++){
        tmp_a[i] = a[argsort[i]];
        tmp_b[i] = b[argsort[i]];
    }

    // コピーして終わり
    memcpy(a, tmp_a, sizeof(float)*len);
    memcpy(b, tmp_b, sizeof(float)*len);

    free(tmp_a);
    free(tmp_b);
}


std::string Torus::get_e_plus(Node *a, int i){
    using namespace std;
    string e_plus = "";
    for (int l = 0; l < n; l++)
        if (l == i)
            e_plus += to_string((a->value[l] + 1) % k) + ",";
        else
            e_plus += to_string(a->value[l]) + ",";
    return e_plus;
}


std::string Torus::get_e_minus(Node *a, int i){
    using namespace std;
    string e_minus = "";
    for (int l = 0; l < n; l++)
        if (l == i) {
            if (a->value[l] == 0)
                e_minus += to_string(k-1) + ",";
            else
                e_minus += to_string(a->value[l] - 1) + ",";
        } else {
            e_minus += to_string(a->value[l]) + ",";
        }
    return e_minus;
}


void Torus::XxYxPreProcess(Node *a, int h, int d, float *p, float *q){
    std::string e_plus, e_minus;

    // pとqを計算する
    for (int i=0; i<n; i++) {
        e_plus = get_e_plus(a, i);
        e_minus = get_e_minus(a, i);
 
        assert(nodes[node_value_index[e_plus]].value[0] == e_plus[0]-48);
        assert(nodes[node_value_index[e_minus]].value[0] == e_minus[0]-48);
        
        assert(hasLink(a, &nodes[node_value_index[e_plus]]));
        assert(hasLink(a, &nodes[node_value_index[e_minus]]));

        p[i] = calc_p_i(a, &nodes[node_value_index[e_plus]], h, d);
        q[i] = calc_p_i(a, &nodes[node_value_index[e_minus]], h, d);
    }
}
 

void Torus::XxYxMainProcess(int h, int d, float* p, float* q,
        float out_p[], int out_xp[], int out_yp[]){

    // O(nlogn)

    using namespace std;

    struct TNode* t_l = NULL;
    struct TNode* t_r = NULL;
 
    for(int i=0; i<n; i++)
        t_r = insert(t_r, q[i]);

    assert(CountLesser(t_r, 1000.0) == n);
    assert(CountLesser(t_r, -1.0) == 0);
   
    for(int i=1; i<=n; i++){
        preOrder(t_l);
        out_xp[i-1] = CountLesser(t_l, p[i-1]);
        cout << "threshold:" << p[i-1] << endl;
        cout << "out:" << out_xp[i-1] << endl << endl;
        t_r = deleteNode(t_r, q[i-1]);

        preOrder(t_r);
        out_yp[i-1] = (i-1-out_xp[i-1]) + CountLesser(t_r, p[i-1]);
        cout << "threshold:" << p[i-1] << endl;
        cout << "out:" << (i-1-out_xp[i-1]) << "+ " << CountLesser(t_r, p[i-1]) << "=" << out_yp[i-1] << endl << endl;
        t_l = insert(t_l, q[i-1]);
    }

    memcpy(out_p, p, n * sizeof(float));
    free(p);
    free(q);
    delete_tree(t_r);
    delete_tree(t_l);
}


void Torus::XpYp(Node *a, int h, int d,
        float out_p[], int out_xp[], int out_yp[]) {

    // pとqの領域を確保する
    float *p = (float*)malloc(n*sizeof(float));
    float *q = (float*)malloc(n*sizeof(float));

    // pとqの値を計算する
    XxYxPreProcess(a, h, d, p, q);

    // ここでp[i]とq[i]は次元順序になっている
    // ので、pの大きさベースで確率の順序に並び替える(破壊的変更)
    sortCouple(p, q, n);

    // xpとypを計算する
    XxYxMainProcess(h, d, p, q, out_p, out_xp, out_yp);
}


void Torus::XqYq(Node *a, int h, int d,
        float out_q[], int out_xq[], int out_yq[]) {

    // pとqの領域を確保する
    float *p = (float*)malloc(n*sizeof(float));
    float *q = (float*)malloc(n*sizeof(float));

    // pとqの値を計算する
    XxYxPreProcess(a, h, d, p, q);

    // ここでp[i]とq[i]は次元順序になっている
    // ので、qの大きさベースで確率の順序に並び替える(破壊的変更)
    sortCouple(q, p, n);

    // xqとyqを計算する
    XxYxMainProcess(h, d, q, p, out_q, out_xq, out_yq);
}


int combination(int n, int r){
    //using namespace std;
    //cout << "n,r,ncr="<< n << ","<<r<<","<<tgamma(n+1)/(tgamma(n-r+1)*tgamma(r+1))<<endl;
    if(r > n)
        return 0;
    return tgamma(n+1)/(tgamma(n-r+1)*tgamma(r+1));
}


// procedure P(a)
void Torus::calcRoutingProbabilities() {
    float *ps, *qs;
    int *xp, *yp, *xq, *yq;
    Node *a;
    float p, p11;
    Node *neighbor;

    ps = (float *)malloc(n*sizeof(float));
    qs = (float *)malloc(n*sizeof(float));
    xp = (int *)malloc(n*sizeof(int));
    yp = (int *)malloc(n*sizeof(int));
    xq = (int *)malloc(n*sizeof(int));
    yq = (int *)malloc(n*sizeof(int));

    /*
     * まず全てのノードについてP(a)_{1,1}を計算
     * First, calculate P(a)_{1,1} for all nodes
     */
//#pragma omp parallel for
    for(int i=0; i<V; i++) {
        a = &nodes[i];

        p11 = 0.0;
        for (int j=0; j<2*n; j++) {
            neighbor = &nodes[nodes[i].neighbors[j]];
            if(not hasFaultyLink(a, neighbor))
                p11 += 1.0;
        }
        p11 /= 2*n;
        //std::cout << "p11=" << p11 << std::endl;
        setProbability(a, 1, 1, p11);
    }

    /*
     * P(a)_{h,d}を計算していく
     * P(a)_{h,d}の計算には, P(a)_{h,d-1}やP(a)_{h-1,d-1}が必要になる.
     * P(a)_{1,1}だけで計算できるP(a)_{1,2}から順番に計算していく(動的計画法)
     *
     * Calculation of P(a)_{h,d}
     * P(a)_{h,d} needs P(a)_{h,d-1} and P(a)_{h-1,d-1}.
     * Calculation starts with P(a)_{1,2}.
     */
    for (int d = 2; d <= diameter; d++)
        for (int h = 1; h <= std::min(n, d); h++) {
            // for all nodes
            for (int i = 0; i < V; i++) {
                a = &(nodes[i]);
                
                XpYp(a, h, d, ps, xp, yp);
                XqYq(a, h, d, qs, xq, yq);

                p = 0.0;
                for(int ii=1; ii<=n; ii++){
                    for(int kk=0; kk<=h-1; kk++){
                        // xp yp xq yq は 0-origin なので ii から 1 引く
                        p += pow(2, kk) * ( 
                             combination(xp[ii-1], kk) * 
                             combination(yp[ii-1], h-kk-1) *
                             ps[ii-1] +
                             combination(xq[ii-1], kk) *
                             combination(yq[ii-1], h-kk-1) *
                             qs[ii-1]);

                        std::cout << "ps[i]=" << ps[ii-1] << std::endl;
                        std::cout << "qs[i]=" << qs[ii-1] << std::endl;
                        std::cout << "xp[i]=" << xp[ii-1] << std::endl;
                        std::cout << "k=" << kk << std::endl;

                        std::cout << "yp[i]=" << yp[ii-1] << std::endl;
                        std::cout << "h-k1=" << h-kk-1 << std::endl;
                        std::cout << "xq[i]=" << xq[ii-1] << std::endl;
                        std::cout << "yq[i]=" << yq[ii-1] << std::endl;
                        
                        std::cout << "c*c=" << combination(xp[ii-1], kk) * 
                             combination(yp[ii-1], h-kk-1)<< std::endl;
                        std::cout << "p*c*c=" <<combination(xp[ii-1], kk) * 
                             combination(yp[ii-1], h-kk-1) *
                             ps[ii-1]<< std::endl << std::endl;
                    }
                }
                p /= pow(2, h) * combination(n, h);
                
                std::cout << "setp=" << p << std::endl;
                setProbability(a, h, d, p);
                std::cout << "h,d=" << h << "," << d << std::endl;

                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
}   


// 1-α(a, b)
bool Torus::hasFaultyLink(Node *a, Node *b) {
    // 対称行列で上三角しか使ってないので入れ替えて補う
    bool ab = static_cast<bool>(F(a -> index, b -> index) != 0);
    bool ba = static_cast<bool>(F(b -> index, a -> index) != 0);

    assert(not (ab and ba));

    return (ab or ba);
}


void Torus::setFaultyLink(Node *a, Node *b) {
    F(a->index, b->index) = 1;
}


bool Torus::hasLink(Node *a, Node *b) {
    bool ab = static_cast<bool>(adjacency(a -> index, b -> index) != 0);
    bool ba = static_cast<bool>(adjacency(b -> index, a -> index) != 0);

    assert(not (ab and ba));

    return (ab or ba);
}


void Torus::setProbability(Node *a, int h, int d, float p) {
    assert(0 < h and h <= n);
    assert(0 < d and d <= diameter);
    assert(0.0 <= p and p <= 1.0);

    // setProbability()によってセットされるのが1回目であることをチェック
    assert(nodes[a->index].P[h-1][d-1] == -1);

    nodes[a->index].P[h-1][d-1] = p;
}


float Torus::getProbability(Node *a, int h, int d) {
    assert(0 <= h and h <= n);
    assert(0 <= d and d <= floor(k/2)*n);

    if(h == 0) {
        assert(d == 0);
        return 1.0;
    }

    // setProbability()によってセットされた値であることをチェック
    assert(nodes[a->index].P[h-1][d-1] != -1);

    return nodes[a->index].P[h-1][d-1];
}


void Torus::printProbabilities() {
    Node* a;
    std::cout << "Node ";
    for (int d = 1; d <= floor(k / 2) * n; d++)
        for (int h = 1; h <= std::min(n, d); h++){
            std::cout << "P_{" << h << "," << d << "} ";
        }
    std::cout << std::endl;

    for(int i=0; i<V; i++) {
        a = &nodes[i];
        std::cout << "(";
        for(int j=0; j<n; j++) {
            std::cout << a->value[j];
            if(j < n-1)
                std::cout << ",";
        }
        std::cout << ") ";

        for (int d = 1; d <= floor(k / 2) * n; d++)
            for (int h = 1; h <= std::min(n, d); h++) {
                std::cout << int(getProbability(a, h, d)*100)/100.0
                          << " ";
            }
        std::cout <<  std::endl;
    }
}


int Torus::distance(Node *a, Node *b) {
    int d = 0;
    for(int i=0; i<n; i++)
        if(abs(a->value[i]-b->value[i]) > floor(k/2))
            d += k-abs(a->value[i]-b->value[i]);
        else
            d += abs(a->value[i]-b->value[i]);

    assert(0 <= d and d <= diameter);
    return d;
}


int Torus::hammingDistance(Node *a, Node *b) {
    int d = 0;
    for(int i=0; i<n; i++)
        if(a->value[i] != b->value[i])
            d++;
    return d;
}


bool Torus::inPre(Node *neighbor, Node *a, Node *b) {
    return distance(neighbor, b) < distance(a, b);
}


bool Torus::inSpr(Node *neighbor, Node *a, Node *b) {
    return not inPre(neighbor, a, b);
}


int Torus::route(Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    int dist_ct, h;
    float p, p_max;
    Node *neighbor, *max_neighbor;

    //std::cout << "c:";
    //for(int i=0; i<n; i++)
    //    std::cout << c->value[i] << " ";
    //std::cout << std::endl;

    // Base case
    if(c->index == t->index)
        return d;

    dist_ct = distance(c, t);
    p_max = -1;
    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(hasLink(c, neighbor));
        if (inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor) ) {
            h = hammingDistance(neighbor, t);
            p = getProbability(neighbor, h, dist_ct - 1);

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
        return route(max_neighbor, t, visited, d + 1);
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(hasLink(c, neighbor));
        if (inSpr(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            h = hammingDistance(neighbor, t);
            p = getProbability(neighbor, h, std::min(dist_ct + 1, diameter));

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
        return route(max_neighbor, t, visited, d + 1);
    }

    return DELIVERY_FAIL;
}


int Torus::brute(Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    Node *neighbor;

    //std::cout << "c:";
    //for(int i=0; i<n; i++)
    //    std::cout << c->value[i] << " ";
    //std::cout << std::endl;

    if(c->index == t->index)
        return d;

    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        if(inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            if(visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
            return brute(neighbor, t, visited, d + 1);
        }
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        if (not hasFaultyLink(c, neighbor)){
            assert(inSpr(neighbor, c, t));
            if (visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
            return brute(neighbor, t, visited, d + 1);
        }
    }

    return DELIVERY_FAIL;
}


int Torus::bfs(Node *c, Node *t, std::unordered_map<int, bool> visited) {
    std::queue<int> q_i, q_d;
    int a_d;
    Node *a, *neighbor;

    visited[c->index] = true;
    q_i.push(c->index);
    q_d.push(0);

    while(not q_i.empty()){

        assert(q_d.size() == q_i.size());

        a = &(nodes[q_i.front()]);
        q_i.pop();
        a_d = q_d.front();
        q_d.pop();

        if(a->index == t->index)
            return a_d;

        for(int i=0; i<2*n; i++) {
            neighbor = &(nodes[a->neighbors[i]]);
            if(not hasFaultyLink(a, neighbor) and visited.count(neighbor->index) == 0){
                visited[neighbor->index] = true;
                q_i.push(neighbor->index);
                q_d.push(a_d+1);
            }
        }
    }

    assert(q_d.empty());
    return DELIVERY_FAIL;
}

Torus::~Torus() {
    adjacency.clear();
    F.clear();

    for(int i=0; i<V; i++){
        free(nodes[i].value);
        for(int j=0; j<n; j++)
            free(nodes[i].P[j]);
        free(nodes[i].P);
        free(nodes[i].neighbors);
    }
    free(nodes);
}
