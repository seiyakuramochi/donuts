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

#include "torus.hpp"

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
        nodes[i].value = new int[n];

        nodes[i].P = new double[diameter];
        std::fill_n(nodes[i].P, diameter, -1);

        nodes[i].neighbors = new int[2*n];
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


void Torus::setRandomFaultyNodes(double p_faulty) {
    // O(k^n)

    int n_faulty = static_cast<int>(p_faulty * V);
    int count_faulty=0;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dice_V(0, V-1);

    bool* is_faulty_node = new bool[V];
    for(int i=0; i<V; i++)
        is_faulty_node[i] = false;

    assert(0.0 <= p_faulty and p_faulty <= 1.0);
    assert(0 <= n_faulty and n_faulty <= V);

    while(count_faulty < n_faulty) {
        Node *a = &(nodes[dice_V(mt)]);
        if(is_faulty_node[a->index])
            continue;

        for (int j=0; j<2*n; j++) {
            Node* neighbor = &nodes[a->neighbors[j]];
            if (hasFaultyLink(a, neighbor))
                continue;
            setFaultyLink(a, neighbor);
            assert(hasFaultyLink(a, neighbor));
        }
        is_faulty_node[a->index] = true;
        count_faulty++;
    }
}


void Torus::setRandomFaultyLinks(double p_faulty) {
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
            assert(hasFaultyLink(a, b));
            count_faulty++;
        }
    }
}


// procedure P(a)
void Torus::calcRoutingProbabilities() {
    Node *a;
    double p, p1;
    Node *neighbor;

    for(int i=0; i<V; i++) {
        a = &nodes[i];

        p1 = 0;
        for (int j=0; j<2*n; j++) {
            neighbor = &nodes[nodes[i].neighbors[j]];
            if(hasFaultyLink(a, neighbor))
                p1 += 1.0;
        }
        p1 /= 2*n;
        setProbability(a, 1, p1);
    }

    assert(n == 3);

    for (int d = 2; d <= diameter; d++)
        // for each nodes in Torus
        for (int i = 0; i < V; i++) {
            a = &(nodes[i]);

            p = 1.0;
            // for each neighbors of node a
            for (int j=0; j<2*n; j++) {
                double r = 1.0;
                neighbor = &nodes[nodes[i].neighbors[j]];
                if(not hasFaultyLink(a, neighbor))
                    for(int h=1; h<=std::min(d, n); h++)
                        r -= h/6 * (1-getProbability(a, d-1));
                p *= r;
            }
            setProbability(a, d, p);
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


void Torus::setProbability(Node *a, int d, double p) {
    assert(0 < d and d <= diameter);
    assert(0.0 <= p); //and p <= 1.0);

    // setProbability()によってセットされるのが1回目であることをチェック
    assert(nodes[a->index].P[d-1] == -1);

    nodes[a->index].P[d-1] = p;
}


double Torus::getProbability(Node *a, int d) {
    assert(0 <= d and d <= floor(k/2)*n);

    if(d == 0){
        return 1.0;
    }

    // setProbability()によってセットされた値であることをチェック
    assert(nodes[a->index].P[d-1] != -1);

    return nodes[a->index].P[d-1];
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
            std::cout << int(getProbability(a, d)*100)/100.0 << " ";
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


int Torus::route(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    int dist_ct, h;
    double p, p_min;
    Node *neighbor, *max_neighbor;

    //std::cout << "c:";
    //for(int i=0; i<n; i++)
    //    std::cout << c->value[i] << " ";
    //std::cout << std::endl;

    // Base case
    if(c->index == t->index)
        return d;

    dist_ct = distance(c, t);
    p_min = 10000;
    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(hasLink(c, neighbor));
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor) ) {
            p = getProbability(neighbor, dist_ct - 1);

            if (p < p_min) {
                p_min = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(p_min < 10000) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
        return route(c, max_neighbor, t, visited, d + 1);
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(hasLink(c, neighbor));
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (inSpr(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            p = getProbability(neighbor, std::min(dist_ct + 1, diameter));

            if (p < p_min) {
                p_min = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(p_min < 10000) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
        return route(c, max_neighbor, t, visited, d + 1);
    }

    return DELIVERY_FAIL;
}


int Torus::brute(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
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
        if(prev and prev->index == neighbor->index){
            continue;
        }

        if(inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            if(visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
            return brute(c, neighbor, t, visited, d + 1);
        }
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (not hasFaultyLink(c, neighbor)){
            assert(inSpr(neighbor, c, t));
            if (visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
            return brute(c, neighbor, t, visited, d + 1);
        }
    }

    return DELIVERY_FAIL;
}


int Torus::bfs(Node *c, Node *t, std::unordered_map<int, bool> visited) {
    std::queue<int> q_i, q_d;
    int a_d;
    Node *a, *neighbor;

    if(c->index == t->index)
        return 0;

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
        delete[] nodes[i].value;\
        delete[] nodes[i].P;
        delete[] nodes[i].neighbors;
    }
    free(nodes);
}
