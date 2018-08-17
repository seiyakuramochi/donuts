#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <random>
#include <stdlib.h>
#include <chrono>
 
#include "torus.hpp"
 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#define N 10000
#define DELIVERY_FAIL (-1)

int main(int argc, char* argv[]){
    using namespace std;

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    double p_faulty = atof(argv[3]);
    random_device rd;
    mt19937 mt(rd());

    double mean_d_route=0, mean_d_brute=0, mean_d_bfs=0;
    int n_brute_success = 0, n_route_success=0, n_bfs_success=0;
    int n_loop = N;

    bool* route_success = new bool[n_loop];
    bool* brute_success = new bool[n_loop];
    bool* bfs_success = new bool[n_loop];

    double* d_routes = new double[n_loop];
    double* d_brutes = new double[n_loop];
    double* d_bfss = new double[n_loop];

    uniform_int_distribution<int> dice(0, static_cast<int>(pow(k, n) - 1));

    // パラレルに実行される 同じデータにアクセスしないように注意
#pragma omp parallel for
    for (int i = 0; i < n_loop; i++) {
        Torus *t = new Torus(n, k);
        t->setRandomFaultyLinks((double)p_faulty);
        t->calcRoutingProbabilities();

        int from = dice(mt);
        int to = dice(mt);

        d_brutes[i] = t->brute(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>), 0);
        brute_success[i] = (d_brutes[i] != DELIVERY_FAIL);

        d_routes[i] = t->route(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>), 0);
        route_success[i] = (d_routes[i] != DELIVERY_FAIL);

        d_bfss[i] = t->bfs(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>));
        bfs_success[i] = (d_bfss[i] != DELIVERY_FAIL);

        delete t;
    }

    // データを集計する
    for(int i=0; i<n_loop; i++){
        if(route_success[i]){
            n_route_success++;
            mean_d_route += d_routes[i];
        }

        if(brute_success[i]){
            n_brute_success++;
            mean_d_brute += d_brutes[i];
        }

        if(bfs_success[i]){
            n_bfs_success++;
            mean_d_bfs += d_bfss[i];
        }
    }

    double p_route_success = n_route_success / double(N);
    double p_brute_success = n_brute_success / double(N);
    double p_bfs_success = n_bfs_success / double(N);

    cout << p_faulty << ", "
         << n << ", "
         << k << ", "
         << p_route_success << ", "
         << p_brute_success << ", "
         << p_bfs_success << ", "
         << mean_d_route / n_route_success << ", "
         << mean_d_brute / n_brute_success << ", "
         << mean_d_bfs / n_bfs_success
         << endl;
}
