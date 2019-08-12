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
#define UNREACHABLE -1
#define DETECT_LOOP -2

int main(int argc, char* argv[]){
    using namespace std;

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    double p_faulty = atof(argv[3]);
    random_device rd;
    mt19937 mt(rd());

    double mean_d_route=0, mean_d_bfs=0;
    int n_route_success=0;
    int n_loop = N;

    bool* route_unreachable = new bool[n_loop];
    bool* route_looping = new bool[n_loop];

    double* d_routes = new double[n_loop];
    double* d_bfss = new double[n_loop];
    double* d_lee =  new double[n_loop];

    uniform_int_distribution<int> dice(0, static_cast<int>(pow(k, n) - 1));

    assert(0.0 <= p_faulty and p_faulty < 1.0);

    // test
    Torus *t = new Torus(n, k);
    t->setFixedFaultyNodes();
    t->calcRoutingProbabilities();
    t->printFaultyLinks();
    t->printProbabilities();
    return 0;
    // ここまでtest

    // パラレルに実行される 同じデータにアクセスしないように注意
#pragma omp parallel for
    for (int i = 0; i < n_loop; i++) {
        //cout << i << endl;

        bool has_non_faulty_route = false;
        while(true){
            // step 1
            Torus *t = new Torus(n, k);
            //t->setRandomFaultyLinks((double)p_faulty);
            t->setRandomFaultyNodes((double)p_faulty);

            // step 2
            int from = dice(mt);
            int to = dice(mt);
            //d_bfss[i] = t->bfs(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>));
            //has_non_faulty_route = (d_bfss[i] != UNREACHABLE);
            //if (not has_non_faulty_route){
            //    delete t;
            //    continue;
           // }

            t->calcRoutingProbabilities();
            d_routes[i] = t->route(0, &(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>), 0);
            route_unreachable[i] = (d_routes[i] == UNREACHABLE);
            route_looping[i] = (d_routes[i] == DETECT_LOOP);

            d_lee[i] = t->distance(&(t->nodes[from]), &(t->nodes[to]));
            delete t;
            break;
        }
    }

    float sum_deviation = 0;
    int n_route_unreachable = 0;
    int n_route_looping = 0;
    // データを集計する
    for(int i=0; i<n_loop; i++){
        if((not route_unreachable[i]) and (not route_looping[i])){
            n_route_success++;
            mean_d_route += d_routes[i];
            if(d_lee[i] > 0){
                sum_deviation += (d_routes[i] - d_lee[i]) / d_lee[i];
            }
        }
        if(route_unreachable[i])
            n_route_unreachable++;
        if(route_looping[i])
            n_route_looping++;

        //assert(d_bfss[i] != UNREACHABLE);
        //mean_d_bfs += d_bfss[i];
    }

    assert(n_route_success == n_loop - n_route_looping - n_route_unreachable);

    //double p_route_success = n_route_success / (double)n_loop;

    cout << p_faulty << ", "
         << n << ", "
         << k << ", "
         << sum_deviation / n_route_success << ","
         //<< (n_route_unreachable)/(double)n_loop << ", "
         //<< n_route_looping/(double)n_loop << ", "
        // << p_brute_success << ", "
       //  << mean_d_route / n_route_success << ", "
       //  << mean_d_brute / n_brute_success << ", "
       //  << mean_d_bfs / n_loop
         << endl;
}
