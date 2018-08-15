#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <random>
#include <stdlib.h>
#include <chrono>
 
#include "torus.h"
 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#define N 1000
#define DELIVERY_FAIL (-1)

int main(int argc, char* argv[]){
    using namespace std;

    double p_faulty;
    double mean_d_route, mean_d_brute, mean_d_bfs;
    int n, k;

    n = atoi(argv[1]);
    k = atoi(argv[2]);
    p_faulty = atof(argv[3]);
    random_device rd;
    mt19937 mt(rd());

    int n_brute_success = 0, n_route_success=0, n_brute_carried=0,  n_route_carried=0;
    int n_loop = 6*N;

    bool* all_success = new bool[n_loop];
    bool* route_carried = new bool[n_loop];
    bool* brute_carried = new bool[n_loop];
    bool* route_success = new bool[n_loop];
    bool* brute_success = new bool[n_loop];

    double* d_routes = new double[n_loop];
    double* d_brutes = new double[n_loop];
    double* d_bfss = new double[n_loop];

    mean_d_route = 0;
    mean_d_brute = 0;
    mean_d_bfs = 0;

    uniform_int_distribution<int> dice(0, static_cast<int>(pow(k, n) - 1));

    int count_success = 0;

#pragma omp parallel for
    for (int i = 0; i < n_loop; i++) {
        //cout << i << endl;
        all_success[i] = false;
        brute_carried[i] = false;
        route_carried[i] = false;

        if(count_success < N+1) {
            // auto start = std::chrono::system_clock::now();
            Torus *t = new Torus(n, k);
            // auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
            // auto dur = end - start;
            // auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            // //std::cout << msec << "milli sec (init) \n";
            // start = std::chrono::system_clock::now();


            t->setRandomFaultyLinks(p_faulty);

            //
            // end = std::chrono::system_clock::now();       // 計測終了時刻を保存
            // dur = end - start;
            // msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            // //std::cout << msec << "milli sec (setrandomfaultylink) \n";
            // start = std::chrono::system_clock::now();


            t->calcRoutingProbabilities();


            // end = std::chrono::system_clock::now();       // 計測終了時刻を保存
            // dur = end - start;
            // msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            // //std::cout << msec << "milli sec (P) \n";
            // start = std::chrono::system_clock::now();


            int from = dice(mt);
            int to = dice(mt);

            int d_brute = t->brute(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>), 0);
            brute_carried[i] = true;
            if (d_brute != DELIVERY_FAIL) {
                brute_success[i] = true;
            } else {
                brute_success[i] = false;
                delete t;
                continue;
            }

            int d_route = t->route(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>), 0);
            route_carried[i] = true;
            if (d_route != DELIVERY_FAIL) {
                route_success[i] = true;
            } else {
                route_success[i] = false;
                delete t;
                continue;
            }

            int d_bfs = t->bfs(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>));
            assert(d_bfs != DELIVERY_FAIL);

            d_bfss[i] = d_bfs;
            d_brutes[i] = d_brute;
            d_routes[i] = d_route;
            all_success[i] = true;
            count_success++;
            //t->printAdjMatrix();
            delete t;

            //
            // end = std::chrono::system_clock::now();       // 計測終了時刻を保存
            // dur = end - start;
            // msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            // //std::cout << msec << "milli sec (routing) \n";

            return 0;
        }
    }

    count_success = 0;
    for(int i=0; i<n_loop; i++){
        if(all_success[i]){
            if(count_success == N)
                break;
            mean_d_bfs += d_bfss[i];
            mean_d_brute += d_brutes[i];
            mean_d_route += d_routes[i];

            count_success++;
        }

        if(route_carried[i]) {
            n_route_carried++;
            if(route_success[i])
                n_route_success++;
        }

        if(brute_carried[i]) {
            n_brute_carried++;
            if(brute_success[i])
                n_brute_success++;
        }

        if(i == n_loop-1)
            cout << "FAIL!" << endl;
    }

    double p_route_success = n_route_success / double(n_route_carried);
    double p_brute_success = n_brute_success / double(n_brute_carried);

    cout << p_faulty << ", "
         << n << ", "
         << k << ", "
         << p_route_success << ", "
         << p_brute_success << ", "
         << mean_d_bfs / N << ", "
         << mean_d_route / N << ", "
         << mean_d_brute / N
         << endl;
}
