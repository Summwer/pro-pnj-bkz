#include <iostream>
#include <map>
#include <vector>
#include <boost/rational.hpp>
#include<bits/stdc++.h>
#include <fplll.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <chrono>


using namespace std;
using namespace boost;
using namespace fplll;
using namespace std::chrono; 

#ifndef ZT
#define ZT mpz_t
#endif

#ifndef FT
#define FT mpfr_t
#endif


#ifndef MAX_DIM
#define MAX_DIM 1000
#endif



#ifndef PREC
#define PREC 1e3
#endif


struct Params{
    double succ_prob = 0.999; //0.999(default)
    int J =  100; //J -- maximal jump value;
    int gap = 1; //gap -- gap of each beta;
    int J_gap = 1; //J_gap -- gap of each jump value;
    //cost_model: 1: gate model; 2: sec model(in qd) with threads=32, gpus = 2; 3: sec model(in dd) with threads=32, gpus = 2;
    int cost_model = 2; 
    int compute_jub = 2; //Compute jump upper bound in 1:pessimistic d4f(B); 2: optimistic d4f(B); 3: d4f(beta)/2; 4: min(d4f(beta),0.1beta)


    bool verbose = false; //print logging or not
    //progressieve_sieve: True: progressieve sieve; False: normal sieve
    
    int threads = 10;
    int max_dim = MAX_DIM; //set the maximal blocksize to find the expected minimal  strategy
    double max_num = 1e3;
    double max_RAM = 1000; //43.58; //=1.5TB , //no-limit, 300
    int max_loop = 5; // 15; //set the maximal loop for one blocksize to find the expected minimal  strategy1

    int method = 1; //1: enumbs estimation; 2: bssa estimation

    bool debug = false; //print debug logging or not.
    bool verification =false; //verify the correctness of strategy

    //params for cost model;
    string list_decoding = "matzov22"; //list_decoding = "agps20" or "matzov22"
    // int est_model = 2; //1, cumulated pnj-BKZ + cumulated proba pump; 2, cumulated pnj-BKZ + succ-fail proba pump;  (deleted)3. cumulated pnj-BKZ + cumulated proba pump with d4f: ||pi_f(s)||<= sqrt(4/3) GH(L_f)

    //enumbs params
    int delta_beta = -1;//10; //Ensure quality(beta,jump) >= quality(beta-delta_beta,1)
    double enumbs_G_prec = 1e-3;//1e-3; //1e-5; //set the precision of enumbs
    double enumbs_slope_prec = 1e-4;//1e-6;
    int beta_start = 50;  
    bool worst_case = false;
    // bool enum_add_G2 = true;
    bool enumbs_min_G = true; //if enumbs_min_G = false, then call a min RAM enumbs
    double min_G_prec = 0.1;
    bool print_Gcums = false; //The increased Gcum in pump

    //bssa params
    bool bssa_tradion = true;
    // bool mul_node  = false;
    // bool beta_gap = 1;

    //params for pnj-bkz
    int theo_pnjbkz_d4f = 2; //the dim4free function for pnjbkz in theoretical cost mode. 1: theo d4f1; 2: theo d4f2; 3: d4f in default g6k
    int practical_pnjbkz_d4f = 3; //the dim4free function for pnjbkz in prectical cost mode. 1: theo d4f1; 2: theo d4f2; 3: d4f in default g6k; 
    

    //params for last pump
    // bool progressive_sieve =  true; 
    int theo_pump_d4f = 2; //the dim4free function for last pump in theoretical cost mode. 1: theo d4f1; 2: theo d4f2; 3: d4f in default g6k
    int practical_pump_d4f = 3; //the dim4free function for last pump in prectical cost mode. 1: theo d4f1; 2: theo d4f2; 3: d4f in default g6k; 4: optimistic d4f(B);
};


struct LWEchal{
    int n;
    double alpha;
    double sigma; //average variance for e
    vector<double> log_rr; //log2(||b_i^*||)
    int q;
    int m;
    vector<Z_NR<ZT>> c = {};
    ZZ_mat<ZT> A;
    ZZ_mat<ZT> B; //primal basis

    int dim = 0;
    FP_NR<FT> dvol = 0.;
    
};




void print_map(map<int, double> mp);
void print_map(map<int, FP_NR<FT>>  mp);
void printf_red(const char *s);
void printf_input(int d, FP_NR<FT> dvol);
void print_vector(vector<double> v,int index_start=0, int index_end=-1);
void print_vector(vector<int> v,int index_start=0, int index_end=-1);
void print_vector(vector<Z_NR<ZT>> v,int index_start=0, int index_end=-1);
void print_matrix(vector<vector<Z_NR<ZT>>> matrix);
void print_vector(vector<FP_NR<FT>> v,int index_start=0, int index_end=-1);
void print_vector(vector<pair<int,int>> v,int index_start=0, int index_end=-1);
void print_vector(vector<pair<double,double>> v,int index_start=0, int index_end=-1);
void print_vector(vector<tuple<double,double,double>> v,int index_start=0, int index_end=-1);
void print_vector(vector<tuple<int,int,int>> v,int index_start=0, int index_end=-1);

//void print_matrix(ZZ_mat<ZT> matrix);

pair<double,double> average_variance(std::map<int,double> D);
vector<int> KeySet(map<int, double> mp);
vector<double> ValueSet(map<int, double> mp);
int draw_from_distribution(std::map<int,double> D, int sample_num=1000);
vector<int> expand_samples(vector<int> samples,vector<double> prob,int sample_num);


// void gen_samples(ZZ_mat<ZT> &matrix, int m, int n, int q);
void build_LWE_lattice(ZZ_mat<ZT> &matrix, ZZ_mat<ZT> A, int q);
void kannan_embedding(ZZ_mat<ZT> &matrix, vector<Z_NR<ZT>> target, int factor=1);

FP_NR<FT> compute_delta(int beta);
double get_current_slope(vector<double> l, int start, int end);
FP_NR<FT> bkzgsa_gso_len( FP_NR<FT> logvol, int i, int d, int beta=-1);//FP_NR<FT> delta,
vector<double> gen_simulated_gso(int d, FP_NR<FT> logvol);

void simulator_test(int d, FP_NR<FT> logvol); //test simulator 

double gaussian_heuristic(vector<FP_NR<FT>> l, int index_start);
double gaussian_heuristic_log2(vector<FP_NR<FT>> l, int index_start);
double gaussian_heuristic_log2(vector<double> l, int index_start);

int dims4free(int beta); //leaky-lwe-estimator
int wrapper_default_dim4free_fun(int beta);
int default_dim4free_fun(int beta);
int theo_dim4free_fun1(int beta);//Params params,
int theo_dim4free_fun2(int beta);
int get_beta_from_sieve_dim(int sieve_dim, int d, int choose_dims4f_fun);
int get_f_for_pnjbkz(Params params, int beta);
int theo_dim4free1_in_B(vector<double> l);
int theo_dim4free2_in_B(vector<double> l);
int jump_upper_bound(Params params, int beta, vector<double> l);

// int get_f_for_pump(Params params, int beta);
int get_f_for_pump(Params params, int beta, vector<double> l);
int get_beta_(Params params,int beta, int jump, int d);

std::map<int,double>  build_centered_binomial_law(int k);