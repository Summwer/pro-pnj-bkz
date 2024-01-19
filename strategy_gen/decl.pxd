# -*- coding: utf-8 -*-

#from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
#from fplll cimport Z_NR, ZZ_mat, FP_NR
#from libgmp cimport mpz_t
#from libmpfr cimport mpfr_t




cdef extern from "framework/enumbs.h" nogil:

    #ctypedef mpz_t ZT
    #ctypedef mpfr_t FT

    cdef struct Params:
        double succ_prob
        int J 
        int gap 
        int J_gap 
        int cost_model
        bool verbose 
   
        int threads 
        int max_dim
        double max_num 
        double max_RAM 
        int max_loop 

        int method 

        bool debug
        bool verification 

        #params for cost model;
        string list_decoding

        #enumbs params
        int delta_beta 
        double enumbs_G_prec
        double enumbs_slope_prec 
        int beta_start 
        bool worst_case 
        bool enumbs_min_G 
        double min_G_prec 
        bool print_Gcums 

        #bssa params
        bool bssa_tradion
        #bool mul_node  
        #bool beta_gap 

        #params for pnj-bkz
        int theo_pnjbkz_d4f 
        int practical_pnjbkz_d4f 
        int compute_jub 

        #params for last pump
        int theo_pump_d4f 
        int practical_pump_d4f


    #cdef cppclass BKZJSim:
    #    BKZJSim(Params params, int d);
        
    #cdef cppclass COST:
    #    COST(Params params);
      

    
    #cdef struct LWEchal:
    #    int n
    #    double alpha
    #    double sigma
    #    vector[double] log_rr
    #   int q
    #    int m
    #    vector[Z_NR[ZT]] c 
    #    ZZ_mat[ZT] A
    #    ZZ_mat[ZT] B

    #    int dim 
    #    FP_NR[FT] dvol 


    #cdef LWEchal* gen_lwechal_instance(int n, double alpha);

#cdef extern from "framework/enumbs.h" nogil:
    cdef cppclass EnumBS:
        ##BKZJSim this->sim;
        ##COST cost;
        #int d
        ##Params params
        
        EnumBS(Params params,int d)

        #EnumBS(int d)
        #EnumBS()
        
        int strategy_size()
        #void enumbs_est_in_parallel()
        void enumbs_est_in_parallel(double* l0)
        #void enumbs_est_in_parallel(int dim, double dvol)
        void get_strategy(long* strategy)
        #void set_threads(int nr)
        #void print_param_setting()