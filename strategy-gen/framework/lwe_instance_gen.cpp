// #include "instance_gen.h"
// #include "utils.h"
#include "lwe_estimation.h"

#include <boost/rational.hpp>



void find_min_m(LWEchal* lwechal){
    
    printf("-------------------------\n");
    printf("Primal attack, TU LWE challenge n=%d, m=%d, sigma=%3.5f alpha=%e, q=%d \n" , lwechal->n, lwechal->m - lwechal->n, lwechal->sigma,lwechal->alpha, lwechal->q);
    cout<<"===========1============"<<endl;
    // tuple<int,int,int> min_cost_param = gsa_params(lwechal);
    cout<<"======================="<<endl;
    // int b = get<0>(min_cost_param), s = get<1>(min_cost_param);
    // int m = get<2>(min_cost_param);
    int d = lwechal->m  + 1;
    // lwechal->m = m;
    lwechal->dim = d;
   
    // printf("Chose %d samples. Predict solution at bkz-%d + svp-%d. \n", m, b, s);
    
    

    primal_lattice_basis(lwechal);

    // ostringstream os;
    // lwechal->B[0].print(os);
	// cout<<os.str()<<endl;

    // double  sigma = alpha * q; 

    ZZ_mat<ZT> U;
    ZZ_mat<ZT> UT;
    MatGSO<Z_NR<ZT>, FP_NR<FT>> M(lwechal->B, U, UT, GSO_DEFAULT);
    M.update_gso();
    lwechal->log_rr.resize(d);
    for(int i = 0; i < d; i++){
        FP_NR<FT> tmp;
        M.get_r(tmp,i,i);
        lwechal->log_rr[i] = log2(tmp.get_d())/2.;
    }
    double slope = get_current_slope(lwechal->log_rr,0,d);
    cout<<"Initial slope = "<<slope<<endl;
    lwechal->dvol = M.get_log_det(0,d)/2. - log(lwechal->sigma)*d;

    printf("dim = %d, dvol = %3.11f\n\n", lwechal->dim, lwechal->dvol.get_d());
   
    // return lwechal;
}



int main(){

    // FP_NR<mpfr_t>::set_prec(32); 
    bool verbosity = true;
    /*----------------Instance Generation-----------------*/
    //Dilithium-I round-3 parameters
    printf("============= Dilithium-I\n");
    int n = 4*256, m = 4*256, q = 8380417, eta = 2;
    // map<int,rational<int>> D_s,D_e;
    // rational<int> one(1);
    // for(int x=-eta; x<=eta; x++){
    //     D_s[x] = one/(2*eta+1);
    //     D_e[x] = one/(2*eta+1);
    // }
    map<int,double> D_s={},D_e={};
    for(int x=-eta; x<=eta; x++){
        D_s[x] = 1./(2*eta+1);
        D_e[x] = 1./(2*eta+1);
    }
    LWEchal* lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);



    //Dilithium-II round-3 parameters
    printf("============= Dilithium-II\n");
    n = 5*256, m = 6*256, q = 8380417, eta = 4;
    D_s={},D_e={};
    // for(int x=-eta; x<=eta; x++){
    //     D_s[x] = one/(2*eta+1);
    //     D_e[x] = one/(2*eta+1);
    // }
    for(int x=-eta; x<=eta; x++){
        D_s[x] = 1./(2*eta+1);
        D_e[x] = 1./(2*eta+1);
    }
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);



    //Dilithium-III round-3 parameters
    printf("============= Dilithium-III\n");
    n = 7*256, m = 8*256, q = 8380417, eta = 2;
    D_s={},D_e={};
    // for(int x=-eta; x<=eta; x++){
    //     D_s[x] = one/(2*eta+1);
    //     D_e[x] = one/(2*eta+1);
    // }
    for(int x=-eta; x<=eta; x++){
        D_s[x] = 1./(2*eta+1);
        D_e[x] = 1./(2*eta+1);
    }
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);



    // Kyber-512 round-3 parameters
    printf("============= Kyber-512\n");
    n = 512, m = 512, q = 3329;
    D_s = build_centered_binomial_law(3);
    D_e = build_centered_binomial_law(2);
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);


    // Kyber-768 round-3 parameters
    printf("============= Kyber-768\n");
    n = 768, m = 768, q = 3329;
    D_s = build_centered_binomial_law(2);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);


    // Kyber-1024 round-3 parameters
    printf("============= Kyber-1024\n");
    n = 1024, m = 1024, q = 3329;
    D_s = build_centered_binomial_law(2);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, verbosity);

    return 0;
}


