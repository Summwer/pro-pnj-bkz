
#include "../framework/est.h"





// params input in main function
// argv[0]: implemented file name
// argv[1]: method -- 1:enumba; 2: bssa
int main(int argc,char **argv){
    Params params = new Params; //J, gap, J_gap, cost_model, verbose,
    params.threads = 95;
    params.cost_model = 1; //sec model;
    // params.progressive_sieve = true;
    params.verbose = true;
    params.debug = false;
    params.worst_case = false;
    params.method = atoi(argv[1]); //1:enumbs;2:bssa
    params.gap = 1;
    params.J = 100; 
    params.J_gap = 3;
    params.enumbs_G_prec = 0.001;
    params.max_loop = 1; 
    params.max_dim = 1500; 
    params.min_G_prec = 0.001;
    if(atoi(argv[2]) == 1)
        params.list_decoding = "agps20"; //"matzov22"
    if(atoi(argv[2]) == 2)
        params.list_decoding = "matzov22"; //"matzov22"
    // params.bssa_tradition = true;  
    

    int n, m, q,eta;
    map<int,double> D_s,D_e;
    LWEchal* lwechal;

    // Kyber-1024 round-3 parameters
    printf("============= Kyber-1024\n");
    n = 1024, m = 1024, q = 3329;
    D_s = build_centered_binomial_law(2);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);

    
    // Kyber-512 round-3 parameters
    printf("============= Kyber-512\n");
    n = 512, m = 512, q = 3329;
    D_s = build_centered_binomial_law(3);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);

    // Kyber-768 round-3 parameters
    printf("============= Kyber-768\n");
    n = 768, m = 768, q = 3329;
    D_s = build_centered_binomial_law(2);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);


    /*----------------Instance Generation-----------------*/
    //Dilithium-I round-3 parameters
    printf("============= Dilithium-I\n");
    n = 4*256, m = 4*256, q = 8380417, eta = 2;
    // map<int,rational<int>> D_s,D_e;
    // rational<int> one(1);
    // for(int x=-eta; x<=eta; x++){
    //     D_s[x] = one/(2*eta+1);
    //     D_e[x] = one/(2*eta+1);
    // }
    D_s={},D_e={};
    for(int x=-eta; x<=eta; x++){
        D_s[x] = 1./(2*eta+1);
        D_e[x] = 1./(2*eta+1);
    }
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);



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
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);

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
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params.verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);



    return 1;
}


