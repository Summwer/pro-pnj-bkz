#include "est.h"


//dvol = sum(log2(||b_i^*||/sigma))
void gsa_est(int dim, FP_NR<FT> dvol, Params params){
    // printf_input(dim,dvol);
    printf("Generate gs-lengths by GSA assumption...");
    vector<double> l = gen_simulated_gso(dim, dvol);
    est(l, params);
}


void est(vector<double> l, Params params){
    int dim = int(l.size());
    switch(params.method){
        case 1:
            call_enumbs(l,params);
            break;
        case 2:
            call_bssa(l,params,50, dim);
            // call_bssa(l,params,50, int(0.9*(double)dim));
            break;
        // case 3:
        //     call_search_tree(l,params);
        //     break;
        default:
            cout<<"Tere's no method named: "<<params.method<<endl;
    }
}



void lwechal_est(int n, double alpha, Params params){
    LWEchal* lwechal = gen_lwechal_instance(n,alpha);
    // int dim = lwechal->dim;
    // FP_NR<FT> dvol = lwechal->dvol;
    // gsa_est(dim, dvol, params);
    double  sigma = lwechal->alpha * lwechal->q;

    printf("After a sigma normalization,");
    vector<double> l = lwechal->log_rr;
    int dim = int(l.size());
    for(int i = 0; i < dim; i++){
        l[i] -= log2(sigma);
    }

    double slope = get_current_slope(l,0,dim);
    printf("Slope = %f\n", slope);

    est(l,params);
    printf("\n\n\n");
}


// vector<LWEchal*> load_lwechallenge(int n, double alpha){
//     vector<LWEchal*> lwechals;
//     for(int i = 0; i < int(lwes.size());i++){
//         int n = lwes[i].first;
//         double alpha = lwes[i].second;
//         lwechals.insert(lwechals.end(), gen_lwechal_instance(n, alpha));
//     }
//     return lwechals;
// }


// solved lwe challenges: (n,alpha,dim,dvol)
// vector<LWEchal*> load_solved_lwechallenges(){
//     vector<pair<int,double>> solved_lwechals = {{45,0.030}, {50, 0.025}, {55, 0.020}, {60, 0.015}, {85, 0.005}, {90, 0.005}};
//     return load_lwechallenges(solved_lwechals);
// }


// // low-dim lwe challenges: (n,alpha,dim,dvol)
// vector<LWEchal> load_low_dim_lwechallenges(){
//         // {40, 0.035,	188,	327.7388246557668}
//         // {40, 0.025,	172,	331.9735338747315},
//         // {45, 0.020,	185,	373.4658972674150},
//         // {50, 0.015,	194,	415.6552752456852},
//         // {55, 0.010,	205,	495.0168620849964},
//         // {60, 0.010,	222,	522.7192487542407},
//         // {70, 0.005,	235,	641.7748006815514},
//         // {75, 0.005,	252,	678.7288625607595}

//     vector<pair<int,double>> low_dim_lwechals = {{40,0.035}, {40, 0.025}, {45, 0.020}, {50, 0.015}, {55, 0.010}, {60, 0.010}, {70, 0.005}, {75, 0.005}};
//     return load_lwechallenges(low_dim_lwechals);
// }



   

// unsolved lwe challenges: (n,alpha,dim,dvol)
// vector<LWEchal> load_unsolved_lwechallenges(){
//         // {40, 0.045,	195,	302.1993617},
//         // {45, 0.035,	211,	357.0995642},
//         // {50, 0.030,	228,	400.4076907},
//         // {55, 0.025,	241,	439.9769224},
//         // {60, 0.020,	254,	494.0253108},
//         // {65, 0.015,	262,	557.6405653},
//         // {75, 0.010,	281,	637.6057085},
//         // {95, 0.005,	323,	836.9696072}
//     vector<pair<int,double>> unsolved_lwechals = {{40,0.045}, {45, 0.035}, {50, 0.030}, {55, 0.025}, {60, 0.020}, {65, 0.015}, {75, 0.010}, {95, 0.005}};
//     return load_lwechallenges(unsolved_lwechals);

// }

