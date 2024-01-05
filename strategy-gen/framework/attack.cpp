#include "attack.h"


void call_enumbs(vector<double> l, Params* params){
    int d = int(l.size());
    EnumBS* enumbs = new EnumBS(params,d);

    
    if(params->threads == 1){
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS)"<<endl;
    }else if (params->threads > 1){
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS in parallel)"<<endl;
    }
    else{
        cerr<<"Set bad threads = "<<params->threads<<"."<<endl;
    }
    printf("v1. beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, G_prec = %e,  slope_prec = %e,  progressive_sieve = True, worst_case = %d, succ_prob = %1.3f, ", params->beta_start, params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_G_prec, params->enumbs_slope_prec,  params->worst_case, params->succ_prob);
    if(params->enumbs_min_G)
        printf("Find minimal time cost strategy, ");
    else
        printf("Find a strategy below the maximal RAM = %f log2(bit), ", params->max_RAM);
    if(params->worst_case)
        printf("worst_case, ");
    else
        printf("average_case, ");
    
    if(params->cost_model == 1){
        printf("theo_pnjbkz_d4f = %d, ", params->theo_pnjbkz_d4f);
        // if(params->est_model !=3)
        printf("theo_pump_d4f = %d, ", params->theo_pump_d4f);
        // else
            // printf("theo_pump_d4f: compute ||pi_f(target_vector)||<= sqrt(4/3) GH(L_f)\n");
        printf("list_decoding = `%s`. \n\n", params->list_decoding.c_str());
    }
    if(params->cost_model == 2){
        printf("practical_pnjbkz_d4f = %d, ", params->practical_pnjbkz_d4f);
        // if(params->est_model !=3)
        printf("practical_pump_d4f = %d.\n ", params->practical_pump_d4f);
        // else
            // printf("practical_pump_d4f: compute ||pi_f(target_vector)||<= sqrt(4/3) GH(L_f)\n");
    }
    
    // if(params->enum_add_G2)
    // printf("Min G2 Strategy. \n");
    // else
    //     printf("Min G Strategy. \n");
    auto start = system_clock::now();
    if(params->threads == 1)
        enumbs->enumbs_est(l);
    else
        enumbs->enumbs_est_in_parallel(l);
    auto finish = system_clock::now();
    duration<double> diff = finish - start;
    cout<<"EnumBS cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
    
}

// void call_search_tree(vector<double> l, Params* params){
//     SearchTree* search_tree = new SearchTree(params);

//     if(params->threads == 1){
//         cout<<" Attack Estimation via simulation + probabilistic model (Blocksize Strategy Search Tree)"<<endl;
//     }else if (params->threads > 1){
//         cout<<" Attack Estimation via simulation + probabilistic model (Blocksize Strategy Search Tree in parallel)"<<endl;
//     }
//     else{
//         cerr<<"Set bad threads = "<<params->threads<<"."<<endl;
//     }
//     // printf("beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, G_prec = %e,  slope_prec = %e,  progressive_sieve = %d, ", params->beta_start, params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_G_prec, params->enumbs_slope_prec,  params->worst_case);
//     if(params->worst_case)
//         printf("worst_case, ");
//     else
//         printf("average_case, ");
    
//     // if(params->enum_add_G2)
//     //     printf("Min G2 Strategy. \n");
//     // else
//     //     printf("Min G Strategy. \n");
//     auto start = system_clock::now();
//     if(params->threads == 1)
//         search_tree->search_tree_est(l);
//     else
//         search_tree->search_tree_est_in_parallel(l);
//     auto finish = system_clock::now();
//     duration<double> diff = finish - start;
//     cout<<"Blocksize Strategy Search Tree cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
// }

void call_bssa(vector<double> l, Params* params, int sbeta, int gbeta){
    int d = int(l.size());
    BSSA* bssa = new BSSA(params,d);
    
    auto start = system_clock::now();
    cout<<"v1.  Attack Estimation via simulation + probabilistic model (BSSA)"<<endl;
    printf("beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, G_prec = %e,  slope_prec = %e,  progressive_sieve = True, worst_case = %d, succ_prob = %1.2f, ", params->beta_start, params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_G_prec, params->enumbs_slope_prec,  params->worst_case, params->succ_prob);
    if(params->enumbs_min_G)
        printf("Find minimal time cost strategy, ");
    else
        printf("Find a strategy below the maximal RAM = %f log2(bit), ", params->max_RAM);
    if(params->worst_case)
        printf("worst_case, ");
    else
        printf("average_case, ");
    if(params->cost_model == 1){
        printf("theo_pnjbkz_d4f = %d, theo_pump_d4f = %d, ", params->theo_pnjbkz_d4f, params->theo_pump_d4f);
        printf("list_decoding = `%s`. \n\n", params->list_decoding.c_str());
    }
    if(params->cost_model == 2){
        printf("practical_pnjbkz_d4f = %d, practical_pump_d4f = %d, ", params->practical_pnjbkz_d4f, params->practical_pump_d4f);
    }
    if(params->bssa_tradion == 1){
        printf("traditonal BSSA \n");
    }
    if(params->bssa_tradion == 0){
        printf("improved BSSA \n");
    }
    // if(params->mul_node)
    //     bssa->bssa_est_mul_node(l, sbeta, gbeta);
    // else
    bssa->bssa_est(l, sbeta, gbeta);
    auto finish = system_clock::now();
    duration<double> diff = finish - start;
    cout<<"BSSA cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
}


