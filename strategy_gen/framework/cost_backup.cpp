#include "cost.h"

// // Return log2 of the number of gates for FindAllPairs according to AGPS20
// FP_NR<FT> COST::agps20_gates(int beta_prime){
//     FP_NR<FT> k = ((double)beta_prime) / 8.;
//     FP_NR<FT> d1,d2,x;
//     if(k != round(k.get_d())){
//         x = k - floor(k);
//         d1 = agps20_gates(8*(int)floor(k).get_d());
//         d2 = agps20_gates(8*((int)floor(k).get_d() + 1));
//         return x * d2 + (1. - x) * d1;
//     }
//     return COST::agps20_gate_data[beta_prime];
// }


// Return log2 of the number of gates for FindAllPairs according to AGPS20
FP_NR<FT> COST::agps20_gates(int beta_prime){
    FP_NR<FT> k = ((double)beta_prime) / 8.;
    FP_NR<FT> d1,d2,x;
    if(k != round(k.get_d())){
        x = k - floor(k);
        d1 = agps20_gates(8*(int)floor(k).get_d());
        d2 = agps20_gates(8*((int)floor(k).get_d() + 1));
        return x * d2 + (1. - x) * d1;
    }
    return COST::agps20_gate_data[beta_prime];
}



// Function C from AGPS20 source code
double caps_vol(int d, double theta){
    /*
    The probability that some v from the sphere has angle at most theta with some fixed u.

    :param d: We consider spheres of dimension `d-1`
    :param theta: angle in radians
    :param: compute via explicit integration

    */
    return alglib::incompletebeta((d-1.)/2., 1./2., pow(sin(theta),2))/2.;

}

// Return log2 of the number of vectors for sieving according to AGPS20
double agps20_vectors(int beta_prime){
    double N = 1./caps_vol(beta_prime, M_PI/3.);
    return log2(N);
}


//cost of bkz with progressive sieve
pair<double,double> COST::theo_bkz_cost(int n, int beta,int J){
    /*Return cost of bkz-beta in theoretical Gate: T/G -- time cost, B -- memory cost*/
    double gates = 0.;
    if(beta <=10)
        return make_pair(0.,0.);
    // int f = get_f_for_pnjbkz(params,beta);
    int f = get_f_for_pnjbkz(params,beta);
    int beta_prime = floor(beta - f);
    if(params.list_decoding == "matzov22")
        gates = log2((((double)n -beta)/J)*COST::C) +  list_decoding_classical_matzov22.first * beta_prime + list_decoding_classical_matzov22.second; //agps20_gates(beta_prime).get_d();
    if(params.list_decoding == "agps20"){
        if(beta_prime < 64 or beta < beta_prime)
            return make_pair(0.,0.);
        else if(beta_prime > 1024)
            gates = log2((((double)n-beta)/J)*COST::C) +  list_decoding_classical_agps20.first * beta_prime + list_decoding_classical_agps20.second;
        else
            gates = log2((((double)n-beta)/J)*COST::C)  + agps20_gates(beta_prime).get_d();
    }
    double bits = log2(8*beta_prime) + agps20_vectors(beta_prime);
    return make_pair(gates, bits);
    
}



pair<double,double> COST::theo_pump_cost(int beta){
    /*Return cost of pump-beta in theoretical Gate: T/G -- time cost, B -- memory cost*/
    double gates = 0.;
    if(beta <=10)
        return make_pair(0.,0.);
    
    if(params.list_decoding == "matzov22")
        gates = log2(COST::C) +  list_decoding_classical_matzov22.first * beta + list_decoding_classical_matzov22.second; //agps20_gates(beta).get_d();
    if(params.list_decoding == "agps20"){
        if(beta < 64)
            return make_pair(0.,0.);
        else if(beta > 1024)
            gates = log2(COST::C) + list_decoding_classical_agps20.first * beta + list_decoding_classical_agps20.second;
        else
            gates = log2(COST::C) + agps20_gates(beta).get_d();
    }
    double bits = log2(8.*beta) + agps20_vectors(beta);
    return make_pair(gates, bits);
    
}








//threads = 32, gpus = 2,  pnj-bkz cost in qd precision type
pair<double,double> COST::get_k1_k2_pnj_qd(int beta,bool sieve){
    double k1,k2;
    if(beta >=0 and beta <10 and not sieve){
        k1 = -1;
        k2 = -1;
    }
    else if(beta>=10 and beta<=42 and not sieve){
        // k1 = 0.03;
        // k2 = -2.317327;
        k1 = 0.006000;
        k2 = 2.390000;
    }
    else if(beta < 50 and not sieve){
        // k1 = 0.202385;
        // k2 = -9.340418;
        k1 = 0.033000;
        k2 = 1.320000;
    }
    else if(beta <= 97 and sieve){
        k1 = 0.056;
        k2 = 7.85;
    }
    else if(beta <= 118 and sieve){
        k1 = 0.215;
        k2 = -7.61;
    }
    else if(beta <= 128 and sieve){
        k1 = 0.314;
        k2 = - 19.24;
    }
    else{
        k1 = 0.368;
        k2 = -26.15;
    }
    return make_pair(k1,k2);
}




//threads = 32, gpus = 2, pump cost paramter in qd precision type
pair<double,double> COST::get_k1_k2_pnj_dd(int beta){
    double k1,k2;
    if(beta >=0 and beta <10){
        k1 = -1 ;
        k2 = -1 ;
    }
    else if(beta>=10 and beta<=50){
        k1 = 0.030;
        k2 = 6.41;
    }
    else if(beta>50 and beta<=65){
        k1 = 0.070;
        k2 = 4.42;
    }
    else if(beta>65 and beta <= 83){
        k1 = 0.063;
        k2 = 5.47;
    }
    else if(beta <= 94){
        k1 = 0.146;
        k2 = -1.44;
    }
     else if(beta <= 118){
        k1 = 0.215;
        k2 = -7.61;
    }
    else if(beta <= 128){
        k1 = 0.314;
        k2 = - 19.24;
    }
    else{
        k1 = 0.368;
        k2 = -26.15;
    }
    return make_pair(k1,k2);
}



//threads = 32, gpus = 2, pump cost paramter in qd precision type
pair<double,double> COST::get_k1_k2_pump_qd(int beta){
    double k1,k2;
    if(beta >=0 and beta <10){
        k1 = -1 ;
        k2 = -1 ;
    }
    else if(beta>=10 and beta<=60){
        k1 = 0.035657;
        k2 = -2.317327;
    }
    else if(beta <= 96){
        k1 = 0.078794;
        k2 = -0.039742;
    }
    else if(beta <= 116){
        k1 = 0.231927;
        k2 = -14.713430;
    }
    else if(beta <= 128){
        k1 = 0.314;
        k2 = -24.21;
    }
    else{
        k1 = 0.368;
        k2 = -31.12;
    }
    
    return make_pair(k1,k2);
}





// //get pump cost in threads = 20
// pair<double,double> COST::practical_pump_cost_qd(int beta){
//     //make sure not use the enum cost 
//     int f = get_f_for_pump(params, beta);
//     // f = dims4free(beta);
//     int beta_prime = beta - f;
//     double secs, bits;
//     pair<double,double> k = get_k1_k2_pump_qd(beta_prime); // threads = 20
//     double k1 = k.first, k2 = k.second;
//     // k = (1/71.)*((1.33)**(beta/10.));
    
//     secs = k1*((double) beta_prime)+k2; 

//     //unit: GB
//     if( beta_prime <= 56)
//     	bits = 2.0311;
//     else if( beta_prime >=57 and  beta_prime <=63)
//     	bits = 8e-5 *  pow(beta_prime,2) - 0.0083*beta_prime + 2.2555;
//     else if( beta_prime >= 64 and beta_prime <= 94)
//     	bits = 2.195202e-6 * pow(beta_prime,4) - 6.297613e-4 * pow(beta_prime,3) +6.803540e-2 * pow(beta_prime,2) - 3.274476 * beta_prime + 61.31963;
//     else if( beta_prime >= 95)
//     	bits = 2.0311 + pow(2,0.1992* beta_prime  - 18.714);

//     //unit: log2(bit)
//     bits = log2(bits * pow(2,33));
//     return make_pair(secs,bits); //n_expected = beta -f , beta = d-llb
// }
    


//get pump cost in threads = 32, gpus = 2
pair<double,double> COST::practical_pump_cost_qd(int beta){
    //make sure not use the enum cost 
    // int f = get_f_for_pump(params, beta);
    // f = dims4free(beta);
    int beta_prime = beta;
    double secs, bits;
    pair<double,double> k = get_k1_k2_pump_qd(beta_prime); // threads = 20
    double k1 = k.first, k2 = k.second;
    // k = (1/71.)*((1.33)**(beta/10.));
    
    secs = k1*((double) beta_prime)+k2; 

    //unit: GB
    if( beta_prime <= 56)
    	bits = 2.0311;
    else if( beta_prime >=57 and  beta_prime <=63)
    	bits = 8e-5 *  pow(beta_prime,2) - 0.0083*beta_prime + 2.2555;
    else if( beta_prime >= 64 and beta_prime <= 94)
    	bits = 2.195202e-6 * pow(beta_prime,4) - 6.297613e-4 * pow(beta_prime,3) +6.803540e-2 * pow(beta_prime,2) - 3.274476 * beta_prime + 61.31963;
    else if( beta_prime >= 95)
    	bits = 2.0311 + pow(2,0.1992* beta_prime  - 18.714);

    //unit: log2(bit)
    bits = log2(bits * pow(2,33));
    return make_pair(secs,bits); //n_expected = beta -f , beta = d-llb
}


//get pnj-BKZ time test in threads = 32, gpus= 2 in qd precision
double COST::practical_bkz_cost_qd(int d,int beta,int jump){
    // int extra_dim4free = 12;
    // beta = beta + extra_dim4free;
    // f = f + extra_dim4free;
    int f = 0;
    bool sieve;
    if(beta < 50){
        sieve = false;
        f = 0;
    }
    else{
        sieve = true;  
        // f = get_f_for_pnjbkz(params, beta);
        f = get_f_for_pnjbkz(params, beta);
    }
    pair<double,double> k = get_k1_k2_pnj_qd(beta-f,sieve); // threads = 20
    double k1 = k.first, k2 = k.second;
    double c3= 0.018, c4 = -2.24;

    // if(beta - f <= 60)
    //     return (k1*(beta-f)+k2) - log2(jump);
    // else
    //if(c3*d+c4 > 1)
    if(c3*d+c4>=1)
        return (k1*(beta-f)+k2) + log2(c3*d+c4) - log2(jump);
    else
        return (k1*(beta-f)+k2) - log2(jump);
}



//threads = 32, gpus = 2, pump cost paramter in dd precision type
pair<double,double> COST::get_k1_k2_pump_dd(int beta){
    double k1,k2;
    if(beta >=0 and beta <10){
        k1 = -1 ;
        k2 = -1 ;
    }
    else if(beta>=10 and beta<=60){
        k1 = 0.035657;
        k2 = -2.317327;
    }
    else if(beta <= 86){
        k1 = 0.064883;
        k2 = - 0.172902;
    }
    else if(beta <= 109){
        k1 = 0.196298;
        k2 = - 11.518570;
    }
    else if(beta <= 124){
        k1 = 0.260275;
        k2 =  - 18.491171;
    }
    else{
        k1 = 0.359131;
        k2 = - 30.787325;
    }
    
    return make_pair(k1,k2);
}



// //get pnj-BKZ time test in threads = 20
// double COST::practical_bkz_cost_dd(int d,int beta,int jump){
//     // int extra_dim4free = 12;
//     // beta = beta + extra_dim4free;
//     // f = f + extra_dim4free;
//     int f = 0;
//     int extra_dim4free =12;
//     // bool sieve;
//     if(beta < 50){
//         // sieve = false;
//         f = 0;
//         return practical_bkz_cost_qd(d, beta, jump);
//     }
//     else{
//         // sieve = true;  
//         // f = get_f_for_pnjbkz(params, beta);
//         f = get_f_for_pnjbkz(params, beta);
//     }
//     if(beta-f<=60)
//         return practical_bkz_cost_qd(d, beta, jump);
//     pair<double,double> k = get_k1_k2_pump_dd(beta-f); // threads = 20
//     double k1 = k.first, k2 = k.second;
//     double c3= 0.018, c4 = -2.24;

//     // if(beta - f <= 60)
//     //     return (k1*(beta-f)+k2) - log2(jump);
//     // else
//     //if(c3*d+c4 > 1)
//     // cout<<log2(d-beta+2*f+extra_dim4free)<<endl;
//     // int indices_num = 2*(floor((f+extra_dim4free)/jump)+1) + (floor((d-beta-extra_dim4free)/jump)+1) + 1;
//     // int indices_num = 2*(floor((f+extra_dim4free)/jump)+1) + (floor((d-beta-extra_dim4free)/jump)+1) + 1;
//     if(c3*d+c4>=1)
//         return (k1*(beta-f)+k2) + log2(c3*d+c4) + log2(d-beta+2*f+extra_dim4free) - log2(jump);
//     else
//         return (k1*(beta-f)+k2) + log2(d-beta+2*f+extra_dim4free) - log2(jump);
// }



//get pnj-BKZ time test in threads = 20
double COST::practical_bkz_cost_dd(int d,int beta,int jump){
    // int extra_dim4free = 12;
    // beta = beta + extra_dim4free;
    // f = f + extra_dim4free;
    int f = 0;
    // int extra_dim4free =12;
    pair<double,double> k;
    // bool sieve;
    if(beta < 50){
        // sieve = false;
        f = 0;
        return practical_bkz_cost_qd(d, beta, jump);
    }
    else{
        // sieve = true;  
        // f = get_f_for_pnjbkz(params, beta);
        f = get_f_for_pnjbkz(params, beta);
    }
    if(beta-f<=60)
        return practical_bkz_cost_qd(d, beta, jump);
    else
        k = get_k1_k2_pnj_dd(beta-f); // threads = 20
    double k1 = k.first, k2 = k.second;
    double c3= 0.018, c4 = -2.24;

    // if(beta - f <= 60)
    //     return (k1*(beta-f)+k2) - log2(jump);
    // else
    //if(c3*d+c4 > 1)
    if(c3*d+c4>=1)
        return (k1*(beta-f)+k2) + log2(c3*d+c4) - log2(jump);
    else
        return (k1*(beta-f)+k2) - log2(jump);
}



//get pump cost in threads = 20
pair<double,double> COST::practical_pump_cost_dd(int beta){
    //make sure not use the enum cost 
    // int f = get_f_for_pump(params, beta);
    // f = dims4free(beta);
    int beta_prime = beta;
    double secs, bits;
    pair<double,double> k = get_k1_k2_pump_dd(beta_prime); // threads = 20
    double k1 = k.first, k2 = k.second;
    // k = (1/71.)*((1.33)**(beta/10.));
    
    secs = k1*((double) beta_prime)+k2; 

    //unit: GB
    if( beta_prime <= 56)
    	bits = 2.0311;
    else if( beta_prime >=57 and  beta_prime <=63)
    	bits = 8e-5 *  pow(beta_prime,2) - 0.0083*beta_prime + 2.2555;
    else if( beta_prime >= 64 and beta_prime <= 94)
    	bits = 2.195202e-6 * pow(beta_prime,4) - 6.297613e-4 * pow(beta_prime,3) +6.803540e-2 * pow(beta_prime,2) - 3.274476 * beta_prime + 61.31963;
    else if( beta_prime >= 95)
    	bits = 2.0311 + pow(2,0.1992* beta_prime  - 18.714);

    //unit: log2(bit)
    bits = log2(bits * pow(2,33));
    return make_pair(secs,bits); //n_expected = beta -f , beta = d-llb
}






// //triple gpu time cost test in the artical G6K-GPU
// double pump_cost_triple_gpu(int beta){
//     double k1 = 0.367, k2 = -37.15;
//     return k1*((double)beta)+k2;
// }

pair<double,double> COST::sieve_cost(int beta,int cost_model){
    if(cost_model == 1)
        return make_pair(theo_pump_cost(beta).first - log2(COST::C),  theo_pump_cost(beta).second) ;
    else if(cost_model == 2){
        return make_pair(practical_pump_cost_qd(beta).first - log2(COST::C), practical_pump_cost_qd(beta).second);
    }
    else if(cost_model == 3){
        return make_pair(practical_pump_cost_dd(beta).first - log2(COST::C), practical_pump_cost_dd(beta).second);
    }
    return make_pair(params.max_num,params.max_num);
}
    
pair<double,double> COST::pump_cost(int beta,int cost_model){
    if(cost_model == 1)
        return theo_pump_cost(beta);
    else if(cost_model == 2){
        return practical_pump_cost_qd(beta);
    }
    else if(cost_model==3){
        return practical_pump_cost_dd(beta);
    }
    return make_pair(params.max_num,params.max_num);
}

pair<double,double> COST::bkz_cost(int d, int beta,int J,int cost_model){
    if(cost_model == 1)
        return theo_bkz_cost(d, beta, J);
    else if(cost_model == 2){
        return make_pair(practical_bkz_cost_qd(d,beta,J), practical_pump_cost_qd(beta).second);
    }
    else if(cost_model == 3){
        return make_pair(practical_bkz_cost_dd(d,beta,J), practical_pump_cost_dd(beta).second);
    }
    return make_pair(params.max_num,params.max_num);
}    




