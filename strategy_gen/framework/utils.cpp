#include "utils.h"


// return all keys in map
vector<int> KeySet(map<int, double> mp)
{
    vector<int> keys;
    for(map<int, double>::iterator it = mp.begin(); it != mp.end(); ++it){
        keys.push_back(it->first);
    }
    return keys;
}

// return all values in map
vector<double> ValueSet(map<int, double> mp)
{
    vector<double> values;

    for(map<int, double>::iterator it = mp.begin(); it != mp.end(); ++it){
        // values.push_back(rational_cast<double>(it->second));
        values.push_back(it->second);
    }
    return values;
}


// Generate samples list such that number of each sample obeys a specific prob * sample_num. 
vector<int> expand_samples(vector<int> samples,vector<double> probs,int sample_num){
    vector<int> expand_samples_list;
    int index = 0;
    int upper_index = 0;

    for(int i=0; i < int(samples.size()); i++){
        if(i < int(samples.size())-1)
            upper_index += sample_num*probs[i];
        else
            upper_index = sample_num;

        expand_samples_list.insert(expand_samples_list.begin()+index, upper_index-index, samples[i]);
        index = upper_index;
    }
    return expand_samples_list;
}


void print_map(map<int, double>  mp){
    cout<<"{\t";
    for (map<int, double> ::iterator it = mp.begin();
		it != mp.end(); it++) {
            cout<<it-> first<<","<<it->second<<"\t";
	}
    cout<<"}"<<endl;
}

void print_map(map<int, FP_NR<FT>>  mp){
    cout<<"{\t";
    for (map<int, FP_NR<FT>> ::iterator it = mp.begin();
		it != mp.end(); it++) {
            cout<<it-> first<<","<<it->second<<"\t";
	}
    cout<<"}"<<endl;
}

void print_vector(vector<double> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end-1; i++) {
        cout << v[i] << ", ";
    }
    cout<<v[index_end-1];
    cout<<"]"<<endl;
}

void print_vector(vector<int> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end-1; i++) {
        cout << v[i] << ", ";
    }
    cout<<v[index_end-1];
    cout<<"]"<<endl;
}

void print_vector(vector<Z_NR<ZT>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end-1; i++) {
        cout << v[i] << ", ";
    }
    cout<<v[index_end-1];
    cout<<"]"<<endl;
}

void print_vector(vector<FP_NR<FT>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end-1; i++) {
        cout << v[i] << " ";
    }
    cout<<v[index_end-1];
    cout<<"]"<<endl;
}

void print_vector(vector<pair<int,int>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << (v[i].first) << ", " << v[i].second << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<pair<double,double>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << (v[i].first) << ", " << v[i].second << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<tuple<double,double,double>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << get<0>(v[i]) << ", " << get<1>(v[i]) << ", " << get<2>(v[i])  << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<tuple<int,int,int>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << get<0>(v[i]) << ", " << get<1>(v[i]) << ", " << get<2>(v[i])  << ") ";
    }
    cout<<"]"<<endl;
}

/*
void print_matrix(ZZ_mat<ZT> matrix){
    cout<< "[";
    for(auto it = matrix.begin(); it != matrix.end(); ++it) {
        cout << *it << " ";
    }
    cout<<"]"<<endl;
}
*/


void print_matrix(vector<vector<Z_NR<ZT>>> matrix){
    int row = matrix.size(), col = matrix[0].size();

    printf("[");
    for(int i = 0; i < row; i++){
        print_vector(matrix[i],0,col);
    }
    printf("]\n");
}

void printf_input(int d, FP_NR<FT> dvol)
{
    printf("\033[0m\033[1;31m dim = %3d, dvol = %3.2f \033[0m\n", d, dvol.get_d());
}


void printf_red(const char *s)
{
    printf("\033[0m\033[1;31m%s\033[0m\n", s);
}

pair<double,double> average_variance(std::map<int,double> D){
    double mu = 0.,s = 0.;

    int v;
    double p;
    
    for(auto it : D){
        v = it.first;
        p = it.second;
        mu += v * p;
        s += v * v * p;
    }

    s -= mu * mu;

    return make_pair(mu,s);
}


int draw_from_distribution(std::map<int,double> D, int sample_num){
    /*draw an element from the distribution D
    ,D, distribution in a dictionnary form
    */

    // double rational_cast<double>
    vector<int> samples = KeySet(D);
    vector<double> probs= ValueSet(D);
    vector<int> expand_samples_list = expand_samples(samples,probs,sample_num);

    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0,sample_num-1);
    

    return expand_samples_list[distrib(gen)];

    //nc::random::choice(a, 3); //Corresponding to numpy.random.choice function in numpy.

    //X = np_random_choice([key for key in D.keys()],
    //                     1, replace=True,
    //                     p=[float(prob) for prob in D.values()])
    //return X[0]
}

//It will occur a double free error!!!
// void gen_samples(ZZ_mat<ZT> &matrix, int m, int n, int q){
//     /*
//     m, number of samples
//     n, dimension for each sample
//     */
    
    
//     Z_NR<ZT> q2;
//     q2=q;
    
//     for (int i = 0; i < m; i++)
//         for (int j = 0; j < n; j++)
//             matrix[i][j].randm(q2);
// }

void build_LWE_lattice(ZZ_mat<ZT> &matrix, ZZ_mat<ZT> A, int q){
    /*
    Builds a (m+n)*(m+n) matrix of the form
    q*I, 0
    -A^T, I
    It corresponds to the LWE matrix
    ,A, a m*n matrix
    ,q, integer
    */

    int m = A.get_rows(), n = A.get_cols(), d = m+n;

    //draw matrix B and define the lattice, B = [[qI,0],[-A^T,I]], A is lattice samples. Each col is one group of sample.
    matrix.gen_zero(d,d);
    for (int i = 0; i < m; i++)
        matrix[i][i] = q;

    for (int i = m; i < d; i++)
        matrix[i][i] = 1;

    for (int i = m; i < d; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j].mul_si(A[j][i-m],-1);
   
}



void kannan_embedding(ZZ_mat<ZT> &matrix, vector<Z_NR<ZT>> target, int factor){
    /*
    Builds a (m+n)*(m+n) matrix of the form
    q*I, 0,  0
    -A^T, I, 0
     b ,  0, 1
    It corresponds to the LWE matrix
    ,A, a m*n matrix
    ,q, integer
    */

    int d = matrix.get_rows()+1;
    matrix.resize(d,d);
    
    for (int i = 0; i < d-1; i++){
        matrix[d-1][i] = target[i];
    }
    
    matrix[d-1][d-1] = factor;
    
}


FP_NR<FT> compute_delta( int beta){
    /*
    Computes delta from the block size k. Interpolation from the following
    data table,
    Source , https,//bitbucket.org/malb/lwe-estimator/
    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
    tor.py?at=master&fileviewer=file-view-default
    ,k, integer
    estimator.py table,
    */

    map<int,FP_NR<FT>> small= {{0,1e20}, {1, 1e20}, {2, 1.021900}, {3, 1.020807}, {4, 1.019713}, {5, 1.018620}, {6, 1.018128}, {7, 1.017636}, {8, 1.017144}, {9, 1.016652}, {10, 1.016160}, {11, 1.015898}, {12, 1.015636}, {13, 1.015374}, {14, 1.015112}, {15, 1.014850}, {16, 1.014720}, {17, 1.014590}, {18, 1.014460}, {19, 1.014330}, {20, 1.014200}, {21, 1.014044}, {22, 1.013888}, {23, 1.013732}, {24, 1.013576}, {25, 1.013420}, {26, 1.013383}, {27, 1.013347}, {28, 1.013310}, {29, 1.013253}, {30, 1.013197}, {31, 1.013140}, {32, 1.013084}, {33, 1.013027}, {34, 1.012970}, {35, 1.012914}, {36, 1.012857}, {37, 1.012801}, {38, 1.012744}, {39, 1.012687}, {40, 1.012631}, {41, 1.012574}, {42, 1.012518}, {43, 1.012461}, {44, 1.012404}, {45, 1.012348}, {46, 1.012291}, {47, 1.012235}, {48, 1.012178}, {49, 1.012121}, {50, 1.012065}};

    if(beta <=50)
        return small[beta];
    else{
        return pow((double)beta/(2.*M_PI*exp(1)) * pow((M_PI*(double)beta),(1./(double)beta)), (1./(2*((double)beta-1.))));
    }
}



//calculate slope in G6K
double get_current_slope(vector<double> l, int start, int end){
    /*
    Return current slope of 2*ln(gs-lengths)

    :param log_rr: vector of log2(gs-lengths) Gram-Schmidt norms
    :param start_row: start row
    :param stopr_row: stop row
    
    */

    for(int i = 0; i < int(l.size()); i++){
        l[i] *= 2*log(2.);
    }
    
    int n = end - start;
    double i_mean = (n - 1) * 0.5 + start;
    double x_mean =  accumulate(l.begin(), l.end(), 0.)/n;
    double v1 = 0.0, v2 = 0.0;
    for(int i = start; i < end; i++){
        v1 += (i - i_mean) * (l[i] - x_mean);
        v2 += (i - i_mean) * (i - i_mean);
     }
    return v1 / v2 ; //round(v1 / v2 * 1e6)/1e6; 
}



FP_NR<FT> bkzgsa_gso_len( FP_NR<FT> logvol, int i, int d,int beta){
    FP_NR<FT> gso_len;
    //  FP_NR<FT> delta, 
    //if(delta >= 1000.)
    FP_NR<FT> delta = compute_delta(beta);
    gso_len.pow_si(delta,(d-1-2*i));
    FP_NR<FT> tmp;
    tmp.exponential(logvol.get_d()/(double)(d));
    gso_len.mul(gso_len,tmp);
    // cout<<"delta:"<<delta<<", "<<"exp(logvol / d)): "<<tmp<<endl;

    return gso_len;
}



// double gaussian_heuristic(vector<FP_NR<FT>> l, int index_start){
//     int n = l.size();

//     double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1))/(log(2.));
//     FP_NR<FT> log_lattice_square_vol = 0.;
//     FP_NR<FT> gh;
//     for (int i = index_start; i < n; ++i)
//     {
//         log_lattice_square_vol += l[i];
//     }

//     //gh.exponential(1);
//     double gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));
//     // cout<<log_lattice_square_vol.get_d()<<","<<log_ball_square_vol<<endl;
//     return gh.get_d();

// }

//Input: log2(square gs-lengths)
//return: gh**2
double gaussian_heuristic_log2(vector<FP_NR<FT>> l, int index_start, int index_end){
    //int d = l.size();
    int n = index_end - index_start;
    assert(n>0);
    double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1));
    double log_lattice_square_vol = 0., gh;
    for (int i = index_start; i < index_end; ++i)
    {
        log_lattice_square_vol += (2*l[i].get_d()*log(2.));
    }
    // cout<<log_lattice_square_vol<<","<<log_ball_square_vol  <<endl;

    //FP_NR<FT> gh;
    //gh.exponential((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));
    gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));

    return gh;

}

//Input: log2(gs-lengths)
double gaussian_heuristic_log2(vector<double> l, int index_start, int index_end){
    // int d = l.size();
    int n = index_end - index_start;
    assert(n>0);

    double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1));
    double log_lattice_square_vol = 0., gh;
    for (int i = index_start; i < index_end; ++i)
    {
        log_lattice_square_vol += (2*l[i]*log(2.));
    }
    gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));

    return gh;


}


//Generate input simulated-gs-lengths
vector<double> gen_simulated_gso(int d, FP_NR<FT> logvol){
    FP_NR<FT> delta, gso_len;
    // delta = compute_delta(2);
    
    vector<FP_NR<FT>> l;
    l.resize(d);
    for(int i = 0; i<d; i++){
        gso_len = bkzgsa_gso_len(logvol, i, d, 2); // delta, beta = 2
        l[i].log(gso_len);
        l[i].div(l[i],log(2));
    }
    // print_vector(l);

    vector<double> l_;
    l_.resize(d);
    for(int i = 0; i<d; i++){
        l_[i] = l[i].get_d();
    }

    
    return l_;
    
}



//=theo_dim4free_fun2
int dims4free(int beta){
    if(beta < 40)
        return 0;
    return max(0,int(ceil((double)beta * log(4./3.) / log((double)beta/(2.*M_PI*exp(1.))))));
    // return max(0,int(((double)beta * log(4./3.) / log((double)beta/(2.*M_PI*exp(1.))))));
}


int wrapper_default_dim4free_fun(int beta){
    if(beta < 40)
        return 0;
    
    // //For lwe-est
    // // int f = int(11.5 + 0.075*beta);
    // // return min(int(((double)beta - 40)/2.), f);

    // //For nist-est
    // int f = int(ceil(11.5 + 0.075*beta));
    int f = int(11.5 + 0.075*beta);
    return int(min(((double)beta - 40)/2., (double) f));
}


int default_dim4free_fun(int beta){
    if(beta < 40)
        return 0;
    
    int f = int(11.5 + 0.075*beta);
    return f;
}

int theo_dim4free_fun1(int beta){
    if(beta < 40)
        return 0;
    int f = max(0,int(ceil((double)beta * log(4./3.) / log((double)beta/(2.*M_PI)))));
    // return int(min(((double)beta - 40)/2., (double) f));
    return f; 
}


int theo_dim4free_fun2(int beta){
    if(beta < 40)
        return 0;
    // return min(int(ceil(((double)beta - 40)/2.)), dims4free(beta));
    int f = max(0, dims4free(beta));
    // return int(min((double)beta - 40)/2., (double) dims4free(beta));
    return f;
}




int get_beta_from_sieve_dim(int sieve_dim, int d, int choose_dims4f_fun){
    int f = 0;

    for(int beta = sieve_dim; beta < d; beta++){
        if(choose_dims4f_fun == 1 )
            f = theo_dim4free_fun1(beta);
        else if(choose_dims4f_fun == 2)
            f = theo_dim4free_fun2(beta);
        if(beta - f >= sieve_dim)
            return beta;
    }
    return d;
}


//d4f(B): pessimistic
int theo_dim4free1_in_B(vector<double> l){
    int d = l.size();
    double gh = gaussian_heuristic_log2(l,0,d);
    for(int f = d-1; f >= 0; f--){
        double ghf = gaussian_heuristic_log2(l,f,d);
        if(ghf * 4/3. >=  gh){
            // cout<<"f = "<<f<<endl;
            return f;
        }
    }
    return 0;
}


// //d4f(B): pessimistic
// int theo_dim4free1_in_B(vector<double> l, int beta){
//     int d = l.size();
//     double gh = gaussian_heuristic_log2(l,0,beta);
//     for(int f = beta-1; f >= 0; f--){
//         double ghf = gaussian_heuristic_log2(l,f,beta);
//         if(ghf * 4/3. >=  gh){
//             // cout<<"f = "<<f<<endl;
//             return f;
//         }
//     }
//     return 0;
// }


//d4f(B): postive 
int theo_dim4free2_in_B(vector<double> l){
    int d = l.size();
    double gh = gaussian_heuristic_log2(l,0,d);
    for(int f = d-1; f >= 0; f--){
        double ghf = gaussian_heuristic_log2(l,f,d);
        // cout<<f<<","<<ghf + log2(4/3.)<<", "<< log2((double) (d-f)/d) + gh<<endl;
        if(ghf * 4/3. >=  ((double) (d-f)/d) * gh){
            // cout<<"f = "<<f<<endl;
            return f;
        }
        // gh = gaussian_heuristic_log2(l,d-dsvp,d);
        // boost::math::chi_squared chisquare(d-f);
        // double psvp = boost::math::cdf(chisquare,ghf);
        // if(psvp>=0.999)
        //     return f;
    }
    return 0;
}


// //d4f(B): postive 
// int theo_dim4free2_in_B(vector<double> l, int beta){
//     int d = l.size();
//     int minf = beta;
//     for(int i = 0; i<d-beta; i++){
//         // cout<<f<<","<<ghf + log2(4/3.)<<", "<< log2((double) (d-f)/d) + gh<<endl;
//         bool has_f = false;
//         for(int f = beta-2; f >= 0; f--){
//             double gh = gaussian_heuristic_log2(l,i,beta+i);
//             double ghf = gaussian_heuristic_log2(l,f+i,beta+i);
//             if(minf > f and ghf >=  ((double) (beta-f)/beta) * gh){
//                 // cout<<"f = "<<f<<endl;
//                 // cout<<"beta = "<<beta<<", " << "f = "<<f<<endl;
//                 // return f;
//                 minf = f;
//                 has_f = true;
//                 break;
//             }
//         }
//         if(not has_f)
//             minf = 0;
//     }
//     return minf;
// }

int accs_2023_d4f(double slope){
    assert(slope < 0);
    int f = (int) floor(log(sqrt(4/3.))/(-1*slope/4.));
    assert(f>=0);
    return f;
}



int jump_upper_bound(Params params, int beta, vector<double> l){
    int jub = 0;
    double slope = get_current_slope(l,0,(int) l.size());
    if(beta<=50)
        return 1;
    // cout<<"params.compute_jub = "<<params.compute_jub<<endl;
    switch(params.compute_jub){
        case 1:
            jub = theo_dim4free1_in_B(l);
            break;
        case 2:
            jub = theo_dim4free2_in_B(l);
            // jub = theo_dim4free2_in_B(l,beta);
            // cout<<"jub = "<<jub<<endl;
            // cout<<"theo_f = "<< floor(get_f_for_pnjbkz(params,beta))<<endl; 
            break;
        case 3:
            jub = floor((double) get_f_for_pnjbkz(params,beta)/2.);       
            break;
        case 4:
            jub = min((double)get_f_for_pnjbkz(params,beta), ceil(0.1* beta));       
            break;
        case 5:
            jub = min(accs_2023_d4f(slope),theo_dim4free_fun1(beta));
            if(params.cost_model == 3)
                jub = min(jub, min((int) floor(params.jub_ratio * beta), (int) floor((double) get_f_for_pnjbkz(params,beta)/2.)));
            if(params.cost_model == 2)
                jub = min(jub,  (int) floor((double) get_f_for_pnjbkz(params,beta)/2.));
            // jub = min(jub, (int) floor((double) get_f_for_pnjbkz(params,beta)/2.));
            break;
        default:
            jub = 0;
            break;
    }
    // cerr<<"jub = "<<jub<<", beta = "<<beta<<endl;
    assert(jub < beta and jub >= 0);
    return jub;
    // return min(jub,(int) floor((double) get_f_for_pnjbkz(params,beta)/2.));

    // return min(jub, min((int) floor(0.1* beta), (int) floor((double) get_f_for_pnjbkz(params,beta)/2.)));
}


int get_f_for_pnjbkz(Params params, int beta){
    if(params.cost_model == 1){
        if(params.theo_pnjbkz_d4f == 1)
            return max(0,theo_dim4free_fun1(beta));
        if(params.theo_pnjbkz_d4f == 2)
            return max(0,theo_dim4free_fun2(beta));
        if(params.theo_pnjbkz_d4f == 3)
            return max(0,wrapper_default_dim4free_fun(beta));
    }
    if(params.cost_model >= 2){
        if(params.practical_pnjbkz_d4f == 1)
            return max(0,theo_dim4free_fun1(beta));
        if(params.practical_pnjbkz_d4f == 2)
            return max(0,theo_dim4free_fun2(beta));
        if(params.practical_pnjbkz_d4f == 3)
            return max(0,wrapper_default_dim4free_fun(beta));
    }
    return 0;
}


int get_f_for_pump(Params params, int beta, vector<double> l){
    if(params.cost_model == 1){
        if(params.theo_pump_d4f == 1)
            return max(0,theo_dim4free_fun1(beta));
        if(params.theo_pump_d4f == 2)
            return max(0,theo_dim4free_fun2(beta));
        if(params.theo_pump_d4f == 3)
            return max(0,wrapper_default_dim4free_fun(beta));
        if(params.theo_pump_d4f == 4)
            return max(0, theo_dim4free2_in_B(l));
    }
    if(params.cost_model >= 2){
        if(params.practical_pump_d4f == 1)
            return max(0,theo_dim4free_fun1(beta));
        if(params.practical_pump_d4f == 2)
            return max(0,theo_dim4free_fun2(beta));
        if(params.practical_pump_d4f == 3)
            return max(0,wrapper_default_dim4free_fun(beta));
        if(params.practical_pump_d4f == 4)
            return max(0, theo_dim4free2_in_B(l));
            // return max(0, theo_dim4free2_in_B(l,beta));
    }
    return 0;
}


int d4f_gap(int beta){
    int gap = 0;
    if(beta<40)
        gap = 0;
    else if(beta<76)
        gap = round((beta-40)/2.-theo_dim4free_fun1(beta));
    else if(beta>0)
        gap = round((11.5+0.075*beta)-theo_dim4free_fun1(beta));
        // return round(gap)
    else
        throw "Wrong: beta<=0!!";
    // else
    //     return round(wrapper_default_dim4free_fun(beta) - theo_dim4free_fun1(beta));
    if(gap<0)
        gap = 0;
    // assert(gap>=0);
    assert(gap < beta);
    return gap;
}

int d4f_gap_accs(int beta, double slope){
    int gap = 0;
    if(beta<40)
        gap = 0;
    else if(beta<75)
        gap = round((beta-40)/2.-min(theo_dim4free_fun1(beta),accs_2023_d4f(slope)));
    else if(beta>0)
        gap = round((11.5+0.075*beta)-min(theo_dim4free_fun1(beta),accs_2023_d4f(slope)-1));
        // return round(gap)
    else
        throw "Wrong: beta<=0!!";
        // else
        //     return round(wrapper_default_dim4free_fun(beta)-min(theo_dim4free_fun1(beta),accs_2023_d4f(slope)));
    // if(gap<0 or gap > beta){
    //     cout<<"slope = "<<slope<<endl;
    //     cout<<"theo_dim4free_fun1(beta) = "<<theo_dim4free_fun1(beta)<<endl;
    //     cout<<"accs_2023_d4f(slope) = "<<accs_2023_d4f(slope)<<endl;
    //     cout<<"gap = "<<gap<<", beta = "<<beta<<endl;
    // }
    if(gap<0)
        gap = 0;
    assert(gap < beta);
    return gap;
}
    


int get_beta_(Params params, int beta, int jump, int d, double slope){
    if(params.sim_d4f == 1){
        int f = get_f_for_pnjbkz(params, beta);
        if(jump <= 2)
            return beta;
        else if(jump >=3 && jump <=4){
            if((params.cost_model >= 2 and params.practical_pnjbkz_d4f == 3) or (params.cost_model == 1 and params.theo_pnjbkz_d4f == 3)){
                return get_beta_from_sieve_dim( beta-f,d,2);
            }else
                return beta;
        }
        else if(jump>=5){
            if((params.cost_model >= 2  and params.practical_pnjbkz_d4f == 1) or (params.cost_model == 1 and params.theo_pnjbkz_d4f == 1))
                return beta;
            else
                return get_beta_from_sieve_dim( beta-f,d,1);
        }
        return beta;
    }
    if(params.sim_d4f == 2){
        if(jump ==  1){
            //"Optimistic d4f value"
            return beta;
        }
        else if(jump/(double) beta <= 0.05){
            //"Theory conservative d4f value"
            // cout<<"beta = "<<beta<<", "<< "beta - d4f_gap(beta) = "<<beta - d4f_gap(beta)<<endl;
          
            return beta - d4f_gap(beta);
        }
        else if(jump/(double) beta > 0.05){
            //d4f_value = "ACCS d4f value"
            // cout<<"beta = "<<beta<<", "<< "beta - d4f_gap_accs = "<<beta - d4f_gap_accs(beta,slope)<<endl;
            return beta - d4f_gap_accs(beta,slope);
        }
    }
    return beta;
}



int factorial(int n) 
{
    if (n == 0)
       return 1;
    return n * factorial(n - 1);
}

//Implement centered_binomial distribution
int binomial(int x, int y){
    /*Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    */
    if(x<y)
        return 0;
    else
        return factorial(x)/factorial(y)/factorial(x-y); 
}



double centered_binomial_pdf(int k, int x){
    /* Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    */
    return (double) binomial(2 * k, x + k) / pow(2., 2 * k);
}

std::map<int,double>  build_centered_binomial_law(int k){
    /* Construct the binomial law as a dictionnary
    :param k: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    */
    std::map<int,double> D = {};
    for(int i=-k; i<=k; i++)
        D[i] = centered_binomial_pdf(k, i);
    return D;
}




void print_param_setting(Params params){
    // cout<<params.bssa_tradition<<endl;
    if(params.method == 1)
        printf("EnumBS Strategy Generation:");
    if(params.method == 2 and params.bssa_tradition == true)
        printf("BSSAv1 Strategy Generation:");
    if(params.method == 2 and params.bssa_tradition == false)
        printf("BSSAv2 Strategy Generation:");
    printf("v1. beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d,  G_prec = %e,  slope_prec = %e,  progressive_sieve = True, worst_case = %d, succ_prob = %1.3f, ", params.beta_start, params.gap, params.J, params.J_gap, params.cost_model, params.max_loop,params.enumbs_G_prec, params.enumbs_slope_prec,  params.worst_case, params.succ_prob);

    if(params.method == 1)
        printf("threads = %d, ", params.threads);
    if(params.enumbs_min_G)
        printf("Find expected minimal time cost strategy, ");
    else
        printf("Find a strategy below the maximal RAM = %f log2(bit), ", params.max_RAM);
    if(params.worst_case)
        printf("worst_case, ");
    else
        printf("average_case, ");
    
    if(params.cost_model == 1){
        printf("theo_pnjbkz_d4f = %d, ", params.theo_pnjbkz_d4f);
        // if(params.est_model !=3)
        printf("theo_pump_d4f = %d, ", params.theo_pump_d4f);
        // else
            // printf("theo_pump_d4f: compute ||pi_f(target_vector)||<= sqrt(4/3) GH(L_f)\n");
        printf("list_decoding = `%s`. \n\n", params.list_decoding.c_str());
    }
    if(params.cost_model >= 2){
        printf("practical_pnjbkz_d4f = %d, ", params.practical_pnjbkz_d4f);
        // if(params.est_model !=3)
        printf("practical_pump_d4f = %d, ", params.practical_pump_d4f);
        // else
            // printf("practical_pump_d4f: compute ||pi_f(target_vector)||<= sqrt(4/3) GH(L_f)\n");
    }
    if(params.sim_d4f == 1)
        cout<<"adjust beta_ as before 2023.1.31 (heuristic ajust), ";
    if(params.sim_d4f == 2)
        cout<<"adjust beta_ as after 2023.1.31 [WWW, asiaccs2023], ";
    cout<<"jump upper bound = ";
    switch(params.compute_jub){
        case 1:
            cout<<"pessimistic d4f(B)."<<endl; 
            break;
        case 2:
            cout<<"optimistic d4f(B)."<<endl;
            break;
        case 3:
            cout<<" d4f(beta)/2."<<endl;   

            break;
        case 4:
            cout<<"min(d4f(beta),0.1beta)."<<endl;   
            break;
        case 5:
            cout<<"min(asiaccs_d4f(beta,slope), pessitive d4f(beta))."<<endl;   
            break;
        default:
            cout<<"No setting"<<endl;
            break;
    }
    cout<<endl;
}