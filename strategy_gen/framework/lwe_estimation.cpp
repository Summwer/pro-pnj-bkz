#include "lwe_estimation.h"

double log_gh_svp(double d, double delta_bkz, int svp_dim, int n, int q){
    /*
    Calculates the log of the Gaussian heuristic of the context in which
    SVP will be ran to try and discover the projected embedded error.

    The volume component of the Gaussian heuristic (in particular the lengths
    of the appropriate Gram--Schmidt vectors) is estimated using the GSA
    [Schnorr03] with the multiplicative factor = delta_bkz ** -2.

    NB, here we use the exact volume of an n dimensional sphere to calculate
    the ``ball_part`` rather than the usual approximation in the Gaussian
    heuristic.

    :param d: the dimension of the embedding lattice = n + m + 1
    :param delta_bkz: the root Hermite factor given by the BKZ reduction
    :param svp_dim: the dimension of the SVP call in context [d-svp_dim:d]
    :param n: the dimension of the LWE secret
    :param q: the modulus of the LWE instance

    */

    double svp_dim_ = double(svp_dim);
    double ball_part = ((1./svp_dim_)*lgamma((svp_dim_/2.)+1))-(.5*log(M_PI));
    double vol_part = ((1./d)*(d-(double)n-1)*log((double)q))+((svp_dim_-d)*log(delta_bkz));
    return ball_part + vol_part;
}


vector<tuple<int,int,int>> decoupler(LWEchal* lwechal){
    /*
    Creates valid (bkz_dim, svp_dim, d) triples, as determined by
    ``primal_parameters`` and determines which succeed in the recovery of the
    embedded error.

    :param decouple: if True the BKZ dimension and SVP dimension may differ
    :param n: the dimension of the LWE secret
    :param samples: maximum number of LWE samples to use for the embedding
        lattice. ``None`` means ``5*n``
    :param q: the modulus of the LWE instance
    :param stddev: the standard deviation of the distribution from which the
        error vector components were uniformly and indepedently drawn
    :param d: find best parameters for dimension ``d`` embedding lattice
    */
    int n = lwechal->n, q = lwechal->q;
    double alpha = lwechal->alpha;
    vector<tuple<int,int,int>> params = {};

    double stddev = alpha*(double)q;
    
    for(int m = n; m < min(5*n+1, lwechal->m+1); m++){
        int beta_bound = min(m+1, 400+default_dim4free_fun(400));
        int svp_bound = min(m+1, 400);
        for( int bkz_block_size = 40; bkz_block_size < beta_bound; bkz_block_size++){
            double delta_0 = compute_delta(bkz_block_size).get_d();

            for(int svp_dim = 40; svp_dim < svp_bound; svp_dim ++){
                double d = double(m + 1);
                FP_NR<FT> rhs = log_gh_svp(d, delta_0, svp_dim, n, q);
                if(rhs - log(stddev) - log(svp_dim)/2. >= 0)
                    params.insert(params.end(),{bkz_block_size, svp_dim, m+1});
            }
        }
    }

    return params;
}



tuple<int,int,int> find_min_complexity(vector<tuple<int,int,int>> params){
    /*
    For each valid and solving triple (bkz_dim, svp_dim, d) determines an
    approximate (!) cost and minimises.

    :param params: a list of all solving (bkz_dim, svp_dim, d) triples

    */

    double min_cost = -1;
    tuple<int,int,int> min_cost_param;
    double expo = .349;

    for(int i = 0; i < int(params.size()); i++) {
        tuple<int,int,int> param  = params[i];
   
        double bkz_block_size = (double) (get<0>(param) - default_dim4free_fun(get<0>(param)));
        double svp_dim = (double) (get<1>(param) - default_dim4free_fun(get<1>(param)));
        double d = (double) get<2>(param);

        double bkz_cost = 2 * d * pow(2, (expo * bkz_block_size));
        double finisher_svp_cost = pow(2, ((expo * svp_dim)));
        double new_cost = bkz_cost + finisher_svp_cost;

        if(min_cost < 0 or new_cost < min_cost){
            min_cost = new_cost;
            min_cost_param = param;
        }
    }
    return min_cost_param;
}



tuple<int,int,int> gsa_params(LWEchal* lwechal){
    /*
    Finds winning parameters (a BKZ reduction dimension and a final SVP call
    dimension) for a given Darmstadt LWE instance (n, alpha).

    lwechal includes the following parameters:

    :param n: the dimension of the LWE secret
    :param m: number LWE samples
    :param alpha: the noise rate of the LWE instance
    :param q: the modulus of the LWE instance. ``None`` means determine by
        reloading the challenge
    :param matrix: LWE samples
    */    
    vector<tuple<int,int,int>> params = decoupler(lwechal);
    // print_vector(params, 0, params.size());

    tuple<int,int,int> min_cost_param = find_min_complexity(params);

    
    return min_cost_param;
}



void primal_lattice_basis(LWEchal* lwechal){
    /*
    Construct primal lattice basis for LWE challenge
    ``(A,c)`` defined modulo ``q``.

    :param A: LWE matrix
    :param c: LWE vector
    :param q: integer modulus
    :param m: number of samples to use (``None`` means all)

    */
    
    int m = lwechal->m, n = lwechal->n, q = lwechal->q, d = lwechal->dim;
    ZZ_mat<ZT> A = lwechal->A;
    ZZ_mat<ZT> B;
    B.resize(m+n+1, m+1);
    for(int i=0; i < m; i++){
        for(int j=0; j < n; j++){
            B[j][i] = A[i][j];
        }
        B[i+n][i] = q;
        B[m+n][i] = lwechal->c[i];
    }
    B[m+n][m] = 1;
    // ostringstream os;
    // B.print(os);
	// cout<<os.str()<<endl;
    lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);

    // B.print(os);
	// cout<<os.str()<<endl;

    //B[:n] = O
    lwechal->B.resize(d,d);
    for(int i = 0; i < d; i++){
        for(int j = 0; j < d; j++)
            lwechal->B[i][j] = B[n+i][j];
    }
    // lwechal->B[0].print(os);
	// cout<<os.str()<<endl;
    
    // return B;
}


LWEchal* gen_lwechal_instance(int n, double alpha){
    LWEchal* lwechal = load_lwe_challenge(n, alpha);
    int q = lwechal->q, m = lwechal->m;
    
    printf("-------------------------\n");
    printf("Primal attack, TU LWE challenge n=%d, m=%d, alpha=%.3f, q=%d \n" , n, m, alpha, q);

    tuple<int,int,int> min_cost_param = gsa_params(lwechal);
    int b = get<0>(min_cost_param), s = get<1>(min_cost_param);
    m = get<2>(min_cost_param);
    int d = m  + 1;
    lwechal->m = m;
    lwechal->dim = d;
   
    printf("Chose %d samples. Predict solution at bkz-%d + svp-%d. \n", m, b, s);
    
    

    primal_lattice_basis(lwechal);

    // ostringstream os;
    // lwechal->B[0].print(os);
	// cout<<os.str()<<endl;

    double  sigma = alpha * q; 

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
    lwechal->dvol = M.get_log_det(0,d)/2. - log(sigma)*d;

    printf("dim = %d, dvol = %3.11f\n\n", lwechal->dim, lwechal->dvol.get_d());
   
    return lwechal;
}


// int main(){
//     int n = 40;
//     double alpha = 0.035;
//     gen_lwechal_instance(n, alpha);
// }