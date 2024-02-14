from strategy_gen.strategy_gen import lwechal_simulation_gsa, lwechal_simulation_actual_l
from g6k.utils.util import load_lwe_instance,load_lwe_challenge
from g6k.utils.lwe_estimation import gsa_params, primal_lattice_basis
from g6k.siever import Siever
from math import log, log2
from fpylll.tools.quality import basis_quality


def strategy_simulation(n,alpha,S,load_lwe = "lwe_instance", simulation="gsa",m=None,float_type = None):
    if(load_lwe == "lwe_instance"):
        A, c, q = load_lwe_instance(n=n, alpha=alpha)
    if(load_lwe == "lwe_challenge" or load_lwe is None):
        A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    
    print("-------------------------")
    print("Primal attack, LWE challenge n=%d, alpha=%.4f" % (n, alpha))

    if m is None:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True, samples = A.nrows)
            (b, s, m) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    else:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True)
            (b, s, _) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    print()

    # no use in having a very small b
    b = max(b, s-65)

    target_norm = 1.5 * (alpha*q)**2 * m + 1
    # target_norm = max( target_norm, 0.98 * full_gh)

    
    # B_=load_lwe_challenge_mid(n=n, alpha=alpha)
    # if B_ is not None:
    #     B = B_
    # else: 
    #     B = primal_lattice_basis(A, c, q, m=m)
    B = primal_lattice_basis(A, c, q, m=m) #debug
    
    params = None

    g6k = Siever(B, params,float_type=float_type)
    print("GSO precision: ", g6k.M.float_type)
    print("||b_1|| = %d, target_norm = %d"  %(g6k.M.get_r(0, 0), target_norm))

    d = g6k.full_n

    g6k.lll(0, g6k.full_n)
    g6k.update_gso(0,d)
    slope = basis_quality(g6k.M)["/"]
    sigma = alpha * q
    dvol = g6k.M.get_log_det(0,d)/2. - log(sigma)*d
    print("Intial Slope = %.5f, dim = %d, dvol = %3.13f\n" %(slope, d, dvol))
    
    log2_rr = [round((log2(g6k.M.get_r(i,i))/2.) - (log2(sigma)),5) for i in range(d)]
    
    
    if(simulation=="actual_l"):
        lwechal_simulation_actual_l(log2_rr, S,len(S),float_type = float_type)
    if(simulation=="gsa"):
        lwechal_simulation_gsa(d, dvol, S,len(S),float_type = float_type)
        
        

n = 40
alpha = 0.030
float_type = "qd"
S = [(73, 8, 1), (89, 9, 1), (117, 10, 1), (119, 10, 1)]
strategy_simulation(n,alpha,S,load_lwe="lwe_challenge", simulation="actual_l",float_type = float_type)


n = 40
alpha = 0.025
float_type = "qd"
S =  [(89, 9, 1), (114, 10, 1)]
strategy_simulation(n,alpha,S,load_lwe="lwe_challenge", simulation="actual_l",float_type = float_type)


n = 45
alpha = 0.020
float_type = "qd"
S =  [(73, 8, 1), (90, 9, 1), (117, 10, 1)]
strategy_simulation(n,alpha,S,load_lwe="lwe_challenge", simulation="actual_l",float_type = float_type)


 
 
n = 50
alpha = 0.015
float_type = "qd"
S =  [(73, 8, 1), (90, 9, 1), (114, 10, 1)]
strategy_simulation(n,alpha,S,load_lwe="lwe_challenge", simulation="actual_l",float_type = float_type)
