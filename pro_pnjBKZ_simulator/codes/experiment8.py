from g6k.utils.util import load_svpchallenge_and_randomize
from g6k.utils.util import load_lwe_challenge, load_matrix_file,load_lwe_challenge_mid
from g6k.utils.lwe_estimation_old import gsa_params, primal_lattice_basis
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour,dim4free_wrapper,default_dim4free_fun
from fpylll import BKZ as fplll_bkz
# from fpylll.util import gaussian_heuristic
from fpylll.algorithms.bkz2 import BKZReduction
from g6k.siever import Siever
from g6k.siever_params import SieverParams
from fpylll.tools.quality import basis_quality
import time
from g6k.algorithms.pump import pump
from math import log,pi,e,lgamma,exp
from scipy.special import chdtr 
import sys

def ball_log_vol(n):
    """
    Return volume of `n`-dimensional unit ball

    :param n: dimension

    """
    return (n/2.) * log(pi) - lgamma(n/2. + 1)

def theo_dim4free_fun1(blocksize):
    """
    Theoretical Dimension-for-free function 1 without e in [Duc18]
    """

    return int(blocksize*log(4/3.)/log(blocksize/2./pi)) 


def theo_dim4free_fun2(blocksize):
    """
    Theoretical Dimension-for-free function 2 with e in [Duc18]
    """

    return int(blocksize*log(4/3.)/log(blocksize/2./pi/e)) 

#calculate gh
def gaussian_heuristic(log_rr):
    """
    Return squared norm of shortest vector as predicted by the Gaussian heuristic.

    :param log_rr: vector of log(squared) Gram-Schmidt norms

    """
    n = len(list(log_rr))
    log_vol = sum([x for x in log_rr])
    log_gh =  1./n * (log_vol - 2 * ball_log_vol(n))
    return exp(log_gh)

def pump_estimation2(log_rr,q, alpha, succ_prob = 0.7, ebeta = 50, goal_margin=1.5):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(log_rr)

    for beta in range(50,d):
        
        GH = gaussian_heuristic(log_rr[d-beta:])
        
        length=(GH/(sigma**2))
        p = chdtr(beta, length) #integral_chi_squared distribution (n,x)
        # print(beta,p)
        if p > succ_prob:
            #print(beta)
            break

    llb = d - beta
    # while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
    #     llb -= 1
    #     if llb < 0:
    #         break
    llb = max(0, llb)
    # f = max(d-llb-beta, 0)
    f = dim4free_wrapper(default_dim4free_fun, beta)

    #pump time
    # log_pre_pump_time = get_log_pre_pump_time(beta - f) 

    return llb,d-llb,f





#test correctness of pump_estimation2
def experiment8(n,alpha,goal_margin,tours,succ_prob,max_sieve_dim,pre_blocksizes,threads,gpus):
    print("Experiment8. Start Pump Estimation Test " )
    file_name = "experiment8_gs_lengths.txt"
    f1 = open(file_name,"w")
    params = SieverParams(threads=threads,gpus=gpus)
    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    print("-------------------------")
    print("Primal attack, LWE challenge n=%d, alpha=%.4f" % (n, alpha))

    try:
        min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        samples=A.nrows, decouple=True)
        (b, s, m) = min_cost_param
    except TypeError:
        raise TypeError("No winning parameters.")

    print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    print()

    d = m + 1

    target_norm = goal_margin * (alpha*q)**2 * m + 1

    # B = primal_lattice_basis(A, c, q, m=m)
    B_=load_lwe_challenge_mid(n=n, alpha=alpha)
    B = B_
    g6k = Siever(B, params)
    slope = basis_quality(g6k.M)["/"]
    print("Intial Slope = %.5f" % slope)
    pump_params = {"down_sieve": True}
    T0 = time.time()
    for blocksize in pre_blocksizes:
        T0_BKZ= time.time()
        if blocksize < 61:
            print("Starting a fpylll BKZ-%d tour. " % (blocksize), end=' ')
            sys.stdout.flush()

            bkz = BKZReduction(g6k.M)
            par = fplll_bkz.Param(blocksize,strategies=fplll_bkz.DEFAULT_STRATEGY,max_loops=1)
            bkz(par)
        else:
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)

        slope = basis_quality(g6k.M)["/"]
        fmt = "slope: %.5f, BKZ_cost: %.3f sec, walltime: %.3f sec"
        print(fmt % (slope, time.time()-T0_BKZ, time.time() - T0))

        g6k.lll(0, g6k.full_n)

        if g6k.M.get_r(0, 0) <= target_norm:
            print("Finished! TT=%.2f sec" % (time.time() - T0))
            print(g6k.M.B[0])    
            log_rr = [log(g6k.M.get_r(i, i)) for i in range(d)]
            f1.write(str([g6k.M.get_r(i, i) for i in range(d)]))
            f1.write('\n')
            f1.close()
            return


    

    log_rr = [log(g6k.M.get_r(i, i)) for i in range(d)]
    f1.write(str([g6k.M.get_r(i, i) for i in range(d)]))
    f1.write('\n')
    llb,beta,f = pump_estimation2(log_rr,q, alpha,succ_prob)
    T0_pump = time.time()
    print("n_max = %d." %(beta - f))
    # print(llb)
    if beta - f < max_sieve_dim:
        for _ in range(tours):
            print()
            print( " %d/%d. Starting svp pump_{%d, %d, %d}" % (_+1,tours,llb, d-llb, f) ) # noqa
            sys.stdout.flush()
                
            pump(g6k, dummy_tracer, llb, d-llb, f,verbose=True,down_sieve=True)#,goal_r0=target_norm * (d - llb)/(1.*d))       
            print()
            f1.write(str([g6k.M.get_r(i, i) for i in range(d)]))
            f1.write('\n')
            T_pump = time.time() - T0_pump
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, T_pump = %.3f sec, walltime: %.3f sec"
            print(fmt % (slope, T_pump,time.time()-T0))

            g6k.lll(0, g6k.full_n)

            if g6k.M.get_r(0, 0) <= target_norm:
                print("Finished! TT=%.2f sec" % (time.time() - T0))
                print(g6k.M.B[0])    
                f1.close()
                return
    
        f1.close()
        raise ValueError("No solution found.")
    else:
        raise ValueError("Exceed max sieve dim = %d." %max_sieve_dim)

      

    
n,alpha,goal_margin,tours,succ_prob,max_sieve_dim,pre_blocksizes,threads,gpus = 40,0.035,1.5,10,0.95,150,[],40,2
experiment8(n,alpha,goal_margin,tours,succ_prob,max_sieve_dim,pre_blocksizes,threads,gpus)