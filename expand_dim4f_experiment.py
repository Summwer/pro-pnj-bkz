


''' 
Goal of expand_dim4f_experiment:
Expand the range of dim4f value, narrow the dimension of sieve

'''
from g6k.algorithms.ducas18 import ducas18
from g6k.algorithms.bkz import default_dim4free_fun,dim4free_wrapper,pump_n_jump_bkz_tour
from g6k.algorithms.pump import pump
from g6k.utils.stats import dummy_tracer
from g6k.siever_params import SieverParams
from g6k.utils.util import load_svpchallenge_and_randomize,load_lwe_challenge,load_lwe_challenge_mid
from g6k.siever import Siever
from fpylll.tools.quality import basis_quality
from fpylll.util import gaussian_heuristic
from math import log,pi
import psutil,os
from g6k.siever import SaturationError
from g6k.utils.lwe_estimation_old import gsa_params, primal_lattice_basis
from scipy.special import chdtr 
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
import time


#保守f1
def theoretical_dim4free_fun1(blocksize):
    return int(blocksize*log(4/3.)/log(blocksize/(2*pi)))

#保守f2
def theoretical_dim4free_fun2(blocksize):
    return int(blocksize*log(4/3.)/log(blocksize/(2*pi*e)))


def load_lwe_lattice(n,alpha):
    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    print("-------------------------")
    print("Primal attack, LWE challenge n=%d, alpha=%.4f" % (n, alpha))
    min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,decouple=True)
    (b, s, m) = min_cost_param
    # B = primal_lattice_basis(A, c, q, m=m)
    B=load_lwe_challenge_mid(n=n, alpha=alpha)
    print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    print("")
    
    return B,q


def pump_estimation1(rr,q,alpha, goal_margin=1.5, ebeta = 50):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param rr: vector of squared Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    #goal_margin = 1 #1.5
    d = len(rr)
    # file.write(d)
    #compute n_expected
    m = d - 1
    stddev = q*float(alpha)
    target_norm = goal_margin * ((stddev**2) * m + 1 )

    for n_expected in range(2, d-2):
        # file.write(n_expected,d)
        x = (target_norm/goal_margin) * n_expected/(1.*d)
        if 4./3 * gaussian_heuristic(rr[d-n_expected:]) > x:
            break
    
    llb = d - ebeta
    # print(d,llb,ebeta)
    while gaussian_heuristic(rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
        llb -= 1
        if llb < 0:
            break
    llb = max(0, llb)
    f = max(d-llb-n_expected, 0)

    return llb,d-llb,f



def pump_estimation2(rr,q, alpha, succ_prob = 0.90, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(rr)

    for beta in range(50,d):
        GH = gaussian_heuristic(rr[d-beta:])
        length=( 4/3. * GH/(sigma**2)) #Condition: shortest vector projection in short vector set 
        # length=(GH/(sigma**2)) #Condition: shortest vector's projection is smaller than projected GH.
        p = chdtr(beta, length) #integral_chi_squared distribution (n,x)
        if p > succ_prob:
            #print(beta)
            break
        
    llb = d - beta
    target_norm =  goal_margin * ((sigma**2) * (d-1) + 1 )
    while gaussian_heuristic(rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
        llb -= 1
        if llb < 0:
            break
    llb = max(0, llb)
    
    return llb,d-llb,d-llb-beta


def test1(params,A,shrinkage,d,alg,dim4free,threads,saturation_ratio,start_n = 50):
    # A, bkz = load_svpchallenge_and_randomize(d)
    g6k = Siever(A, params)
    slope = basis_quality(g6k.M)["/"]
    gh = gaussian_heuristic([g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)])
    
    print("Test1. Start %s_sieve(%d,%d,%d), try to find the shortest vector whose norm less than %.3f GH=%d." %(alg,0,dim4free,d,shrinkage,gh))
    print("Intial Slope = %.5f." % (slope))
    g6k.lll(0, d)
    g6k.initialize_local(0, d - start_n, d)
    g6k.shrink_db(0)

    print("\r %3d: ↑%3d  , RAM cost:%.4f GB   " % (d-dim4free,g6k.r-g6k.l,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
    # with g6k.temp_params(dh_dim4free=0):
    T0 = time.time()
    try:
        g6k(alg=alg)
    except SaturationError:
        pass

    while g6k.l > 0:
        # Extend the lift context to the left
        g6k.extend_left(1)
        if g6k.l >= dim4free:
            print("\r %3d: ↑%3d  , RAM cost:%.4f GB   " % (d-dim4free,g6k.r-g6k.l,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
            try:
                g6k(alg=alg)
            except SaturationError:
                pass
    print()
    db = list(g6k.itervalues())
    print("T(sieve) = %.3f s. Finish sieve, db_size = %d" %(time.time()-T0,len(db)))
    print("------------------------------------------")
    
    print("Start check short vectors...")
    #check whether it exists the shortest vector in db
    i = 0
    found = 0
    short_vectors = {}
    min_v = []
    min_l = float("inf")
    for x in db:
        v = A.multiply_left(x)
        l = sum(v_**2 for v_ in v)
        if l < shrinkage * gh:
            # print()
            # print(l/gh, v)
            found += 1
            short_vectors[round(l/gh,3)] = v
            if abs(v[-1]) == 1:
                print(v)
                break
        print("\r %d/%d" %(i,len(db)), end= '')
        if l < min_l:
            min_l = l
            min_v = v
        i += 1
    print()
    if found ==0:
        print("Fail to find the shortest vector.")
    else:
        print("Find the satisfied short vector = %.3f GH: "%min(short_vectors.keys()))
    print(min_v)
    # print()
    print("==============================================\n\n")
    
#对照组，看直接一个pump能不能解出来短向量解。
def control_group(params,A,d,dim4free,alg,threads,saturation_ratio,start_n = 30):
   
    pump_params = {"down_sieve": True,"saturation_error":"ignore"}
    # A, bkz = load_svpchallenge_and_randomize(d)
    g6k = Siever(A, params)
    slope = basis_quality(g6k.M)["/"]
    gh = gaussian_heuristic([g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)])
    print("CG. Start pump(%s_sieve(%d,%d,%d)), GH=%d." %(alg,0,dim4free,d,gh))
    print("Intial Slope = %.5f." % (slope))
    T0 = time.time()
    pump(g6k, dummy_tracer, 0, d, dim4free, verbose=True,**pump_params)
    
    print()
    print("T(pump) = %.3f s. The shortest vector = %.3f GH:" %(time.time()-T0,round(g6k.M.get_r(0,0)/gh,3)))
    print(g6k.M.B[0])
    shrinkage = round(g6k.M.get_r(0,0)/gh,3) + 0.001
    # print()
    
    print("=================================================\n\n")
    return shrinkage
        
        
        
alg,threads,saturation_ratio,gpus = "gpu",48,0.2,2
# alg,threads,gpus,saturation_ratio,jump = "gpu",32,2,0.5,4
if gpus ==0:
    params = SieverParams(threads=threads,saturation_ratio=saturation_ratio)
else:
    params = SieverParams(threads=threads,gpus=gpus,saturation_ratio=saturation_ratio)
    
# dim4free = dim4free_wrapper(default_dim4free_fun,d)
# A, _ = load_svpchallenge_and_randomize(d, seed = 0)


# shrinkage = control_group(params,A,d,dim4free,alg,threads,saturation_ratio)
# print("\n\n\n")
# A, _ = load_svpchallenge_and_randomize(d, seed = 0)
# test1(params,A,shrinkage,d,alg,dim4free,threads,saturation_ratio)


#LWE instance 
(n,alpha) = (65,0.015)
succ_prob = 0.1
A,q = load_lwe_lattice(n,alpha)
d = A.ncols
g6k = Siever(A, params)
slope = basis_quality(g6k.M)["/"]
# rr = [g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)]
# blocksizes = [86, 88, 107, 107, 107]
# blocksizes = [86]
blocksizes = []
print("Estimation succ_prob = %.2f, saturation_ratio = %.2f, pre-process blocksize strategy: " %(succ_prob,saturation_ratio),end='')
print(blocksizes)
print("Intial Slope = %.5f\n" % slope)


pump_params = {"down_sieve": True,"saturation_error":"ignore"}
T0 = time.time()
for blocksize in blocksizes:
    print("Starting a pnjBKZ-%d tour. " % (blocksize))
    T0_BKZ = time.time()
    T_pumps = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=jump,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun=default_dim4free_fun,
                                            pump_params=pump_params)

    slope = basis_quality(g6k.M)["/"]
    T_BKZ = time.time()-T0_BKZ
    fmt = "slope: %.5f, walltime: %.3f sec"
    print(fmt % (slope, T_BKZ))     
rr = [g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)]   

kappa,beta,dim4free = pump_estimation1(rr,q,alpha)
print("estimation1: (%d,%d,%d), sieve dim = %d" %(kappa,dim4free,d,d-kappa-dim4free))  
kappa,beta,dim4free = pump_estimation2(rr,q, alpha,succ_prob)
dim4free = max(d-kappa-146, dim4free)
print("estimation2: (%d,%d,%d), sieve dim = %d" %(kappa,dim4free,d,d-kappa-dim4free))  

# print(kappa,beta,f)
# target_norm =  1.5 * (((alpha*q)**2) * (d-1) + 1 )
# gh = gaussian_heuristic(rr)
# shrinkage = target_norm / gh
shrinkage = 1
for _ in range(5):
    shrinkage = control_group(params,A,d,dim4free+kappa,alg,threads,saturation_ratio)
    # shrinkage = control_group(params,A,d,d//2,alg,threads,saturation_ratio)
    if shrinkage < 1:
        break
    
# test1(params,A,shrinkage,d,alg,dim4free+kappa,threads,saturation_ratio)
print("walltime: %.3f s" %(time.time()-T0))