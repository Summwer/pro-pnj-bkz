#!/usr/bin/env python
# -*- coding: utf-8 -*-


from locale import MON_10
from math import e, lgamma, log, pi

from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL
from fpylll.tools.bkz_simulator import simulate
from fpylll.util import gaussian_heuristic

from g6k.algorithms.bkz import default_dim4free_fun
from g6k.utils.util import load_lwe_challenge
from pro_pnjBKZ_simulator.codes.usvp_simulator_v2 import gen_bkzstrategy_min_cost_v2
# from pro_pnjBKZ_simulator.codes.usvp_simulator_v3 import gen_bkzstrategy_min_cost_v3
from numpy import floor

from lwe_estimation_old import gsa_params
import time
from strategy_test import read_preprocess_data
import os


def delta_0f(k):
    """
    Auxiliary function giving root Hermite factors. Small values
    experimentally determined, otherwise from [Chen13]

    :param k: BKZ blocksize for which the root Hermite factor is required

    """
    small = (( 2, 1.02190),  # noqa
             ( 5, 1.01862),  # noqa
             (10, 1.01616),
             (15, 1.01485),
             (20, 1.01420),
             (25, 1.01342),
             (28, 1.01331),
             (40, 1.01295))

    k = float(k)
    if k <= 2:
        return (1.0219)
    elif k < 40:
        for i in range(1, len(small)):
            if small[i][0] > k:
                return (small[i-1][1])
    elif k == 40:
        return (small[-1][1])
    else:
        return (k/(2*pi*e) * (pi*k)**(1./k))**(1/(2*(k-1.)))


def log_gh_svp(d, delta_bkz, svp_dim, n, q):
    """
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

    """
    d = float(d)
    svp_dim = float(svp_dim)
    ball_part = ((1./svp_dim)*lgamma((svp_dim/2.)+1))-(.5*log(pi))
    vol_part = ((1./d)*(d-n-1)*log(q))+((svp_dim-d)*log(delta_bkz))
    return ball_part + vol_part





def primal_lattice_basis(A, c, q, m=None):
    """
    Construct primal lattice basis for LWE challenge
    ``(A,c)`` defined modulo ``q``.

    :param A: LWE matrix
    :param c: LWE vector
    :param q: integer modulus
    :param m: number of samples to use (``None`` means all)

    """
    if m is None:
        m = A.nrows
    elif m > A.nrows:
        raise ValueError("Only m=%d samples available." % A.nrows)
    n = A.ncols

    B = IntegerMatrix(m+n+1, m+1)
    for i in range(m):
        for j in range(n):
            B[j, i] = A[i, j]
        B[i+n, i] = q
        B[-1, i] = c[i]
    B[-1, -1] = 1

    B = LLL.reduction(B)
    assert(B[:n] == IntegerMatrix(n, m+1))
    B = B[n:]

    return B


def determineM(beta_bound, svp_bound, m, stddev, q):
    for bkz_block_size in range(40, beta_bound):
            delta_0 = delta_0f(bkz_block_size)
            svp_dims = range(40, svp_bound)
            for svp_dim in svp_dims:
                d = float(m + 1)
                rhs = log_gh_svp(d, delta_0, svp_dim, n, q)
                if rhs - log(stddev) - log(svp_dim)/2. >= 0:
                    return True
    return False

def SamplesRangeDeterminationAlgorithm(n,alpha):
    A, _, q = load_lwe_challenge(n, alpha)
    samples = A.nrows
    stddev = alpha*q

    mrange = []
    ms = range(n, min(5*n+1, samples+1))

    for m in ms:
        beta_bound = min(m+1, 200+default_dim4free_fun(200)+1)
        svp_bound = min(m+1, 200)
        if determineM(beta_bound, svp_bound, m, stddev, q):
            mrange.append(m)
    return mrange



def our_lwe_estimation(n,alpha,tau=40):
    mrange = SamplesRangeDeterminationAlgorithm(n,alpha)
    #print(mrange)
    (beta, svp_dim, d) = gsa_params(n, alpha,decouple=True) #lwe estimation function in default G6K
    print((beta, svp_dim, d))
    m0 = d - 1
    print(n,alpha,m0)
    m_max, m_min = max(mrange), min(mrange)
    mincost = None
    minblocksize_strategy = []
    minblocksize_strategy, mincost, m0 = our_lwe_estimation_sub_function(m0,minblocksize_strategy, mincost,n,alpha,tau,m_max,m_min)
    
    return minblocksize_strategy, mincost, m0


def our_lwe_estimation_sub_function(m0,minblocksize_strategy, mincost,n,alpha,tau,m_max,m_min):
    succ_prob = 0.80
    A, c, q = load_lwe_challenge(n, alpha)
    # data_dir = "pro_pnjBKZ_simulator/strategy_pre_process"
    # file_name = data_dir + "/%d-%s.log" %(n,str(alpha)[2:])
    # data = read_preprocess_data(file_name)
    # (ebeta,svp_dim,target_norm,log_rr0,preprocess_cost) = data
    # print(log_rr0)
    # print(d)
    m_all = A.nrows #m_all
    stddev = alpha*q

    # m0
    B = primal_lattice_basis(A, c, q, m=m0)
    B = LLL.reduction(B)
    M = GSO.Mat(B)
    M.update_gso()
    log_rr0 = [log(M.get_r(i, i)) for i in range(M.B.nrows)]
    
    file_path = "pro_pnjBKZ_simulator/estimate_m_v2_%.2f" %succ_prob
    try:
        os.mkdir(file_path)
    except FileExistsError:
        pass

    file_path = "pro_pnjBKZ_simulator/estimate_m_v2_%.2f/%d" %(succ_prob,m0)

    try:
        os.mkdir(file_path)
    except FileExistsError:
        pass

    d = len(log_rr0)
    goal_margin = 1.5
    MAX_TIME = float("inf")
    
    jump = 1
    blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen = gen_bkzstrategy_min_cost_v2(file_path,n,str(alpha),log_rr0, d, jump, goal_margin,MAX_TIME,q,0) 
    if mincost == None or mincost > total_cost:
        minblocksize_strategy = blocksize_strategy
        mincost = total_cost
        
    print(minblocksize_strategy, mincost,bkz_cost, pump_cost,m0)
    if tau == 0:
        return minblocksize_strategy, mincost, m0
    
    m1 = m0
    jump = 1

    for m in (max(n,m0-tau), min(m_all,m0+tau)):
        if m <= m_max and m >= m_min:
            d = m+1
            target_norm = stddev**2 * m + 1
            B = primal_lattice_basis(A, c, q, m = m)
            B = LLL.reduction(B)
            M = GSO.Mat(B)
            M.update_gso()
            log_rr0 = [log(M.get_r(i, i)) for i in range(M.B.nrows)]

            blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen = gen_bkzstrategy_min_cost_v2(file_path,n,str(alpha),log_rr0, d, jump, goal_margin,MAX_TIME,q,0)
            if mincost == None or mincost > total_cost:
                minblocksize_strategy = blocksize_strategy
                mincost = total_cost
                m1 = m
    if m1 == m0:
        tau = int(floor(tau/2.0))

    m0 = m1
    return our_lwe_estimation_sub_function(m0,minblocksize_strategy, mincost,n,alpha,tau,m_max,m_min)

# n = 90
# alpha = 0.005
lwes = [(80, '0.005'),(90,'0.005'),(50,'0.025'),(55,'0.020'),(45,'0.030'),(60,'0.015'),(40,'0.035')]
for (n,alpha) in lwes:
    alpha = float(alpha)
    T0 = time.time()
    print("Start..........................")
    minblocksize_strategy, mincost, m0 = our_lwe_estimation(n,alpha,40)

    print(minblocksize_strategy, mincost, m0)

    print("Time Cost for finding optimized m0: %f s" %(time.time()-T0))