#!/usr/bin/env python
# -*- coding: utf-8 -*-
####
#
#   Copyright (C) 2018-2021 Team G6K
#
#   This file is part of G6K. G6K is free software:
#   you can redistribute it and/or modify it under the terms of the
#   GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version.
#
#   G6K is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with G6K. If not, see <http://www.gnu.org/licenses/>.
#
####


"""
LWE Challenge Solving Command Line Client
"""

# from __future__ import absolute_import
# from __future__ import print_function
import copy
import re
import sys
import time

from collections import OrderedDict # noqa
from math import log

from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality
from fpylll.util import gaussian_heuristic as gh
from numpy import True_
from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL

from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from g6k.algorithms.pump import pump
from g6k.siever import Siever
from g6k.utils.util import load_svpchallenge_and_randomize,load_lwe_challenge
from g6k.siever_params import SieverParams
from g6k.utils.lwe_estimation import gsa_params, primal_lattice_basis
from g6k.utils.cost import dim4free_wrapper, default_dim4free_fun, practical_pump_cost,dims4free, get_k1_k2_pump,theo_dim4free_fun1
from g6k.utils.stats import dummy_tracer

from os import mkdir
import os

def gen_original_gs(n,alpha, goal_margin):
    """
    Run the primal attack against Darmstadt LWE instance (n, alpha).

    :param n: the dimension of the LWE-challenge secret
    :param alpha: the noise rate of the LWE-challenge
    """

    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    
    #print("-------------------------")
    #print("Primal attack, LWE challenge n=%d, alpha=%.4f" % (n, alpha))
    m = None
    if m is None:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True)
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
    #print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    #print()

    # no use in having a very small b
    b = max(b, s-65)
    
    
    B = primal_lattice_basis(A, c, q, m=m) #debug

    params = None
    g6k = Siever(B, params)
    #print("GSO precision: ", g6k.M.float_type)

    d = g6k.full_n

    g6k.lll(0, g6k.full_n)
    #slope = basis_quality(g6k.M)["/"]
    # print("Intial Slope = %.5f\n" % slope)

    log_rr0 = [log(g6k.M.get_r(i, i)) for i in range(d)]


    target_norm = goal_margin * (alpha*q)**2 * m + 1

    return b,s,target_norm,log_rr0,0,q
    
def gen_svp_instance(d):
    params = SieverParams()
    A, bkz = load_svpchallenge_and_randomize(d)
    g6k = Siever(A, params)

    return [g6k.M.get_r(i, i) for i in range(d)]

# threads = 32, gpus = 2, Known maxT, find the dimension(beta-f).
def get_pump_dim(d,T_max):
    for beta in range(min(d,145)):
        f = dim4free_wrapper(dims4free,beta)
        T_pump = 2**(practical_pump_cost(beta)[0])
        if T_pump > T_max:
            beta = beta -1
            f = dim4free_wrapper(dims4free,beta)
            llb = d - beta
            if beta > 0:
                return llb,beta,f
            
    for f in range(dim4free_wrapper(dims4free,min(d,145)),0,-1):
        beta_prime = beta - f
        k1, k2 = get_k1_k2_pump(beta_prime) # threads = 20
        # k = (1/71.)*((1.33)**(beta/10.))
        T_pump = 2**(k1*beta_prime+k2)
        if T_pump > T_max:
            f = f+1
            llb = d - beta
            if beta > 0:
                return llb,beta,f
    return 0,min(d,144),0


def load_lwe_challenge_mid(n=40, alpha=0.005):
    """
    Load LWE challenge from file or website.

    :param n: LWE dimension
    :param alpha: the *standard deviation* of the secret is alpha*q

    """
    alpha = int(round(alpha * 1450))
    start = "lwechal_midmat"

    if not os.path.isdir(start):
        os.mkdir(start)
    
    end = "{n:03d}-{alpha:03d}-midmat.txt".format(n=n, alpha=alpha)
    filename = os.path.join(start, end)
    try:
        data = open(filename, "r").readlines()
    except FileNotFoundError:
        return None
    n, m, q = [int(x) for x in [data[0], data[1], data[2]]]
    c_index = 3 if data[3].startswith("[") else 4
    #A = eval(",".join([s_.replace(" ", ", ") for s_ in data]))
    B = eval(",".join([s_.replace(" ", ", ") for s_ in data[c_index:]]))
    B = IntegerMatrix.from_matrix(B)
 
    
    return B,q


def load_svp_midmat(d):
    """
    Load svp challenge midmat from file.

    :param d: svp dimension.
    
    """

    start = "svpchallenge"

    if not os.path.isdir(start):
        os.mkdir(start)
    
    end = "{d:03d}-midmat.txt".format(d=d)
    filename = os.path.join(start, end)
    try:
        data = open(filename, "r").readlines()
    except FileNotFoundError:
        return False
    B = eval(",".join([s_.replace(" ", ", ") for s_ in data[0:]]))
    B = IntegerMatrix.from_matrix(B)
 
    
    return B


#for svp instance
def store_svp_midmat(d,g6k):
    filename = 'svpchallenge/%03d-midmat.txt' % (d)
    fn = open(filename, "w")
    fn.write('[')
    for i in range(g6k.M.B.nrows):
        fn.write('[')
        for j in range(g6k.M.B.ncols):
            fn.write(str(g6k.M.B[i][j]))
            if j<g6k.M.B.ncols-1:
                fn.write(' ')
        if i < g6k.M.B.nrows-1:
            fn.write(']\n')
    fn.write(']]')
    fn.close()
    
    
def BKZ_Pump_comparison(n,alpha, d, prebeta, bkz_beta = None , pump_dsvp = None , pump_f = None):
    
    params = SieverParams(threads = 32, gpus = 2)
    pump_params = {"down_sieve": True}
    
    B = load_svp_midmat(d)
    if(not B):
        # A, bkz = load_svpchallenge_and_randomize(d)
        A, c, q = load_lwe_challenge(n=n, alpha=alpha)
        # min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                            # samples=A.nrows, decouple=True)
        # (b, s, m) = min_cost_param
        m = d - 1
        B = primal_lattice_basis(A, c, q, m=m)
        g6k = Siever(B, params)
        g6k.lll(0, g6k.full_n)
        for blocksize in range(10,prebeta):
            bkz = BKZReduction(g6k.M)
            par = fplll_bkz.Param(blocksize,strategies=fplll_bkz.DEFAULT_STRATEGY,max_loops=1)
            bkz(par)
        store_svp_midmat(d,g6k)
    else:
        g6k = Siever(B, params)
    
    slope = basis_quality(g6k.M)["/"]
    print( "-------------------------")
    print( "n=%d, alpha=%.4f, dim = %d" % (n, alpha, d))
    # print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    print()
    print("Intial Slope = %.5f\n" % slope)
    print("ln(||b0||) = %f. \n" %log(g6k.M.get_r(0, 0)))
    
    # #--------------------------bkz test---------------------------# 
    A = load_svp_midmat(d)
    if(bkz_beta == None):
        beta =  d // 2 #+ 10
    else:
        beta = bkz_beta
    T0 = time.time()
    print("Starting a pnjBKZ-%d tour." % (beta))
    pump_n_jump_bkz_tour(g6k, dummy_tracer, beta, verbose=True, extra_dim4free=12,
                                     dim4free_fun=default_dim4free_fun, pump_params=pump_params)
    T_BKZ = time.time() - T0
    print("walltime: %f sec." % T_BKZ)
    g6k.lll(0, g6k.full_n)
    print("ln(||b0||) = %f. \n" %log(g6k.M.get_r(0, 0)))
    
    
    #--------------------------pump test---------------------------#
    
    B = load_svp_midmat(d)
    g6k = Siever(B, params)
    slope = basis_quality(g6k.M)["/"]
    
    print("Loaded challenge dim %d\n" % d)
    print("Intial Slope = %.5f\n\n" % slope)
    d = g6k.full_n
    
    # if get_pump_dim(d,T_BKZ) != None:
    # n_max = min(d,int(58 + 2.85 * log(T_BKZ * params.threads)/log(2.)))
    # f = dim4free_wrapper(default_dim4free_fun,n_max)
    # llb = d - n_max
    if(pump_dsvp == None or pump_f == None):
        print(pump_dsvp,d)
        llb,n_max,f = get_pump_dim(d,T_BKZ)
    else:
        llb,n_max,f = max(0,d - pump_dsvp), min(pump_dsvp,d), pump_f
    print("Starting svp pump_{%d, %d, %d}." % (llb, n_max, f)) # noqa
    T0 = time.time()
    pump(g6k, dummy_tracer, llb, n_max, f, verbose=True)
    T_pump = time.time() - T0
    print("walltime: %f sec." % T_pump)
    g6k.lll(0, g6k.full_n)
    print("ln(||b0||) = %f. \n" %log(g6k.M.get_r(0, 0)))
    # file.close()



# dim = 60
# bkz_beta = 28
# pump_dsvp = dim
# f = dim4free_wrapper(dims4free, pump_dsvp)
# # f = dim4free_wrapper(theo_dim4free_fun1, pump_dsvp)
# # f = dim4free_wrapper(default_dim4free_fun, pump_dsvp)
# BKZ_Pump_comparison(dim, bkz_beta , pump_dsvp, f)



# dim = 80
# bkz_beta = 30
# pump_dsvp = dim
# f = dim4free_wrapper(dims4free, pump_dsvp)
# # f = dim4free_wrapper(theo_dim4free_fun1, pump_dsvp)
# # f = dim4free_wrapper(default_dim4free_fun, pump_dsvp)
# BKZ_Pump_comparison(dim, bkz_beta , pump_dsvp, f)



# (n, alpha, d) = (50, 0.020, 206)
# (n, alpha, d) = (40, 0.025, 172)
# (n, alpha, d) = (80, 0.005, 271)

# n_alpha_d = [(80,0.005,271)]
# (n, alpha, d) = (40,0.040,191)
# prebeta = 30
# bkz_beta = 80
# pump_dsvp = None
# f = None
# # f = dim4free_wrapper(dims4free, pump_dsvp)
# # f = dim4free_wrapper(theo_dim4free_fun1, pump_dsvp)
# # f = dim4free_wrapper(default_dim4free_fun, pump_dsvp)
# BKZ_Pump_comparison(n,alpha, d, prebeta, bkz_beta , pump_dsvp, f)

(n, alpha, d) = (40, 0.025, 172)
(n, alpha, d) = (40, 0.005, 90)

prebeta = 50
bkz_beta = 80
pump_dsvp = None
f = None
# f = dim4free_wrapper(dims4free, pump_dsvp)
# f = dim4free_wrapper(theo_dim4free_fun1, pump_dsvp)
# f = dim4free_wrapper(default_dim4free_fun, pump_dsvp)
#(80, 0.005, 271), 
# for (n,alpha,d) in [(40, 0.025, 172),(60,0.010,222)]:
# for (n,alpha,d) in [(55,0.020,231), (45,0.030,202), (40,0.040,191)]:
for (n,alpha,d) in [(40, 0.020, 163), (45, 0.010, 166), (50,0.010,184), (55,0.005,184), (60,0.005,199), (65,0.005,219)]:
    BKZ_Pump_comparison(n,alpha, d, prebeta, bkz_beta , pump_dsvp, f)