# -*- coding: utf-8 -*-
"""
BKZ Tours.
"""
import sys
from .pump import pump
#from .pump_cpu import pump
from .workout import workout
import six
import psutil  
import os
from math import log,pi
import time
from fpylll.util import gaussian_heuristic
# from ...pro_pnjBKZ_simulator.codes.util

try:
    basestring
except NameError:
    basestring = str


def get_current_slope(r, start_row=0, stop_row=-1):
    """
    A Python re-implementation of ``MatGSO.get_current_slope``.

        >>> from fpylll import IntegerMatrix, GSO, LLL, FPLLL
        >>> FPLLL.set_random_seed(1337)
        >>> A = IntegerMatrix.random(100, "qary", bits=30, k=50)
        >>> _ = LLL.reduction(A)
        >>> M = GSO.Mat(A); _ = M.update_gso()
        >>> from fpylll.tools.quality import get_current_slope
        >>> M.get_current_slope(0, 100)  # doctest: +ELLIPSIS
        -0.085500625...
        >>> get_current_slope(M.r(), 0, 100) # doctest: +ELLIPSIS
        -0.085500625...

    """
    x = [log(r[i]) for i in range(start_row, stop_row)]
    n = stop_row - start_row
    i_mean = (n - 1) * 0.5 + start_row
    x_mean = sum(x)/n
    v1, v2 = 0.0, 0.0
    for i in range(stop_row - start_row):
        v1 += (i - i_mean) * (x[i] - x_mean)
        v2 += (i - i_mean) * (i - i_mean)
    return v1 / v2


def dim4free_wrapper(dim4free_fun, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize)
    return int(min((blocksize - 40)/2, dim4free))


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)

# def theo_dim4free_fun(blocksize):
#     """
#     Theoretical Dimension-for-free function in [Duc18]
#     """

#     return int(blocksize*log(4/3.)/log(blocksize/2./pi)) 



def naive_bkz_tour(g6k, tracer, blocksize, dim4free_fun=default_dim4free_fun,
                   extra_dim4free=0, workout_params=None, pump_params=None):
    """
    Run a naive BKZ-tour: call ``workout`` as an SVP oracle consecutively on
    each block.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param blocksize: dimension of the blocks
    :param dim4free_fun: number of dimension for free as a function of beta (function, or string e.g. `lambda x: 11.5+0.075*x`)
    :param extra_dim4free: increase the number of dims 4 free (blocksize is increased, but not sieve dimension)
    :param workout_params: parameters to pass to the workout
    :param pump_params: parameters to pass to the pump

    """
    if workout_params is None:
        workout_params = {}

    if "dim4free_min" in workout_params:
        raise ValueError("In naive_bkz, you should choose dim4free via dim4free_fun.")

    d = g6k.full_n

    if isinstance(dim4free_fun, basestring):
        dim4free_fun = eval(dim4free_fun)

    dim4free = dim4free_wrapper(dim4free_fun, blocksize) + extra_dim4free
    blocksize += extra_dim4free

    for kappa in range(d-3):
        beta = min(blocksize, d - kappa)
        lost_dim = blocksize - beta
        f = max(dim4free - lost_dim, 0)

        workout(g6k, tracer, kappa, beta, f, pump_params=pump_params, **workout_params)
        g6k.lll(0, d)


def pump_n_jump_bkz_tour(g6k, tracer, blocksize, jump=1,
                         dim4free_fun=default_dim4free_fun, extra_dim4free=0,
                         pump_params=None, goal_r0=0., verbose=False):
    """
    Run a PumpNjump BKZ-tour: call Pump consecutively on every (jth) block.

    :param g6k: The g6k object to work with
    :param tracer: A tracer for g6k
    :param blocksize: dimension of the blocks
    :param jump: only call the pump every j blocks
    :param dim4free_fun: number of dimension for free as a function of beta (function, or string
        e.g. `lambda x: 11.5+0.075*x`)
    :param extra_dim4free: increase the number of dims 4 free (blocksize is increased, but not sieve
        dimension)
    :param pump_params: parameters to pass to the pump
    """
    if pump_params is None:
        pump_params = {"down_sieve": False}

    if "dim4free" in pump_params:
        raise ValueError("In pump_n_jump_bkz, you should choose dim4free via dim4free_fun.")

    d = g6k.full_n
    g6k.shrink_db(0)
    g6k.lll(0,d)
    g6k.update_gso(0,d)

    if isinstance(dim4free_fun, six.string_types):
        dim4free_fun = eval(dim4free_fun)
    

    # file_name = "80-005-gpu-32-thread-gs-lengths-%d-%d.txt" %(g6k.M.B.ncols,blocksize)


    dim4free = dim4free_wrapper(dim4free_fun, blocksize) + extra_dim4free
    blocksize += extra_dim4free

    indices  = [(0, blocksize - dim4free + i, i) for i in range(0, dim4free, jump)]
    indices += [(i, blocksize, dim4free) for i in range(0, d - blocksize, jump)]
    indices += [(d - blocksize + i, blocksize - i, dim4free - i) for i in range(0, dim4free, jump)]

    pump_params["down_stop"] = dim4free+3


    max_RAM_cost = 0
  
    for (kappa, beta, f) in indices:
        if verbose:
            print("\r k:%d, b:%d, f:%d , RAM cost: %.4f GB" % (kappa, beta, f, psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
            sys.stdout.flush()
            # print()
            
            RAM_cost = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
            if max_RAM_cost < RAM_cost:
                max_RAM_cost = RAM_cost
        
        # T0 = time.time()
        pump(g6k, tracer, kappa, beta, f, **pump_params)
        # rr = [g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)]
        # File.write(str(rr))
        # File.write('\n')
        # T_pump = time.time()-T0
        # T_pumps.append(T_pump)
        
        # slopes.append(get_current_slope(rr,kappa+f,kappa+beta))
        # ghs.append(gaussian_heuristic(rr[kappa+f:kappa+beta]))
        # print("k:%d, b:%d, f:%d  , RAM cost: %.4f GB, Pump Cost: %.4f s " % (kappa, beta, f,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024,T_pump))

        # raise RuntimeError("Debug")
        g6k.lll(0, d)
        if g6k.M.get_r(0, 0) <= goal_r0:
            return max_RAM_cost

    if verbose:
        
        print("\r k:%d, b:%d, f:%d  , RAM cost: %.4f GB " % (d-(blocksize-dim4free), blocksize-dim4free, 0,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024), end=' ')
        # print("k:%d, b:%d, f:%d  , RAM cost: %.4f GB, Pump Cost: %.4f s " % (d-(blocksize-dim4free), blocksize-dim4free, 0,psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024,T_pump))
        sys.stdout.flush()

        RAM_cost = psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024
        if max_RAM_cost < RAM_cost:
            max_RAM_cost = RAM_cost

    
    if verbose:
        print("\r k:%d, b:%d, f:%d " % (d-(blocksize-dim4free), blocksize-dim4free, 0))
        
        sys.stdout.flush()

    pump_params["down_stop"] = blocksize - dim4free
    
    T_0 = time.time()
    pump(g6k, tracer, d-(blocksize-dim4free), blocksize-dim4free, 0, **pump_params)
    # rr = [g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)]
    # File.write(str(rr))
    # File.write('\n')
    # T_pump = time.time()-T0
    # T_pumps.append(T_pump)
    # slopes.append(get_current_slope(rr,d-(blocksize-dim4free),d))
    # # ghs.append(gaussian_heuristic(rr[d-(blocksize-dim4free):d]))
    if verbose:
        print('')
        sys.stdout.flush()

    # print(T_pumps)
    # File.close()
    return max_RAM_cost
    # return T_pumps,slopes,ghs #max_RAM_cost,
