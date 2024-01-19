
from g6k.utils.util import load_svpchallenge_and_randomize
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from g6k.siever import Siever
from g6k.siever_params import SieverParams
from fpylll.tools.bkz_stats import dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour, default_dim4free_fun, dim4free_wrapper
from g6k.algorithms.pump import pump
from fpylll.util import gaussian_heuristic
import time
import psutil
import os


#d4f(B): postive 
def theo_dim4free2_in_B(rr):
    gh = gaussian_heuristic(rr)
    d = len(rr)
    for f in range(d-1,-1,-1):
        ghf = gaussian_heuristic(rr[f:])
        if(ghf * 4/3. >=  ((d-f)/d) * gh):
            return f
    return 0

def one_tour_pnjbkz_cost_test(blocksize, jump, n = 180, gpus = 2, threads = 32, extra_dim4free = 12, dim4free_fun = "default_dim4free_fun", pump_params = {"down_sieve": True }):
    A, bkz = load_svpchallenge_and_randomize(n, s=0, seed=0)

    params = SieverParams(gpus = gpus, threads = threads, saturation_ratio = 0.375, db_size_factor = 2.77, max_nr_buckets = 0)
    g6k = Siever(A, params)
    T0 = time.time()
    if blocksize < 50:
        bkz = BKZReduction(g6k.M)
        par = fplll_bkz.Param(blocksize, strategies=fplll_bkz.DEFAULT_STRATEGY,
                                  max_loops=1)
        bkz(par)
        RAM_cost = 0.
    else:
        RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=jump,
                                     extra_dim4free=extra_dim4free,
                                     dim4free_fun=dim4free_fun,
                                     pump_params=pump_params)
    return time.time()-T0, RAM_cost

def one_pump_cost(sieve_dim, n=180,gpus = 2, threads = 32, llb = 0):
    A, bkz = load_svpchallenge_and_randomize(n, s=0, seed=0)
    params = SieverParams(gpus = gpus, threads = threads, saturation_ratio = 0.375, db_size_factor = 2.77, max_nr_buckets = 0)
    g6k = Siever(A, params)
    
    d = g6k.full_n
    f  = d - sieve_dim
    T0 = time.time()
    _ , RAM_cost = pump(g6k, dummy_tracer, llb, d-llb, f, down_sieve=True)
    
    return time.time()-T0, RAM_cost

def PnjBKZCostTest(n=180,tours = 1, gpus = 2, threads = 32):
    print("Start PnjBKZ cost test...")
    AvgTs = {}
    AvgRAMs = {}
    A, bkz = load_svpchallenge_and_randomize(s=0, seed=0)
    params = SieverParams(gpus = gpus, threads = threads, saturation_ratio = 0.375, db_size_factor = 2.77, max_nr_buckets = 0)
    g6k = Siever(A, params)
    rr = [g6k.M.get_r(i, i) for i in range(g6k.full_n)]
    print("{0: <10} {1: <10} {2: <10} {3: <10} {4: <10}".format("beta", "jump","sieve_dim","AvgT/s","AvgRAM/GB"))
    for jump in range(1,theo_dim4free2_in_B(rr)):
        for blocksize in range(51,136,2):
            Ts = []
            RAMs = []
            for _ in range(tours):
                TR = one_tour_pnjbkz_cost_test(blocksize, jump, n = n, gpus = gpus, threads = threads)
                Ts.append(TR[0])
                RAMs.append(TR[1])
            AvgTs[(blocksize,jump)] = round(sum(Ts)/tours,2)
            AvgRAMs[(blocksize,jump)] = round(sum(RAMs)/tours,2)
            sieve_dim = blocksize - dim4free_wrapper(default_dim4free_fun, blocksize)
            print("{0: <10} {1: <10} {2: <10} {3: <10} {4: <10}".format(blocksize,jump, sieve_dim,AvgTs[(blocksize,jump)], AvgRAMs[(blocksize,jump)]))
            # print("Cost of PnjBKZ-(%d,%d), sieve_dim = %d: %.2f sec, RAM cost:%.4f GB" %(blocksize,jump, sieve_dim,AvgTs[(blocksize,jump)], AvgRAMs[(blocksize,jump)]))
    print("AvgTs:",end="")
    print(AvgTs)
    print("AvgRAMs:",end="")
    print(AvgRAMs)
    return AvgTs, AvgRAMs



def PumpCostTest(n=180,tours = 1, gpus = 2, threads = 32):
    print("Start Pump cost test...")
    AvgTs = {}
    AvgRAMs = {}
    
    print("{0: <10} {1: <10} {2: <10} {3: <10} {4: <10}".format("kappa", "dim","sieve_dim","AvgT/s","AvgRAM/GB"))
    # for sieve_dim in range(51,136,2):
    for sieve_dim in range(129,136,2):
        Ts = []
        RAMs = []
        for _ in range(tours):
            TR = one_pump_cost(sieve_dim, n = n, gpus = gpus, threads = threads)
            Ts.append(TR[0])
            RAMs.append(TR[1])
        AvgTs[sieve_dim] = round(sum(Ts)/tours,2)
        AvgRAMs[sieve_dim] = round(sum(RAMs)/tours,2)
        print("{0: <10} {1: <10} {2: <10} {3: <10} {4: <10}".format(0, n,sieve_dim,AvgTs[sieve_dim], AvgRAMs[sieve_dim]))
        # print("AvgCost of Pump-{0,%d,%d}: %.2f sec, RAM cost:%.4f GB" %(n,n-sieve_dim,AvgTs[sieve_dim], AvgRAMs[sieve_dim]))
    print("AvgTs:",end="")
    print(AvgTs)
    print("AvgRAMs:",end="")
    print(AvgRAMs)
    return AvgTs, AvgRAMs

    
if __name__ =="__main__":
    PumpCostTest( gpus = 2, threads = 32)
    PnjBKZCostTest( gpus = 2, threads = 32)