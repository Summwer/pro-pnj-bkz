
from g6k.utils.util import load_svpchallenge_and_randomize
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from g6k.siever import Siever, SieverParam
from fpylll.tools.bkz_stats import dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from g6k.algorithms.pump import pump
from fpylll.util import gaussian_heuristic
import time



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
    params = {"gpus": gpus, "threads": threads, "saturation_ratio": 0.375, "db_size_factor": 2.77, "max_nr_buckets": 0 }
    params = SieverParam(params)
    g6k = Siever(A, params)
    T0 = time.time()
    if blocksize < 50:
        bkz = BKZReduction(g6k.M)
        par = fplll_bkz.Param(blocksize, strategies=fplll_bkz.DEFAULT_STRATEGY,
                                  max_loops=1)
        bkz(par)
    else:
        pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=jump,
                                     extra_dim4free=extra_dim4free,
                                     dim4free_fun=dim4free_fun,
                                     pump_params=pump_params)
    return time.time()-T0

def one_pump_cost(sieve_dim, n=180,gpus = 2, threads = 32, llb = 0):
    A, bkz = load_svpchallenge_and_randomize(n, s=0, seed=0)
    params = {"gpus": gpus, "threads": threads, "saturation_ratio": 0.375, "db_size_factor": 2.77, "max_nr_buckets": 0 }
    params = SieverParam(params)
    g6k = Siever(A, params)
    
    d = g6k.full_n
    f  = d - sieve_dim
    T0 = time.time()
    pump(g6k, dummy_tracer, llb, d-llb, f)
    
    return time.time()-T0

def PnjBKZCostTest(n=180,tours = 10, gpus = 2, threads = 32):
    print("Start PnjBKZ cost test...")
    AvgTs = {}
    A, bkz = load_svpchallenge_and_randomize(s=0, seed=0)
    params = {"gpus": gpus, "threads": threads, "saturation_ratio": 0.375, "db_size_factor": 2.77, "max_nr_buckets": 0 }
    params = SieverParam(params)
    g6k = Siever(A, params)
    rr = [g6k.M.get_r(i, i) for i in range(g6k.full_n)]
    for jump in range(1,theo_dim4free2_in_B(rr)):
        for blocksize in range(10,136,2):
            Ts = []
            for _ in tours:
                Ts.append(one_tour_pnjbkz_cost_test(blocksize, jump, n = n, gpus = gpus, threads = threads))
            AvgTs[(blocksize,jump)] = Ts/tours
        print("(%d,%d): %.2f sec" %(blocksize,jump,AvgTs[(blocksize,jump)]))
    print(AvgTs)
    return AvgTs



def PumpCostTest(n=180,tours = 10, gpus = 2, threads = 32):
    print("Start Pump cost test...")
    AvgTs = {}
    for sieve_dim in range(50,136,2):
        Ts = []
        for _ in tours:
            Ts.append(one_pump_cost(sieve_dim, n = n, gpus = gpus, threads = threads))
        AvgTs[sieve_dim] = Ts/tours
        print("%d: %.2f sec" %(sieve_dim,AvgTs[sieve_dim]))
    print(AvgTs)
    return AvgTs
    
if __name__ =="__main__":
    # PnjBKZCostTest( gpus = 2, threads = 32)
    PumpCostTest( gpus = 2, threads = 32)