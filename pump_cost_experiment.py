from g6k.utils.util import load_svpchallenge_and_randomize
from g6k.utils.util import load_lwe_challenge, load_matrix_file
# from g6k.utils.lwe_estimation_old import gsa_params, primal_lattice_basis
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from g6k.siever import Siever
from g6k.siever_params import SieverParams
from fpylll.tools.quality import basis_quality
import time
from g6k.algorithms.pump import pump
from pro_pnjBKZ_simulator.codes.util import  dim4free_wrapper, default_dim4free_fun,theo_dim4free_fun1
from fpylll.util import gaussian_heuristic

def f_0(blocksize):
    return 0

#test for saturation change.
def test1(d, blocksize, saturation_ratio, threads):
    d = d
    blocksizes = [blocksize]
    tours = 1
    saturation_ratio = saturation_ratio
    threads = threads
    gpus = 2
    dim4free = dim4free_wrapper(default_dim4free_fun,blocksize)
    extra_dim4free = 12
    pump_params = {"down_sieve": True,"saturation_error":"ignore"}#,"start_up_n": blocksize - dim4free}
    print("Start pump cost test with saturation_ratio = %.2f:" %saturation_ratio)
    
    file_name = "[20221011]pump_cost_in_pnjbkz_and_progressive_sieve/test1/svp-%d-%d-xf=%d-f=%d-%d.txt" % (d,blocksize,extra_dim4free,dim4free,int(saturation_ratio * 100))
    
    f = open(file_name,"w")
    bkzCost = {}
    RAMCost = {}
    for blocksize in blocksizes:
        T_tmp = []
        RAM_tmp = []
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            # f.write("Starting a pnjBKZ-%d tour. " % (blocksize))
            params = SieverParams(threads=threads,gpus=gpus,saturation_ratio=saturation_ratio)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f, beta-f = %d" % (slope,blocksize-dim4free))
            # f.write("Intial Slope = %.5f \n" % slope)
            T_pumps = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=extra_dim4free,
                                            dim4free_fun=default_dim4free_fun,
                                            pump_params=pump_params)
            
            # f.write("Pump Cost in pnjBKZ: ")
            f.write(str(T_pumps))   
            T_BKZ = time.time()-T0_BKZ
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            # f.write(fmt % (slope, T_BKZ) +"\n")
    
        f.flush()

    f.close()


#Test the relationship between pump time cost and the slope of reduced basis.
def test2(d, blocksize, threads):
    blocksizes = [blocksize]
    tours = 1
    saturation_ratio = 0.5
    threads = threads
    gpus = 2
    dim4free = dim4free_wrapper(default_dim4free_fun,blocksize)
    extra_dim4free = 12
    pump_params = {"down_sieve": True,"saturation_error":"ignore"}#,"start_up_n": blocksize - dim4free}
    print("Test2. Start pump cost test with saturation_ratio = %.2f, d = %d" %(saturation_ratio,d))
    
    file_name = "[20221011]pump_cost_in_pnjbkz_and_progressive_sieve/test2/svp-%d-%d-xf=%d-f=%d.txt" % (d,blocksize,extra_dim4free,dim4free)
    
    f = open(file_name,"w")
    for blocksize in blocksizes:
        T_tmp = []
        RAM_tmp = []
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            # f.write("Starting a pnjBKZ-%d tour. " % (blocksize))
            params = SieverParams(threads=threads,gpus=gpus,saturation_ratio=saturation_ratio)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f, beta-f = %d" % (slope,blocksize-dim4free))
            # f.write("Intial Slope = %.5f \n" % slope)
            T_pumps,slopes,ghs = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=extra_dim4free,
                                            dim4free_fun=default_dim4free_fun,
                                            pump_params=pump_params)
            
            # f.write("Pump Cost in pnjBKZ: ")
            f.write(str(T_pumps))   
            f.write("\n")
            f.write(str(slopes))
            f.write("\n")
            f.write(str(ghs))
            T_BKZ = time.time()-T0_BKZ
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            # f.write(fmt % (slope, T_BKZ) +"\n")
    
        f.flush()

    f.close()

    
    
def test3(d,saturation_ratio,threads):
    tours = 1
    threads = threads
    gpus = 2
    # dim4free = dim4free_wrapper(default_dim4free_fun,d)
    pump_params = {"down_sieve": True,"saturation_error":"ignore"}#,"start_up_n": blocksize - dim4free}
    print("Test3. Test relation between saturation_ratio = %.2f (d = %d) and the lengths of short vectors." %(saturation_ratio,d))
    
    file_name = "[20221011]pump_cost_in_pnjbkz_and_progressive_sieve/test3/pump(%d,%d,%d)-%d.txt" % (0,d,0,int(saturation_ratio*100))
    f = open(file_name,"w")
    params = SieverParams(threads=threads,gpus=gpus,saturation_ratio=saturation_ratio)
    A, bkz = load_svpchallenge_and_randomize(d,seed = 0)
    g6k = Siever(A, params)
    slope = basis_quality(g6k.M)["/"]
    print("Intial Slope = %.5f." % (slope))
    
    print("Start sieve in dimension %d." % (d))
    g6k.shrink_db(0)
    g6k.lll(0, d)
    g6k.initialize_local(0, 0, d)

    #generate short vectors by sieving
    try:
        g6k()
    except SaturationError:
        pass
    print("Finish sieve.")
    print("-------------------------------------")
    print("Start recover short vectors...")
    db = list(g6k.itervalues())
    i = 0
    short_vectors = {}
    gh = gaussian_heuristic([g6k.M.get_r(_,_) for _ in range(g6k.M.B.ncols)])
    for x in db:
        v = A.multiply_left(x)
        l = sum(v_**2 for v_ in v)
        length = round(l/gh,3) 
        if length not in short_vectors:
            short_vectors[length] = 1
        else:
            short_vectors[length] += 1
        print("\r %d/%d" %(i,len(db)), end= '')
        i += 1
    print()
    print(short_vectors)
    f.write(str(short_vectors))
    f.write("\n")
    print("==========================================")
    



(d, blocksize,tours, saturation_ratio, threads,gpus) = (180,100,10,3,32,2)
for blocksize in range(80,90,10):
    # test2(d, blocksize, threads)
    test1(d, blocksize, saturation_ratio, threads)
# test2(d,blocksize,tours,saturation_ratio,threads,gpus)

# (d, threads) = (80,32)
# for saturation_ratio in [0.1*_ for _ in range(1,11)]:
    # test3(d,saturation_ratio,threads)