from g6k.utils.util import load_svpchallenge_and_randomize
from g6k.utils.util import load_lwe_challenge, load_matrix_file
from g6k.utils.lwe_estimation_old import gsa_params, primal_lattice_basis
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from g6k.siever import Siever
from g6k.siever_params import SieverParams
from fpylll.tools.quality import basis_quality
import time
from g6k.algorithms.pump import pump
from pro_pnjBKZ_simulator.codes.util import gaussian_heuristic, dim4free_wrapper, default_dim4free_fun,theo_dim4free_fun1

#T--(beta-f)
def pnjbkz_cost_test1():
    d = 180
    blocksizes = [_ for _ in range(120,160,10)]
    pump_params = {"down_sieve": True}
    tours = 2 #test times
    print("Start Pump Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("1. fixed d = %d, threads = %d, gpus = %d " %(d, 32, 2))
    file_name = "PnjBKZcost_test_result_test1.txt"
    f = open(file_name,"a")
    bkzCost = {}
    RAMCost = {}
    for blocksize in blocksizes:
        T_tmp = []
        RAM_tmp = []
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            params = SieverParams(threads=32,gpus=2)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f" % slope)
            RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            RAM_tmp.append(RAM_cost)
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[blocksize] = T_tmp
        RAMCost[blocksize] = RAM_tmp
        f.write("beta = %d, " %blocksize)    
        f.write("time cost(s): "+str(T_tmp)+",")
        f.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)
    print("RAM cost:")
    print(RAMCost)



#T--d
def pnjbkz_cost_test2():
    ds = [_ for _ in range(150,200,5)]
    # ds = [_ for _ in range(185,300,10)]
    blocksize = 90
    pump_params = {"down_sieve": True}
    tours = 2 #test times
    print("Start Pump Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("1. fixed beta = %d, threads = %d, gpus = %d " %(blocksize, 32, 2))
    file_name = "PnjBKZcost_test_result_test2_90.txt"
    f = open(file_name,"w")
    bkzCost = {}
    RAMCost = {}
    for d in ds:
        T_tmp = []
        RAM_tmp = []
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour for dim=%d. " % (blocksize,d))
            params = SieverParams(threads=32,gpus=2)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f" % slope)
            RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            RAM_tmp.append(RAM_cost)
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[d] = T_tmp
        RAMCost[d] = RAM_tmp
        f.write("d = %d, " %d)    
        f.write("time cost(s): "+str(T_tmp)+",")
        f.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)
    print("RAM cost:")
    print(RAMCost)


#T--d
def pnjbkz_cost_test2_lwe(n,alpha,blocksize,gap,ds):
    
    # blocksize = 80
    pump_params = {"down_sieve": True}
    tours = 1 #test times
    print("Start Pump Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("1. fixed beta = %d, threads = %d, gpus = %d " %(blocksize, 32, 2))
    file_name = "PnjBKZcost_test_result_test2_lwe_80_005_%d.txt" %blocksize
    # file_name = "PnjBKZcost_test_result_test2_lwe_80_test_pnj_BKZ_one_pump.txt"
    f = open(file_name,"w")
    bkzCost = {}
    RAMCost = {}
    for d in ds:
        T_tmp = []
        RAM_tmp = []
        params = SieverParams(threads=32,gpus=2)
        A, c, q = load_lwe_challenge(n=n, alpha=alpha)
        m = d - 1
        B = primal_lattice_basis(A, c, q, m=m)
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
        
            g6k = Siever(B, params)
            slope = basis_quality(g6k.M)["/"]
            print("Starting a pnjBKZ-%d tour for dim=%d. " % (blocksize,B.nrows))
            print("Intial Slope = %.5f" % slope)
            RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            RAM_tmp.append(RAM_cost)
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[d] = T_tmp
        RAMCost[d] = RAM_tmp
        f.write("d = %d, " %d)    
        f.write("time cost(s): "+str(T_tmp)+",")
        f.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)
    print("RAM cost:")
    print(RAMCost)




#T--basis quality
def pnjbkz_cost_test3():
    d = 180
    blocksizes = [_ for _ in range(100,120,10)]
    pump_params = {"down_sieve": True}
    tours = 10 #test times
    print("Start Pump Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("Test3. fixed d = %d, threads = %d, gpus = %d " %(d, 32, 2))
    file_name = "PnjBKZcost_test_result_test3.txt"
    f = open(file_name,"a")
    bkzCost = {}
    RAMCost = {}
    for blocksize in blocksizes:
        T_tmp = []
        RAM_tmp = []
        A, bkz = load_svpchallenge_and_randomize(d)
        params = SieverParams(threads=32,gpus=2)
        g6k = Siever(A, params)
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f" % slope)
            RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            RAM_tmp.append(RAM_cost)
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[blocksize] = T_tmp
        RAMCost[blocksize] = RAM_tmp
        f.write("beta = %d, " %blocksize)    
        f.write("time cost(s): "+str(T_tmp)+",")
        f.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)
    print("RAM cost:")
    print(RAMCost)



#T-- down_sieve = False
def pnjbkz_cost_test4():
    d = 180
    blocksizes = [_ for _ in range(80,120,10)]
    pump_params = {"down_sieve": False}
    tours = 10 #test times
    print("Start Pump Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("Test4. fixed d = %d, threads = %d, gpus = %d " %(d, 32, 2))
    file_name = "PnjBKZcost_test_result_test4.txt"
    f = open(file_name,"a")
    bkzCost = {}
    RAMCost = {}
    for blocksize in blocksizes:
        T_tmp = []
        RAM_tmp = []
        for _ in range(tours):
            print()
            T0_BKZ = time.time()
            print("Starting a pnjBKZ-%d tour. " % (blocksize))
            params = SieverParams(threads=32,gpus=2)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f" % slope)
            RAM_cost = pump_n_jump_bkz_tour(g6k, dummy_tracer, blocksize, jump=1,
                                                verbose=True,
                                            extra_dim4free=12,
                                            dim4free_fun="default_dim4free_fun",
                                            pump_params=pump_params)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            RAM_tmp.append(RAM_cost)
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[blocksize] = T_tmp
        RAMCost[blocksize] = RAM_tmp
        f.write("beta = %d, " %blocksize)    
        f.write("time cost(s): "+str(T_tmp)+",")
        f.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)
    print("RAM cost:")
    print(RAMCost)



#T--(beta-f), plain bkz2.0
def pnjbkz_cost_test5():
    d = 180
    blocksizes = [_ for _ in range(60,70,1)]
    pump_params = {"down_sieve": True}
    tours = 5 #test times
    print("Start Enum Cost Test Experiment, test %d times for each blocksizes" %(tours))
    print("Test5. fixed d = %d, threads = %d, gpus = %d " %(d, 32, 2))
    file_name = "PnjBKZcost_test_result_test5.txt"
    f = open(file_name,"w")
    bkzCost = {}
    for blocksize in blocksizes:
        T_tmp = []
        for _ in range(tours):
            T0_BKZ = time.time()
            print("Starting a BKZ-%d tour. " % (blocksize))
            params = SieverParams(threads=32,gpus=2)
            A, bkz = load_svpchallenge_and_randomize(d)
            g6k = Siever(A, params)
            slope = basis_quality(g6k.M)["/"]
            
            print("Intial Slope = %.5f" % slope)
            bkz = BKZReduction(g6k.M)
            par = fplll_bkz.Param(blocksize,strategies=fplll_bkz.DEFAULT_STRATEGY,max_loops=1)
            bkz(par)
            
            T_BKZ = time.time() - T0_BKZ
            T_tmp.append(T_BKZ) # each time cost 
            
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, T_BKZ))
            
        bkzCost[blocksize] = T_tmp
        f.write("beta = %d, " %blocksize)    
        f.write("time cost(s): "+str(T_tmp)+".\n")
        f.flush()

    f.close()

    print("pnjBKZ cost:")
    print(bkzCost)


def gen_pump_data():
    #same d
    d = 180
    tours = 2
    # pump_params = {"down_sieve": True}
    print("Test6. Start Pump Simulator Test " )
    # file_name = "PnjBKZcost_test_result_test6.txt"
    # file_name2 = "PnjBKZcost_test_result_test6_pumpData.txt"
    file_name = "PnjBKZcost_test_result_test6_g6k_down_sieve.txt"
    file_name2 = "PnjBKZcost_test_result_test6_pumpData_g6k_down_sieve.txt"
    f1 = open(file_name,"w")
    f2 = open(file_name2,"w")
    #simulate pump_estimation2
    pumpCost = {}
    RAMCost = {}
    for llb in range(d-50,0,-5):
        T_tmp = []
        RAM_tmp = []
        beta = d - llb
        # f = dim4free_wrapper(theo_dim4free_fun1, beta)
        f = dim4free_wrapper(default_dim4free_fun, beta)

        params = SieverParams(threads=32,gpus=2)
        A, bkz = load_svpchallenge_and_randomize(d)
        g6k = Siever(A, params)
        slope = basis_quality(g6k.M)["/"]
        print("Intial Slope = %.5f" % slope)
                
        T0_pump = time.time()
        f2.write(str([g6k.M.get_r(i, i) for i in range(d)]))
        if beta - f < 131:
            for _ in range(tours):
                print()
                print( "Starting svp pump_{%d, %d, %d}" % (llb, d-llb, f) ) # noqa
                
                _,RAM_pump = pump(g6k, dummy_tracer, llb, d-llb, f, verbose=True,down_sieve=True)       
                f2.write(str([g6k.M.get_r(i, i) for i in range(d)]))
                T_pump = time.time() - T0_pump
                T_tmp.append(T_pump) # each time cost 
                RAM_tmp.append(RAM_pump)   
                slope = basis_quality(g6k.M)["/"]
                fmt = "slope: %.5f, walltime: %.3f sec, max_RAM_cost: %.3f GB"
                print(fmt % (slope, T_pump,RAM_pump))
                
            pumpCost[beta] = T_tmp
            RAMCost[beta] = RAM_tmp
            f1.write("beta = %d, " %beta)    
            f1.write("time cost(s): "+str(T_tmp)+".\n")
            f1.write("RAM cost(GB): "+str(RAM_tmp)+".\n")
            f1.flush()

    f1.close()
    f2.close()

    print("Pump cost:")
    print(pumpCost)
    print("RAM cost:")
    print(RAMCost)




# pnjbkz_cost_test1() #test1: T--(beta-f)
# pnjbkz_cost_test2() #test2: T--d
# gen_pump_data()
# pnjbkz_cost_test2_lwe(80,0.005,0) #test2: T--d
# pnjbkz_cost_test2_lwe(80,0.005,90) #test2: T--d
# pnjbkz_cost_test2_lwe(80,0.005,70,10) #test2: T--d
gap = 20
ds = list(range(220,270,gap))
pnjbkz_cost_test2_lwe(80,0.005,110,20,ds) #test2: T--d
ds = list(range(160,270,gap))
pnjbkz_cost_test2_lwe(80,0.005,120,30,ds) #test2: T--d
# pnjbkz_cost_test2_lwe(80,0.005,110) #test2: T--d
# pnjbkz_cost_test2_lwe(80,0.005,180) #test2: T--d
# pnjbkz_cost_test3() #test3: T--basis quality
# pnjbkz_cost_test4() #test4: T--down_sieve = False
# pnjbkz_cost_test5() #test5: T--(beta-f), plain bkz2.0
# pnjbkz_cost_test1() #test1: T--(beta-f) 120-160
# gen_pump_data()