

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
from decl cimport Params, test_lwechal_from_gsa, test_lwechal_from_actual_l
from numpy import zeros, float64, int64

#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class EnumBS(object):


    def __init__(self, dim, float_type):
        cdef Params params
        
        params.method = 1

        if(float_type == "dd"):
            params.cost_model = 3
            params.sim_d4f = 2
        else:
            params.cost_model = 2
            params.sim_d4f = 1
        self._core1 = new EnumBS_c(params, dim) #params, 
        #self._core1.print_param_setting()
        self.have_strategy_gen = True

    
    
    #def __call__(self, dim, dvol):
    def __call__(self, l):
        """
        Call EnumBS strategy selction algorithm 
        Input:  -- dim: dimension,
                -- dvol: log(det(lattice)) (dim,dvol) will generate the normlized gs-lengths [2*log2(gs-lengths)], (gs-lengths - sigma)
        Output: blocksize strategy 

        """

        dim = len(l)
        cdef np.ndarray l0 = zeros(dim, dtype=float64)
        for i in range(dim):
            l0[i] = l[i]

        sig_on()
        self._core1.enumbs_est_in_parallel(<double*> l0.data);
        #self._core1.enumbs_est_in_parallel(dim, dvol);
        sig_off()
        self.have_strategy_gen = True
    

    def get_strategy(self):
    
        assert(self.have_strategy_gen)
        strategy_size = self._core1.strategy_size()
     
        cdef np.ndarray strategy = zeros((strategy_size,3), dtype=int64)

        sig_on()
        self._core1.get_strategy(<long*> strategy.data)
        sig_off()
        
        S = []
        for i in range(strategy_size):
            S.append(tuple(strategy[i]))
    
        return S
    
    def get_target_slope(self):
        return self._core1.get_target_slope()

    
    def __dealloc__(self):
        del self._core1
    
    

    

cdef class BSSA(object):
    def __init__(self, dim, version,float_type):
        cdef Params params
        params.method = 2
        if(version == "v2"):
            params.bssa_tradition = 0
        if(version == "v1"):
            params.bssa_tradition = 1
        
        if(float_type == "dd"):
            params.cost_model = 3
            params.sim_d4f = 2
        else:
            params.cost_model = 2
            params.sim_d4f = 1
        self._core2 = new BSSA_c(params,dim) #params, 
        #self._core2.print_param_setting()
        self.have_strategy_gen = False

    
    
    #def strategy_gen(self, dim, dvol):
    def __call__(self, l, sbeta =None, gbeta = None):
        """
        Call BSSA strategy selction algorithm 
        Input:  -- dim: dimension,
                -- dvol: log(det(lattice)) (dim,dvol) will generate the normlized gs-lengths [2*log2(gs-lengths)], (gs-lengths - sigma)
        Output: blocksize strategy 

        """

        dim = len(l)
        cdef np.ndarray l0 = zeros(dim, dtype=float64)
        for i in range(dim):
            l0[i] = l[i]

        if sbeta is None:
            sbeta = 50
        if gbeta is None:
            gbeta = dim

        sig_on()
        self._core2.bssa_est(<double*> l0.data, sbeta, gbeta)
        sig_off()
        self.have_strategy_gen = True
    
    
    def get_strategy(self):
    
        assert(self.have_strategy_gen)
        strategy_size = self._core2.strategy_size()
     
        cdef np.ndarray strategy = zeros((strategy_size,3), dtype=int64)

        sig_on()
        self._core2.get_strategy(<long*> strategy.data)
        sig_off()
        
        S = []
        for i in range(strategy_size):
            S.append(tuple(strategy[i]))
    
        return S

    

    
    
    def __dealloc__(self):
        del self._core2
    
def lwechal_simulation_gsa(dim, dvol, S, strategy_size,float_type = None):
    cdef Params params
    params.verbose = True
    if(float_type == "dd"):
        params.cost_model = 3
        params.sim_d4f = 2
    else:
        params.cost_model = 2
        params.sim_d4f = 1

    cdef np.ndarray strategy = zeros((strategy_size,3), dtype=int64)

    for i in range(strategy_size):
        strategy[i][0] = S[i][0]
        strategy[i][1] = S[i][1]
        strategy[i][2] = S[i][2]
    


    sig_on()
    test_lwechal_from_gsa(params, dim, dvol, <long*> strategy.data, strategy_size)
    sig_off()



def lwechal_simulation_actual_l(l, S, strategy_size,float_type = None):
    cdef Params params
    params.verbose = True
    if(float_type == "dd"):
        params.cost_model = 3
        params.sim_d4f = 2
    else:
        params.cost_model = 2
        params.sim_d4f = 1

    dim = len(l)
    cdef np.ndarray l0 = zeros(dim, dtype=float64)
    for i in range(dim):
        l0[i] = l[i]

    cdef np.ndarray strategy = zeros((strategy_size,3), dtype=int64)

    for i in range(strategy_size):
        strategy[i][0] = S[i][0]
        strategy[i][1] = S[i][1]
        strategy[i][2] = S[i][2]

    sig_on()
    test_lwechal_from_actual_l(params, <long*> strategy.data, strategy_size, <double*> l0.data, dim)
    sig_off()