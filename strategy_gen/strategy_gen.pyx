

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
from decl cimport Params
from numpy import zeros, float64, int64

#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class EnumBS(object):
    def __init__(self, dim):
        cdef Params params
        params.method = 1
        self._core1 = new EnumBS_c(params, dim) #params, 
        #self._core1.print_param_setting()
        self.have_strategy_gen = True

    
    
    #def strategy_gen(self, dim, dvol):
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

    
    def __dealloc__(self):
        del self._core1
    
    

    

cdef class BSSA(object):
    def __init__(self, dim):
        cdef Params params
        params.method = 2
       
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
    