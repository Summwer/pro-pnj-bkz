

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
from decl cimport Params
from numpy import zeros, float64, int64


#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class EnumBS(object):
    def __init__(self, dim):
        cdef Params params
        self._core = new EnumBS_c(params, dim)

        #self._core.print_param_setting()
    
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
        self._core.enumbs_est_in_parallel(<double*> l0.data);
        #self._core.enumbs_est_in_parallel(dim, dvol);
        sig_off()
        self.have_strategy_gen = True
    
    def get_strategy(self):
    
        assert(self.have_strategy_gen)
        strategy_size = self._core.strategy_size()
     
        cdef np.ndarray strategy = zeros((strategy_size,3), dtype=int64)

        sig_on()
        self._core.get_strategy(<long*> strategy.data)
        sig_off()
        
        #S = []
        #for i in range(strategy_size):
        #    S.append(tuple(strategy[i]))
    
        #return S

    def __dealloc__(self):
        del self._core
    