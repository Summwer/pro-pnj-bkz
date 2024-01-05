

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
from decl cimport Params
from numpy import zeros, float64

#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class EnumBS(object):
    def __init__(self, dim):
        cdef Params params
        self._core = new EnumBS_c(&params, dim)

    def strategy_gen(self, l):
        """
        Call EnumBS strategy selction algorithm 
        Input: [2*log2(gs-lengths)]
        Output: blocksize strategy 

        """
        dim = len(l)
        cdef np.ndarray l0 = zeros(dim, dtype=float64)
        for i in range(dim):
            l0[i] = l[i]

        sig_on()
        self._core.enumbs_est_in_parallel(<double*>l0.data);
        sig_off()
        
        return 