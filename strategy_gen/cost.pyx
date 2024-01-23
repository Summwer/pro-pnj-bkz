

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
#from decl cimport Params
from numpy import zeros, float64, int64

#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class COST(object):
    def __init__(self):
        self._core = new COST_c() #params, 

    def practical_bkz_cost_dd(self,d,beta,jump):
        self._core.practical_bkz_cost_dd( d,beta, jump)
       