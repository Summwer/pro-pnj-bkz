

from cysignals.signals cimport sig_on, sig_off
cimport numpy as np
#from decl cimport Params
from numpy import zeros, float64, int64

#function about EnumBS selection algorithm, to select an appropriate (blocksize,jump) strategy for reduction.
cdef class TEST(object):
    def __init__(self):
        self._core = new TEST_c() #params, 

    def print_info(self):
        self._core.print_info()