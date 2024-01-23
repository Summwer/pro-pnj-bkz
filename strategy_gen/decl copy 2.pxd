# -*- coding: utf-8 -*-




cdef extern from "framework/cost.h" nogil:
    cdef cppclass COST:
        COST()
        double practical_bkz_cost_dd(int d,int beta,int jump);
        