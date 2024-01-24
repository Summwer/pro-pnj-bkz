# -*- coding: utf-8 -*-
"""
enumbs class & BSSA class
"""
from decl cimport EnumBS as EnumBS_c
from decl cimport BSSA as BSSA_c

cdef class EnumBS(object):
    cdef EnumBS_c *_core1
    cdef object have_strategy_gen
    

cdef class BSSA(object):
    cdef BSSA_c *_core2
    cdef object have_strategy_gen
    