# -*- coding: utf-8 -*-
"""
enumbs class
"""
from decl cimport EnumBS as EnumBS_c


cdef class EnumBS(object):
    cdef EnumBS_c *_core
    cdef object have_strategy_gen
    

