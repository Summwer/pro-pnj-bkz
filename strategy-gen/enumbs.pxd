# -*- coding: utf-8 -*-
"""
enumbs class
"""
from decl cimport EnumBS as EnumBS_c

cdef class EnumBS(object):
    cdef EnumBS_c *_core
    #cdef public object M
    #cdef SieverParams _params
    #cdef object initialized
    #cdef object dual_hash_l
