# -*- coding: utf-8 -*-
"""
cost class
"""
from decl cimport COST as COST_c


cdef class COST(object):
    cdef COST_c *_core
   