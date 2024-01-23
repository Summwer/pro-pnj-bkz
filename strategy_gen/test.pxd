# -*- coding: utf-8 -*-
"""
cost class
"""
from decl cimport TEST as TEST_c


cdef class TEST(object):
    cdef TEST_c *_core
   