# -*- coding: utf-8 -*-




cdef extern from "framework/test.h" nogil:
    cdef cppclass TEST:
        TEST()
        void print_info();
