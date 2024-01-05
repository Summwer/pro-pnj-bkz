#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import subprocess

# compile framework/libest.so if not done already
subprocess.check_call("make -C framework",shell=True)

# read actual values of all build variables from framework/Makefile
makefile_defs = subprocess.getoutput("make -C framework printvariables | grep '='").splitlines()

def read_from_makefile(field):
    global makefile_defs
    data = [line for line in makefile_defs if line.startswith(field)][0]
    data = '=' .join(data.split('=')[1:])
    data = data.strip()
    data = [arg for arg in data.split(' ') if arg.strip()]
    return data

extra_compile_args = read_from_makefile("CXXFLAGS")
extra_link_args = read_from_makefile("LDFLAGS") + read_from_makefile("LIBADD")

kwds = {
    "language": "c++",
    "extra_compile_args": extra_compile_args,
    "extra_link_args": extra_link_args,
    "libraries": ["gmp", "pthread", "est"],
    "include_dirs": [numpy.get_include()]
    }

extensions = [
    Extension("enumbs", ["enumbs.pyx"], **kwds),
    # Extension("g6k.siever_params", ["g6k/siever_param.pyx"], **kwds)
]

setup(
    name="enumbs",
    version="0.0.1",
    ext_modules=cythonize(extensions, compiler_directives={'binding': True,
                                                           'embedsignature': True,
                                                           'language_level': 2}),
    packages=[],
)
