#!/usr/bin/env python
# -*- coding: utf-8 -*-
####
#
#   Copyright (C) 2018-2021 Team G6K
#
#   This file is part of G6K. G6K is free software:
#   you can redistribute it and/or modify it under the terms of the
#   GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version.
#
#   G6K is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with G6K. If not, see <http://www.gnu.org/licenses/>.
#
####


"""
LWE Challenge Solving Command Line Client
"""

from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time

from collections import OrderedDict # noqa
from math import log

from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality
from fpylll.util import gaussian_heuristic as gh
from numpy import True_

from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from g6k.algorithms.pump import pump
#from g6k.algorithms.pump_cpu import pump
from g6k.siever import Siever
from g6k.utils.cli import parse_args, run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from g6k.utils.util import load_lwe_challenge,load_lwe_challenge_mid

# from g6k.utils.lwe_estimation import gsa_params, primal_lattice_basis
from g6k.utils.lwe_estimation_old import gsa_params, primal_lattice_basis
# from six.moves import range

#from pro_pnj_bkz_optimization import *


def gen_original_gs(n,alpha, goal_margin):
    """
    Run the primal attack against Darmstadt LWE instance (n, alpha).

    :param n: the dimension of the LWE-challenge secret
    :param alpha: the noise rate of the LWE-challenge
    """

    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    
    #print("-------------------------")
    #print("Primal attack, LWE challenge n=%d, alpha=%.4f" % (n, alpha))
    m = None
    if m is None:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True)
            (b, s, m) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    else:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True)
            (b, s, _) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    #print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    #print()

    # no use in having a very small b
    b = max(b, s-65)
    
    
    B = primal_lattice_basis(A, c, q, m=m) #debug

    params = None
    g6k = Siever(B, params)
    #print("GSO precision: ", g6k.M.float_type)

    d = g6k.full_n

    g6k.lll(0, g6k.full_n)
    #slope = basis_quality(g6k.M)["/"]
    # print("Intial Slope = %.5f\n" % slope)

    log_rr0 = [log(g6k.M.get_r(i, i)) for i in range(d)]


    target_norm = goal_margin * (alpha*q)**2 * m + 1

    return b,s,target_norm,log_rr0,0
    
