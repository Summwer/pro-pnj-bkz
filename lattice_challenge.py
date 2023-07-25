#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lattice Challenge Solver Command Line Client (SIS problem)
https://www.latticechallenge.org/
"""

import copy,sys
import logging
import pickle as pickler
from collections import OrderedDict

from fpylll.util import gaussian_heuristic

from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k.utils.cli import parse_args, run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import load_latticechallenge, load_matrix_file, db_stats
from fpylll import BKZ as BKZ_FPYLLL
from fpylll.tools.bkz_stats import dummy_tracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour
import time
from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality

def sis_kernel(arg0, params=None, seed=None):
    logger = logging.getLogger('sis')

    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)
    load_matrix = params.pop("load_matrix")
    pre_bkz = params.pop("pre_bkz")
    pump_params = pop_prefixed_params("pump", params)
    workout_params = pop_prefixed_params("workout", params)
    verbose = params.pop("verbose")
    if verbose:
        workout_params["verbose"] = True
    challenge_seed = params.pop("challenge_seed")
    high_prec = params.pop("high_prec")
    trace = params.pop("trace")
    dim4free_fun = params.pop("bkz/dim4free_fun")


    A, goal_r0, q = load_latticechallenge(n)
    if verbose:
        print("Loaded challenge dim=%d, goal_r0=%d, q=%d" %(n,goal_r0,q))

    g6k = Siever(A, params)

    if trace:
        tracer = SieveTreeTracer(g6k, root_label=("lattice-challenge", n), start_clocks=True)
    else:
        tracer = dummy_tracer

    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

    if verbose:
        print("gh = %f, goal_r0/gh = %f, r0/gh = %f" % (gh, goal_r0/gh, sum([x*x for x in A[0]])/gh))

    #flast = workout(g6k, tracer, 0, n, goal_r0=goal_r0,
     #               pump_params=pump_params, **workout_params)
    
    blocksizes = [50,60,70]
    
    for blocksize in blocksizes:
        T0 = time.time()
            
        # BKZ tours

        if blocksize < 100:
            if verbose:
                print("Starting a fpylll BKZ-%d tour. " % (blocksize), end='')
                sys.stdout.flush()
            bkz = BKZReduction(g6k.M)
            par = fplll_bkz.Param(blocksize,
                                          strategies=fplll_bkz.DEFAULT_STRATEGY,
                                          max_loops=1)
            bkz(par)

        else:
            if verbose:
                print("Starting a pnjBKZ-%d tour. " % (blocksize))

            pump_n_jump_bkz_tour(g6k, tracer, blocksize, jump=1,
                                         verbose=verbose,
                        
                                         #dim4free_fun=dim4free_fun,
                                         goal_r0=goal_r0,
                                         pump_params=pump_params)

        if verbose:
            slope = basis_quality(g6k.M)["/"]
            fmt = "slope: %.5f, walltime: %.3f sec"
            print(fmt % (slope, time.time() - T0))
                
            norm = sum([x*x for x in A[0]])
            logger.info("norm %.1f ,hf %.5f" % (norm**.5, (norm/gh)**.5))

        g6k.lll(0, g6k.full_n)


        if g6k.M.get_r(0, 0) <= goal_r0:
            break

    
    if verbose:
        logger.info("sol %d, %s" % (n, A[0]))

    


def sis():
    """
    Run a Workout until 1.05-approx-SVP on matrices with dimensions in ``range(lower_bound, upper_bound, step_size)``.
    """
    description = sis.__doc__

    args, all_params = parse_args(description,
                                  load_matrix=None,
                                  pre_bkz=None,
                                  verbose=True,
                                  challenge_seed=0,
                                  workout__dim4free_dec=2,
                                  bkz__dim4free_fun="default_dim4free_fun",
                                  trace=True)

    stats = run_all(sis_kernel, all_params.values(),
                    lower_bound=args.lower_bound,
                    upper_bound=args.upper_bound,
                    step_size=args.step_size,
                    trials=args.trials,
                    workers=args.workers,
                    seed=args.seed)

    inverse_all_params = OrderedDict([(v, k) for (k, v) in all_params.items()])

    for (n, params) in stats:
        stat = stats[(n, params)]
        if stat[0] is None:
            logging.info("Trace disabled")
            continue

        if len(stat) > 0:
            cputime = sum([float(node["cputime"]) for node in stat])/len(stat)
            walltime = sum([float(node["walltime"]) for node in stat])/len(stat)
            flast = sum([float(node["flast"]) for node in stat])/len(stat)
            avr_db, max_db = db_stats(stat)
            fmt = "%48s :: m: %1d, n: %2d, cputime :%7.4fs, walltime :%7.4fs, flast : %2.2f, avr_max db: 2^%2.2f, max_max db: 2^%2.2f" # noqa
            logging.info(fmt % (inverse_all_params[params], params.threads, n, cputime, walltime, flast, avr_db, max_db))
        else:
            logging.info("Trace disabled")

    if args.pickle:
        pickler.dump(stats, open("hkz-sis-%d-%d-%d-%d.sobj" %
                                 (args.lower_bound, args.upper_bound, args.step_size, args.trials), "wb"))


if __name__ == '__main__':
    sis()
