#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SVP Challenge Solver Command Line Client
"""

import copy
import logging
import pickle as pickler
from collections import OrderedDict
from fpylll.util import gaussian_heuristic
from g6k.algorithms.workout import workout
from g6k.siever import Siever
from g6k.utils.cli import parse_args, run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer
from g6k.utils.util import load_svpchallenge_and_randomize, load_matrix_file, db_stats
from fpylll import BKZ as BKZ_FPYLLL
from fpylll.tools.bkz_stats import dummy_tracer

def asvp_kernel(arg0, params=None, seed=None):
    logger = logging.getLogger('asvp')

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

    if load_matrix is None:
        A, bkz = load_svpchallenge_and_randomize(n, s=challenge_seed, seed=seed)
        if verbose:
            print("Loaded challenge dim %d" % n)
        if pre_bkz is not None:
            par = BKZ_FPYLLL.Param(pre_bkz, strategies=BKZ_FPYLLL.DEFAULT_STRATEGY, max_loops=1)
            bkz(par)

    else:
        A, _ = load_matrix_file(load_matrix, doLLL=False, high_prec=high_prec)
        if verbose:
            print("Loaded file '%s'" % load_matrix)


    g6k = Siever(A, params, seed=seed)

    if trace:
        tracer = SieveTreeTracer(g6k, root_label=("svp-challenge", n), start_clocks=True)
    else:
        tracer = dummy_tracer

    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

    gamma = params.pop("gamma")

    if gamma is None:
        goal_r0 = (1.05**2) * gh
    else:
        goal_r0 = gamma**2 * gh

    if verbose:
        print("gh = %f, goal_r0/gh = %f, r0/gh = %f" % (gh, goal_r0/gh, sum([x*x for x in A[0]])/gh))

    flast = workout(g6k, tracer, 0, n, goal_r0=goal_r0,
                    pump_params=pump_params, **workout_params)
    if verbose:
        logger.info("sol %d, %s" % (n, A[0]))

    norm = sum([x*x for x in A[0]])
    if verbose:
        logger.info("norm %.1f ,hf %.5f" % (norm**.5, (norm/gh)**.5))

    tracer.exit()


    if hasattr(tracer, "trace"):
        stat = tracer.trace
        stat.data["flast"] = flast
        return stat
    else:
        return None


def asvp():
    """
    Run a Workout until 1.05-approx-SVP on matrices with dimensions in ``range(lower_bound, upper_bound, step_size)``.
    """
    description = asvp.__doc__

    args, all_params = parse_args(description,
                                  load_matrix=None,
                                  pre_bkz=None,
                                  verbose=True,
                                  challenge_seed=0,
                                  workout__dim4free_dec=2,
                                  trace=True)

    stats = run_all(asvp_kernel, all_params.values(),
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
        pickler.dump(stats, open("hkz-asvp-%d-%d-%d-%d.sobj" %
                                 (args.lower_bound, args.upper_bound, args.step_size, args.trials), "wb"))


if __name__ == '__main__':
    asvp()
