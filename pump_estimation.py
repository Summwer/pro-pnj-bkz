# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function
from math import pi,exp,log,sqrt,lgamma
from math import pi,e, lgamma, log, pi, log2
from scipy.special import chdtr 
from fpylll.util import gaussian_heuristic
from g6k.utils.cost import  dim4free_wrapper, theo_dim4free_fun2
from cost import get_pump_time

#pump estimation with chi-squared probability
def pump_estimation(rr,q, alpha, succ_prob = 0.999):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param rr: vector of squared Gram-Schmidt norms
    :param q: LWE param, to compute sigma 
    :param alpha: LWE param, to compute sigma, alpha*q = sigma

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(rr)
    PSC = 0.
    pre_psvp = 0.
    for beta in range(30,d):
        GH = gaussian_heuristic(rr[d-beta:])
        length=(GH/(sigma**2))
        psvp = chdtr(beta, length)
        if(pre_psvp >= psvp):
            continue
        if(psvp >= succ_prob):
            break
        f = dim4free_wrapper(theo_dim4free_fun2, beta)
        PSC +=  get_pump_time(beta-f, d) * (psvp-pre_psvp)
        pre_psvp = psvp
    return log2(PSC), beta



