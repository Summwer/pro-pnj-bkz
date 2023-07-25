# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time
import math
import matplotlib.pyplot as plt
import os

from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma
from math import pi,e, lgamma, log, pi
from .gs_length_simulator_v2 import calculate_sim_log_gs_lengths_test as sim_log_gs_lengths_above_45
from .gs_length_simulator_v2 import calculate_sim_log_gs_lengths_test1 as sim_log_gs_lengths_below_45
from .util import default_dim4free_fun,theo_dim4free_fun1, theo_dim4free_fun2, dim4free_wrapper, get_pre_pnj_time, get_current_slope, predict_slope_with_beta,get_beta_from_sieve_dim

             
def log_gs_lengths_simulator(log_GS_lengths, beta, loop, d, jump):
    if beta > 45:
        return sim_log_gs_lengths_above_45(log_GS_lengths, beta, loop, d, jump)
    else: 
        return sim_log_gs_lengths_below_45(log_GS_lengths, beta, loop, d)




#Simulate pnj-bkz(beta) for one round
def simulatepnjbkzloop(log_rr0, beta, loop, d, jump):
    # file.write(beta)
    # extra_dim4free = 12
    extra_dim4free = 12
    f = dim4free_wrapper(default_dim4free_fun, beta)
    if jump <=2:
        beta_ = beta
    elif jump>=3 and jump <=4:
        beta_ = get_beta_from_sieve_dim(beta-f,d,theo_dim4free_fun2)
    elif jump>=5:
        beta_ = get_beta_from_sieve_dim(beta-f,d,theo_dim4free_fun1)

    sim_log_gs_lengths = log_gs_lengths_simulator(log_rr0, beta_, loop, d, jump) 
    # file.write(beta, loop,get_current_slope(sim_log_gs_lengths,0,d))

    beta = beta + extra_dim4free
    
    f = f + extra_dim4free
    
    
    #pnj-BKZ time
    pre_pnj_time = loop * get_pre_pnj_time(d, beta, f, jump)
    

    return sim_log_gs_lengths,pre_pnj_time


#Simulate pnj-bkz(beta) for sufficient rounds = loop-1
def PPNJBKZSimulate(log_rr0,beta,d,jump):
    # square_error = 1*len(log_rr0)
    loop = 0
    slope  = get_current_slope(log_rr0,0,d)
    slope1 = slope - 1
    
    # if beta<50:
        # tem_condition = ((abs(slope1) - abs(slope))>0.001)
    # elif beta >=50:
    predict_slope = predict_slope_with_beta(beta)
    tem_condition = (((abs(slope1) - abs(slope))>0.001) and (abs(slope)>abs(predict_slope)))
    
    while(tem_condition): # and 
        log_rr1 =  log_gs_lengths_simulator(log_rr0, beta, loop, d, jump) 
        log_rr = log_gs_lengths_simulator(log_rr0, beta, loop+1, d, jump) 
        loop+=1
        
        slope1 = get_current_slope(log_rr1,0,d)
        slope = get_current_slope(log_rr,0,d)
        # file.write(beta,slope1,slope,loop)
        
        if beta<50:
            tem_condition = ((abs(slope1) - abs(slope))>0.001)
        elif beta >=50:
            predict_slope = predict_slope_with_beta(beta)
            tem_condition = (((abs(slope1) - abs(slope))>0.001) and (abs(slope)>abs(predict_slope)))
    # raise TypeError("Debug")
    return log_rr,loop,slope



