#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time
import math
import matplotlib.pyplot as plt
import os
import shutil

from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma
from .util import  get_current_slope,show_pump_gs_figure,dim4free_wrapper,default_dim4free_fun,theo_dim4free_fun1
from math import sqrt



'''v2: Add the jump value into the simulator'''

#read GS length 
def read_gs_lengths(filename):
    f = open(filename,'r')
    #f = open("rr_set_bkz_50.txt",'r')
    GS_lengths = f.read()
    GS_lengths = GS_lengths.split(']')
    log_GS_lengths = []
   
    for i in range(len(GS_lengths)):
        gs = GS_lengths[i]
        if gs != '':
            gs = gs.replace('[','')
            gs = gs.split(', ')
            GS_lengths[i] = [float(_) for _ in gs]
            log_GS_lengths.append([log(_) for _ in GS_lengths[i]])
    return log_GS_lengths,GS_lengths

def read_blocksizes(file_name):
    f = open(file_name,'r')
    #f = open("50-0010/slopes-%s-50-0010.txt" %(pnj),'r') 
    blocksizes = f.read()
    f.close()
    blocksizes = blocksizes.split(' ')
    blocksizes.remove('')
    blocksizes = [int(_) for _ in blocksizes]
    return blocksizes

def read_times(file_name):
    f = open(file_name,'r')
    #f = open("50-0010/slopes-%s-50-0010.txt" %(pnj),'r') 
    times = f.read()
    f.close()
    times = times.split(' ')
    times.remove('')
    times = [float(_) for _ in times]
    return times

def read_file(file_name):
    #n=50,alpha=0.01,k=1,f=0,beta=10-66
    
    log_GS_lengths,GS_lengths = read_gs_lengths(file_name)
    
    return log_GS_lengths,GS_lengths,file_name


#beta>45
#pump has a scoring down equation, where score = ln(||b_i^*||theta^i*||pi(v_i)|| ), theta = 1.04
def calculate_sim_log_gs_lengths_for_pump(log_rr0,d,llb,f):
    
    # print(db_size_base)
    N=1
    r = d
    beta = r - llb
    #dim4f = dim4free_wrapper(theo_dim4free_fun1, beta)

    db_size_base = 0 #log(sqrt(4/3. * (beta / (beta-dim4f))))
    # print(db_size_base)
    if beta == 0:
        return log_rr0
    l = []
    for i in range(d):
        temp = log_rr0[i]/2
        l.append(temp)
    #print(len(l))
    l_ = []
    
    randomSquaredAverages = [2.95747208 , 2.905419231 , 2.808209731 , 2.704122291 , 2.589140986 , 2.473617892 , 2.363752978 , 2.254288677 , 2.160395606 , 2.056417068 , 1.962571394 , 1.86329581 , 1.781471788 , 1.691206302 , 1.607855519 , 1.526258286 , 1.45048576 , 1.377169955 , 1.304800173 , 1.238714637 , 1.172164783 , 1.11055388 , 1.04872671 , 0.9933346386 , 0.9385826564 , 0.8858378414 , 0.8346947761 , 0.7869433841 , 0.743196743 , 0.6999110693 , 0.6612589966 , 0.620890347 , 0.5840212199 , 0.5494461008 , 0.5153120689 , 0.4851862344 , 0.4542255962 , 0.426615261 , 0.4011004324 , 0.3761683243 , 0.3521434119 , 0.3305054658 , 0.3109962405 , 0.2976205498 , 0.2897254705]
    #print(len(randomSquaredAverages))
    rk = []
    for i in range(len(randomSquaredAverages)):
        RK = log(math.sqrt(randomSquaredAverages[i]))
        rk.append(RK)
    #print(rk)
    cd = []
    for i in range(beta):
        CD = lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2
        cd.append(CD)
    #print(cd)
    #print(N)
   
    for tours in range(N):
        flag = True
        for k in range(d - beta):
            l_.append(l[k])

        for k in range(d - beta ,d - 45):
            beta_ = d - k
            f = d
            sumf = 0
            sumk = 0
            for i in range(f):
                sumf = sumf + l[i]
            for i in range(k):
                sumk = sumk + l_[i]
            logV = sumf - sumk
            if flag == True:
                if logV / beta_ + cd[beta_-1] < l[k]:
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
                    flag = False
                else:
                    l_k = l[k]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                l_k = db_size_base + logV / beta_ + cd[beta_-1]
                if tours == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
                
            #print(l_)
        k_ = min (len(randomSquaredAverages), beta)
        #last 45 norms
        sumf=0
        sumk=0
        for i in range(d):
            sumf = sumf + l[i]
        for i in range(d - 45):
            sumk = sumk + l_[i]
        logV = sumf - sumk
        for k in range(d - k_, d):
            l_k =  db_size_base + logV/k_ + rk[k + 45 - d]
            if tours == 0:
                l_.append(l_k)
            else:
                l_[k] = l_k
                #Set l[i] as l_[i] 
        for i in range(d):
            l[i] = l_[i]
    return [_*2 for _ in l_]


def compute_square_error(list1,list2,flag = 1):
    square_error = 0
    l = min(len(list1),len(list2))
    L = max(len(list1),len(list2))
    for i in range(L):

        if i == 0 and list2[1]-list2[0] > 1 and flag == 1:
            #remove the error point.
            continue
        elif i < l :
            square_error += (list1[i]-list2[i])**2
        elif i>=l:
            try:
                square_error += list1[i]**2
            except IndexError:
                square_error += list2[i]**2
            
    return square_error/sum(_**2 for _ in list1)
             



def gen_pump_figures(dir):
    file_name = "../../PnjBKZcost_test_result_test6_pumpData_theo_down_sieve.txt"
    log_GS_lengths,GS_lengths,_ = read_file(file_name) #read files
    d = len(log_GS_lengths[0])
    

    # try:
    #     shutil.rmtree(dir, ignore_errors=False, onerror=None)
    # except FileNotFoundError:
    #     pass
    
    i = 0
    tours = 2
    log_GS_lengths_dic = {}
    log_GS_lengths_dic[(0,d,0,0)] = log_GS_lengths[0]
    for llb in range(d-50,0,-5):
        T_tmp = []
        RAM_tmp = []
        beta = d - llb
        f = dim4free_wrapper(theo_dim4free_fun1, beta)
        if beta - f < 131:
            for t in range(tours):
                if 2*tours*(i+1) >= len(log_GS_lengths):
                    break
                # show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],0,d,0,t)
                show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],log_GS_lengths[2*tours*i + 2*t + 1],llb,beta,f,t)

                log_GS_lengths_dic[(llb,beta,f,t)] = log_GS_lengths[2*tours*i + 2*t + 1]

                sim_gs = calculate_sim_log_gs_lengths_for_pump(log_GS_lengths[2*tours*i + 2*t],d,llb,d)
                # sim_gs = calculate_sim_log_gs_lengths_for_pump_test(log_GS_lengths[2*tours*i + 2*t],d,llb,d)
    
                show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t+1],sim_gs,llb,beta,f,2)
    
            i += 1

        

    return log_GS_lengths_dic




def gen_pump_figures_more_rounds(dir):

    file_name = "../../PnjBKZcost_test_result_test6_pumpData_g6k_down_sieve.txt"
    log_GS_lengths,GS_lengths,_ = read_file(file_name) #read files
    d = len(log_GS_lengths[0])
    

    try:
        shutil.rmtree(dir, ignore_errors=False, onerror=None)
    except FileNotFoundError:
        pass
    
    i = 0
    tours = 2
    log_GS_lengths_dic = {}
    log_GS_lengths_dic[(0,d,0,0)] = log_GS_lengths[0]
    for llb in range(d-50,0,-5):
        T_tmp = []
        RAM_tmp = []
        beta = d - llb
        f = dim4free_wrapper(theo_dim4free_fun1, beta)
        if beta - f < 131:
            for t in range(tours):
                if (tours+1)*(i+1) >= len(log_GS_lengths):
                    break
                # show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],0,d,0,t)
                show_pump_gs_figure(dir,log_GS_lengths[(tours+1)*i],log_GS_lengths[(tours+1)*i + t +1],llb,beta,f,t+1)

                # log_GS_lengths_dic[(llb,beta,f,t+1)] = log_GS_lengths[tours*i + t]

                sim_gs = calculate_sim_log_gs_lengths_for_pump(log_GS_lengths[(tours+1)*i],d,llb,d)
    
                show_pump_gs_figure(dir,log_GS_lengths[(tours+1)*i+t+1],sim_gs,llb,beta,f,11+t)

                # print(11+t)
    
            i += 1

        

    # return log_GS_lengths_dic



#=========================================input==============================================

if __name__ == "__main__":
    dir = "pump_simulator_g6k/"
    log_GS_lengths_dic = gen_pump_figures_more_rounds(dir)
  
    

    # gen_pump_figures_more_rounds("pump_simulator_theo_more_rounds/")
