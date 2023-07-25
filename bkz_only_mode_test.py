#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 09:53:07 2022

@author: summer

#Implement strategy v1/v2 for all lwe challenges for different jump(1,2,3) and goal_margin(1-1.5)
#store them in strategy_test folder,
#Sort the strategy and time cost out of strategy_test.xls
"""
from multiprocessing import Pool
import requests

import sys
import xlwt
import os
from pro_pnjBKZ_simulator.codes.usvp_simulator_v2 import gen_bkzstrategy_min_cost_v2
from pro_pnjBKZ_simulator.codes.usvp_simulator_v1 import gen_strategy_v1
from pro_pnjBKZ_simulator.codes.usvp_simulator_default import usvp_simulator_default
from pro_pnjBKZ_simulator.codes.usvp_simulator_bkz_only_mode import usvp_simulator_bkz_only_mode
from gen_original_gs import gen_original_gs

def gen_lwes():
    ns = [40,45,50,55,60,65,70,75,80,85,90]
    alphas = ["0.005","0.010","0.015","0.020","0.025","0.030","0.035"]
    
    unsolved_lwes = []
    for i in range(7):
        if(i<5):
            unsolved_lwes.append((ns[i],alphas[-1-i]))
        if(i==5):
            unsolved_lwes.append((75,"0.010"))
        if(i==6):
            unsolved_lwes.append((90,"0.005"))
    
    solved_lwes = []
    for i in range(6):
        if(i==0):
            for j in range(ns.index(90)):
                solved_lwes.append((ns[j],alphas[0]))
        if(i==1):
            for j in range(ns.index(75)):
                solved_lwes.append((ns[j],alphas[1]))
        if(i>=2):
            k = i-2
            for j in range(ns.index(60)-k):
                solved_lwes.append((ns[j],alphas[i]))
    return sorted(solved_lwes+unsolved_lwes)






    
    
def read_preprocess_data(file_name):
    f = open(file_name,'r')
    data = f.readlines()
    try:
        if data[-4][0] == '[':
            preprocess_cost = float(data[-5].split(' ')[-2])
            ebeta = int(data[2].split('-')[1:][0].split(" ")[0])
            svp_dim = int(data[2].split('-')[1:][1].split("\n")[0])
            target_norm = float(data[-3])
            log_rr0 = data[-4]
            log_rr0 = log_rr0.replace("]\n","")
            log_rr0 = log_rr0.replace("[","")
            log_rr0 =[ float(_) for _ in log_rr0.split(", ")]
        
            return ebeta,svp_dim,target_norm,log_rr0,round(preprocess_cost/3600.,2)
    except IndexError:
        pass
    
    return None


def read_q(filename,n,alpha):
    if not os.path.isfile(filename):
        url = ("https://www.latticechallenge.org/lwe_challenge/challenges/"
               "LWE_{n:d}_{alpha:03d}.txt")
        url = url.format(n=n, alpha=alpha)
        r = requests.get(url)
        m = "Cannot retrieve challenge; server response was: '%s'. \n URL was: %s" % (r.reason, url)
        if not r.status_code == 200:
            raise ValueError(m)
        fn = open(filename, "w")
        fn.write(r.text)
        fn.close()

    data = open(filename, "r").readlines()

    return int(data[2])

    

def gen_strategy_for_lwe(lwe):
    (n,alpha) = lwe

    data_dir = "pro_pnjBKZ_simulator/strategy_pre_process"
    goal_margin0 = 1
    #lwes = gen_lwes() #[(40, '0.030'),(40, '0.035'),(45, '0.015'),(45, '0.020'),(45, '0.025'),(50, '0.005'),(50, '0.010'),(50, '0.015'),(50, '0.020'),(55, '0.005'),(55, '0.010'),(55, '0.015'),(60, '0.005'),(60, '0.010'),(60, '0.015'),(65, '0.005'),(65, '0.010'),(70, '0.005'),(70, '0.010'),(75, '0.005'),(80, '0.005'),(85, '0.005'),(90, '0.005')]

    MAX_TIME = float("inf")
    threads = 20
    it = 0
    col_id = 2
    # preprocess_blocksizes = [51,52,53,54,55,60,62,65]
    preprocess_blocksizes = []

    
    (n,alpha) = lwe
    file_name = data_dir + "/%d-%s.log" %(n,alpha[2:])
    data = read_preprocess_data(file_name)
    # if data==None:
    #     work_book.save(file_path+'/all_strategies-%d-%s.xls' %(n,alpha[2:]))
    #     return
    # else:
    #     (ebeta,svp_dim,target_norm,log_rr0,preprocess_cost) = data
    #     d = len(log_rr0)
    
    file_path = "pro_pnjBKZ_simulator/default-strategy-test_pro_085"

    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)

    q = read_q("pro_pnjBKZ_simulator/lwechallenge/%d-%s.log" %(n,alpha[2:]),n,int(round(float(alpha) * 1000)))

    time_cost_all = []

    goal_margin = 1.5
    data = gen_original_gs(n,float(alpha), goal_margin)
    (ebeta,svp_dim,target_norm,log_rr0,preprocess_cost) = data
    d = len(log_rr0)

    
    #default strategy
    blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen = usvp_simulator_default(file_path+"/g6k-default",n,alpha,log_rr0,target_norm,ebeta,goal_margin,q)
    time_cost_all.append(total_cost)






if __name__ == '__main__':
    
    
    
    # lwes =[(80, '0.005')] #[(70,'0.010')]
    #lwes = [(70, '0.010'),(75, '0.005'),(80, '0.005'),(85, '0.005'),(90, '0.005')] #[(40, '0.030'),(40, '0.035'),(45, '0.015'),(45, '0.020'),(45, '0.025'),(50, '0.005'),(50, '0.010'),(50, '0.015'),(50, '0.020'),(55, '0.005'),(55, '0.010'),(55, '0.015'),(60, '0.005'),(60, '0.010'),(60, '0.015'),(65, '0.005'),(65, '0.010'),(70, '0.005'),],
    
    #gen_strategy_for_lwe(lwes[5])
    lwes = [(80, '0.005'),(90,'0.005'),(50,'0.025'),(55,'0.020'),(45,'0.030'),(60,'0.015'),(40,'0.035')]
    # lwes =[(80, '0.005')] #[(70,'0.010')]
    with Pool(1) as p:
       print(p.map(gen_strategy_for_lwe, lwes))

    # for lwe in lwes:
    #     gen_strategy_for_lwe(lwe)
