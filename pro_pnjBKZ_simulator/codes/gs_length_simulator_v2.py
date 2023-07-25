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

from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma
from .util import  get_current_slope


'''v2: Add the jump value into the simulator'''

#read GS length 
def read_gs_lengths(filename):
    f = open(filename,'r')
    #f = open("rr_set_bkz_50.txt",'r')
    GS_lengths = f.read()
    GS_lengths = GS_lengths.split('\n')
    GS_lengths.remove('')
    log_GS_lengths = []
    for i in range(len(GS_lengths)):
        gs = GS_lengths[i]
        gs = gs.split(' ')
        gs.remove('')
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

def read_file(n,alpha_,jump,Pumpdown,d,beta,tours):
    #n=50,alpha=0.01,k=1,f=0,beta=10-66
    dir = "gs-lengths-simulator/"+"n=%d,alpha=%s,jump=%d,Pumpdown=%s,Blocksize=%d,Tours=%d,d=%d,/" %(n,alpha_,jump,Pumpdown,beta,tours,d)

    log_GS_lengths,GS_lengths = read_gs_lengths(dir+"rr_set.txt")
    
    return log_GS_lengths,GS_lengths,dir

def calculate_sim_log_gs_lengths_test(log_rr0,beta,N,d,jump):
    if N == 0 or beta == 0:
        return log_rr0
    l = []
    for i in range(d):
        temp = log_rr0[i]/2
        l.append(temp)
    #print(len(l))
    l_ = []
    J = jump
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
            beta_ = min(beta, d - k)
            if k % J == 0:
                f = min(k + beta , d)
                #之前的 f = min(k + beta, d)
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
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                f = min(k - (k % J) + beta, d)
                sumf = 0
                sumk = 0
                for i in range(f):
                    sumf = sumf + l[i]
                for i in range(k):
                    sumk = sumk + l_[i]
                logV = sumf - sumk
                if flag == True:
                    if logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1] < l[k]:
                        l_k = logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1]
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

                    l_k = logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
        for k in range(d - beta,d - 45):
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
                l_k = logV / beta_ + cd[beta_-1]
                if tours == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
                
            #print(l_)
        k_ = min (len(randomSquaredAverages), beta)
        #last 45 norms
        sumf=0;
        sumk=0;
        for i in range(d):
            sumf = sumf + l[i]
        for i in range(d - 45):
            sumk = sumk + l_[i]
        logV = sumf - sumk
        for k in range(d - k_, d):
            l_k = logV/k_ + rk[k + 45 - d]
            if tours == 0:
                l_.append(l_k)
            else:
                l_[k] = l_k
                #Set l[i] as l_[i] 
        for i in range(d):
            l[i] = l_[i]
    return [_*2 for _ in l_]


#beta<45
def calculate_sim_log_gs_lengths_test1(log_rr0,beta,N,d):
    if N == 0 or beta == 0:
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
    k_ = min (len(randomSquaredAverages), beta)
    cd = []
    for i in range(beta):
        CD = lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2
        cd.append(CD)
    #print(cd)
    #print(N)
    
    for j in range(N):
        flag = True
        for k in range(d):
            beta_ = min(beta, d - k)
            f = min(k + beta, d)
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
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
                    flag = False
                else:
                    l_k = l[k]
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                l_k = logV / beta_ + cd[beta_-1]
                if j == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
            #print(l_)
        #Set l[i] as l_[i]         
        for i in range(d):
            l[i] = l_[i]
    return [_*2 for _ in l_]

def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error,beta,N):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(4,17) #set range of y_ticks
    # plt.xlim(-5,210)
    t = 0
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    plt.title("n=%d, alpha = %s, dimension = %d, Blocksize = %d, Current Tours= %d,square error = %f" %(n,alpha_,dimension,beta,N,square_error))
    
    try:
        os.mkdir(dir+"gs-lengths-gh simulator/")
    except FileExistsError:
        pass
    #plt.savefig(dir+"n=%d,alpha=%s,d=%d.png" %(n,alpha_,dimension))
    plt.show()

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
             
        
#=========================================input==============================================

if __name__ == "__main__":
    pnj = 'pnj-bkz'
    n, d, alpha_, jump, Pumpdown, beta, tours = 75,252,"005",6,"True",70,20
    
    log_GS_lengths,GS_lengths,dir = read_file(n,alpha_,jump,Pumpdown,d,beta,tours) #read files
    
    # n = 60 alpha = 0.005
    #n, d, alpha_, jump, Pumpdown, beta, tours = 60,199,"005",6,"True",74,20
    #dir =""
    #log_rr = [16.38126336180708, 16.38126336180708, 16.38126336180708, 16.300063311805157, 16.28526080671212, 16.200879750843132, 16.15017550609299, 16.10473373239724, 15.968768175423325, 15.881426496560731, 16.015954257952952, 15.919246988953013, 15.901707041273523, 15.897310931136548, 15.795467006051574, 15.697834303536137, 15.67380491215057, 15.710831477860223, 15.653236282434943, 15.540878979853252, 15.439630336855728, 15.248749637622526, 15.373931586873098, 15.310691130773993, 15.310245320407665, 15.260739755646755, 15.239873079375265, 15.136690430270765, 15.112866195935238, 14.960413876531286, 14.988798360157938, 14.887601953670924, 14.835314021634499, 14.740827301642604, 14.772745527438804, 14.747484161214011, 14.60072371551572, 14.584518732477365, 14.631460411158667, 14.453781433735925, 14.405698948828103, 14.299887950434854, 14.297701728487871, 14.269258888449704, 14.201803258789703, 14.216181916387727, 14.142455503742164, 14.054753339250167, 14.014977864850877, 13.937392461905613, 13.945388391863256, 13.89868538828191, 13.64607450256329, 13.596591354922982, 13.477312053812806, 13.600131070817847, 13.66530379835224, 13.607505879395763, 13.449806601720288, 13.368478539656865, 13.414859661668794, 13.430751617243944, 13.268338136790458, 13.17074021343514, 13.222920646886905, 13.171178346566697, 13.08604101276061, 13.006647623201616, 12.924014695034833, 12.916250954218858, 12.896041620344317, 12.816725320484995, 12.78597383459337, 12.739807095663899, 12.610353307753739, 12.676334409200425, 12.522333357255262, 12.521093214582766, 12.465671400458946, 12.43283219599929, 12.364339040880038, 12.125443447076156, 12.301927660698386, 12.21533092225634, 12.155975358030101, 12.028654108330292, 12.035603925256515, 12.013867424567065, 11.903645644670366, 11.877939752537738, 11.82139730614894, 11.772196674757932, 11.688136616775845, 11.677384537802421, 11.612938597326933, 11.525938600738865, 11.513169032131058, 11.47510404519387, 11.433425933994934, 11.384010353774725, 11.248266760092198, 11.258129701316111, 11.195494016061206, 11.042522746499731, 11.00899242661159, 11.140358662124928, 11.009565352275786, 10.834394844959434, 10.798412534414124, 10.592103595838813, 10.864436651909639, 10.793119249756627, 10.608658806321781, 10.578695059922936, 10.523042440219383, 10.524027024145184, 10.45024795621493, 10.38384853987705, 10.334772258490984, 10.419193135470763, 10.252955052937754, 10.210188550218604, 10.184420774630135, 10.151727852313975, 10.102693952932901, 9.965960530460874, 9.910912238942737, 9.918441661072206, 9.810536841637772, 9.710036225598149, 9.751077010714273, 9.669016610499906, 9.662444075240122, 9.576760038364517, 9.52424726928894, 9.607089495710003, 9.503476014033453, 9.345034307902516, 9.330021478139194, 9.289309372576522, 9.127023211823092, 9.08636070517182, 9.150220947761113, 9.100061164925867, 8.995966186149214, 8.922174204612228, 8.868473290852979, 8.871145437682658, 8.72359878381206, 8.834131101504463, 8.753121563351332, 8.7366663697295, 8.587158253424828, 8.583562987980487, 8.56248150005601, 8.54497650985077, 8.483112456402486, 8.393118930043975, 8.395148404782423, 8.33086754068521, 8.29582030141405, 8.169050640249454, 8.107128457516257, 8.0954990327283, 8.056440926303823, 8.014822446981105, 7.99047219946272, 7.881332607646005, 7.898052328107874, 7.646447204021777, 7.814131953892247, 7.699039206476208, 7.641009249429498, 7.550913001236983, 7.628978491254942, 7.486340431907946, 7.44645250793637, 7.353983226169353, 7.379739076961872, 7.386251822289362, 7.172224347107424, 7.175189404160399, 7.10675906143906, 6.974970508201688, 6.9290080313663145, 6.828733669048952, 6.759574349620633, 6.760428756655113, 6.719518114593199, 6.668640216873729, 6.593835410289552, 6.561739557540937, 6.321799874696559, 6.481096626028254, 6.3381905537090955, 6.376563048663048, 6.129768192148478, 6.317944524997385, 6.143269783462458]
      
    
    #N = k
    N = 1
    #log_GS_lengths = [ log_rr for _ in range(N+1)]
    
    if beta>44:
        sim_log_gs_lengths = calculate_sim_log_gs_lengths_test(log_GS_lengths[0], beta, N, d, jump)
    else:
        sim_log_gs_lengths = calculate_sim_log_gs_lengths_test1(log_GS_lengths[0], beta, N, d)
    for i in range(N+1):
        for j in range(d):
            log_GS_lengths[i][j] = log_GS_lengths[i][j]
    square_error = compute_square_error(sim_log_gs_lengths,log_GS_lengths[N],1)
    show_gs_slope_figure(dir,log_GS_lengths[N],sim_log_gs_lengths,n,d,alpha_,square_error,beta,N)

    print(get_current_slope(sim_log_gs_lengths,0,d))
