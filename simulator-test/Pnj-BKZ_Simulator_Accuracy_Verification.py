# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import os


from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma

#read GS length
def read_gs_lengths(filename):
    f = open(filename,'r')
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
    blocksizes = f.read()
    f.close()
    blocksizes = blocksizes.split(' ')
    blocksizes.remove('')
    blocksizes = [int(_) for _ in blocksizes]
    return blocksizes

def read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d):
    
    dir = "gs-lengths-simulator/n=%d,alpha=%s,jump=%d,Pumpdown=%s,Blocksize=%d,Tours=%d/n=%d,alpha=%s,jump=%d,Pumpdown=%s,Blocksize=%d,Tours=%d,test_number=%dd=%d,/" %(n,alpha_,jump,Pumpdown,Blocksize,Tours, n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d)
    
    

    log_GS_lengths,GS_lengths = read_gs_lengths(dir+"rr_set.txt")
    
    return log_GS_lengths,GS_lengths,dir

#When jump=1, pnj-bkz simulator degenerates to CN11 simulator
def pnjBKZ_simulator(log_rr0,beta,N,d,jump):
    l = []
    for i in range(d):
        temp = log_rr0[i]
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
    print(N)
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
    return l_


def CN11_simulator(log_rr0,beta,N,d):
    l = []
    for i in range(d):
        temp = log_rr0[i]
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
    print(N)
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
    return l_

def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(15, 10), dpi=100)

    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[log_gs_length[_] for _ in range(t,dimension) ],marker="o",c='none',edgecolors='b')
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="x",c='r')
    #plt.title("square error = %f" %(square_error),x=0.9,y=0.875)
    plt.rcParams.update({'font.size':20})
    plt.legend(['Experimental','jump simulator'])
   

    plt.savefig(dir+"n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.show()
    plt.close()
    
def show_gs_slope_figure2(dir,log_gs_length,sim_log_gs_lengths,CN_11_mid_log_gs,n,dimension,alpha_,square_error1,square_error2,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(15, 10), dpi=100)

    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[log_gs_length[_] for _ in range(t,dimension) ],marker="^",c='none',edgecolors='k')
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="x",c='r')
    plt.scatter([_+1 for _ in range(t,len(CN_11_sim)) ],[CN_11_sim[_] for _ in range(t,len(CN_11_sim)) ],marker="o",c='none',edgecolors='b')
    plt.xlabel("Index")
    plt.ylabel("log GS norm")
    plt.legend(['Experimental','pnj-BKZ simulator','BKZ2.0 simulator'])
  

    plt.savefig(dir+"Comparing figure, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    #plt.savefig(dir+"Comparing figure, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(dimension,Blocksize,Jump,N+1))
    plt.show()
    plt.close()
    
def show_gs_slope_figure3(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error1,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(15, 10), dpi=100)

    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[log_gs_length[_] for _ in range(t,dimension) ],marker="^",c='none',edgecolors='k')
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="x",c='b')
    plt.xlabel("Index", fontsize=24)
    plt.ylabel("log GS norm", fontsize=24)
    plt.legend(['Experimental','pnj-BKZ simulator'], fontsize=24)
    #plt.title("n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d, prediction error = %f" %(n,alpha_,dimension,Blocksize,Jump,N+1,square_error1))
    plt.title("Dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d,SimError = %f" %(dimension,Blocksize,Jump,N+1,square_error1), fontsize=20)
    plt.tick_params(labelsize=22)

    plt.savefig(dir+"Comparing figure, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.show()
    plt.close()
    
def show_gs_slope_figure4(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error1,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(8, 10), dpi=100)
    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[sim_log_gs_lengths[_]/log_gs_length[_] for _ in range(t,dimension) ],s=12,marker="*",c='k',edgecolors='k')
    plt.ylim(0.9,1.1)
    
    #plt.axhline(1,color='k')
    plt.plot([1,dimension],[1,1],c='k')
    plt.plot([1,dimension],[1.05,1.05],c='k',linestyle='--')
    plt.plot([1,dimension],[0.95,0.95],c='k',linestyle='--')
    #plt.fill_between([1,dimension],[1.05,1.05],[0.95,0.95],facecolor='gray',edgecolor='k',alpha=0.3)
    
    plt.legend(['Real ||bi*||/ Sim ||bi*||'])
    #plt.legend(['Ratio','Ratio=1','Ratio within [0.95,1.05]'])
    plt.xlabel("Index")
    plt.ylabel("Ratio")
    #plt.legend(['Experimental','pnj-BKZ simulator'])
    #plt.title("n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.title("Dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d" %(dimension,Blocksize,Jump,N+1), fontsize=14)

    plt.savefig(dir+"Ratio figure, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.show()
    plt.close()
       
        
def ratio_txt(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error1,Blocksize,Jump,N,test_number):
    f = open(dir+"/ratio, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.txt"%(n,alpha_,dimension,Blocksize,Jump,N+1), "w")
    for _ in range(0,dimension):
        f.write(str(sim_log_gs_lengths[_]/log_gs_length[_])+' ')
    f.close()
    
def Ratio_figure(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error1,Blocksize,Jump,N,test_number):
    f = open(dir+"/ratio, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.txt"%(n,alpha_,dimension,Blocksize,Jump,1), "r")
    ratio1_ = f.read()
    ra = ratio1_.split(' ')
    ratio1 = []
    for _ in ra:
        #print(_)
        if _ !='':
            ratio1.append(float(_))
    #ratio1_ = [float(_) for _ in ra]
        
    #print(ratio1)
    Current_Tours1 = 5
    f0 = open(dir+"/ratio, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.txt"%(n,alpha_,dimension,Blocksize,Jump,Current_Tours1), "r")
    ratio4_ = f0.read()
    ra4 = ratio4_.split(' ')
    ratio4 = []
    for _ in ra4:
        #print(_)
        if _ !='':
            ratio4.append(float(_))
    
    Current_Tours2 = 10
    f1 = open(dir+"/ratio, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.txt"%(n,alpha_,dimension,Blocksize,Jump,Current_Tours2), "r")
    ratio8_ = f1.read()
    ra8 = ratio8_.split(' ')
    ratio8 = []
    for _1 in ra8:
        #print(_)
        if _1 !='':
            ratio8.append(float(_1))
    
    x = [_ for _ in range(len(ratio1))]
    plt.figure(figsize=(16, 9), dpi=100)
    plt.xlabel('Index', fontsize=24)
    plt.ylabel('Ratio', fontsize=24)
    plt.scatter(x, ratio1,s=24,marker="*",c='k',edgecolors='k', label=r"$\mathrm{Real}(  \ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }) / \mathrm{Sim}  (\ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }), \mathrm{Tours}=1$")
    plt.scatter(x, ratio4,s=16,marker="^",c='none',edgecolors='b', label=r"$\mathrm{Real}(  \ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }) / \mathrm{Sim}  (\ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }), \mathrm{Tours}=%d$"%Current_Tours1)
    plt.scatter(x, ratio8,s=16,marker="o",c='r',edgecolors='r', label=r"$\mathrm{Real}(  \ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }) / \mathrm{Sim}  (\ln_{}{\left \| \left \| \mathbf{b}_{i}^* \right \|  \right \| }), \mathrm{Tours}=%d$"%Current_Tours2)
    #plt.plot(x, ratio1, 'ok', label='Real ||bi*||/ Sim ||bi*|| Tours=1')
    #plt.plot(x, ratio8, '^b', label='Real ||bi*||/ Sim ||bi*|| Tours=8')
    plt.ylim(0.9,1.1)
    dimension=len(ratio1)
    plt.plot([1,dimension],[1,1],c='k')
    plt.plot([1,dimension],[1.05,1.05],c='k',linestyle='--')
    plt.plot([1,dimension],[0.95,0.95],c='k',linestyle='--')
    plt.plot([1,dimension],[1.025,1.025],c='k',linestyle=':')
    plt.plot([1,dimension],[0.975,0.975],c='k',linestyle=':')
    plt.tick_params(labelsize=20)

    plt.legend(fontsize=20, loc = 'upper right')
    plt.title("Dimension = %d, Blocksize = %d, Jump = %d" %(dimension,Blocksize,Jump), fontsize=20)
    plt.savefig(dir+"Ratio with different tours, Dimension = %d, Blocksize = %d, Jump = %d.png"%(dimension,Blocksize,Jump))
    
    
 
    
    
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

    return square_error
             
        
#=========================================input==============================================

if __name__ == "__main__":
    n, d, alpha_, jump, Pumpdown, Blocksize, Tours, test_number  = 75,252,"005",5,"True",95,12,1
    
    # n, d, alpha_, jump, Pumpdown, Blocksize, Tours, test_number  = 75,252,"005",9,"True",95,12,1

    #Blocksizes=[95,95,95,95,95,95,95,95]
    #Sim_Blocksizes=[95,95,95,95,95,95,95,95]
    
    #Blocksizes=[95,95,95,95,95,95,95,95,95,95,95,95]
    #Sim_Blocksizes=[95,95,95,95,95,95,95,95,95,95,95,95]
    
    Blocksizes=[95,95,95,95,95,95,95,95,95,95,95,95]
    Sim_Blocksizes=[95,95,95,95,95,95,95,95,95,95,95,95]
    
    #Blocksizes=[85,85,85,85,85,85,85,85,85,85,85,85]
    #Sim_Blocksizes=[85,85,85,85,85,85,85,85,85,85,85,85]
    
    d4f_value = "Theory 1 d4f"
    for i,k in enumerate(Sim_Blocksizes):
        if k>73 and k<89:
            Sim_Blocksizes[i] = k - 4
        elif k>89:
            Sim_Blocksizes[i] = k - 7
    print("d4f_value = %s" % d4f_value)

    #log_GS_lengths,GS_lengths,dir = read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d) #read
    
    #Total_test_number represents the number of pnj-BKZ reduction experiments
    Total_test_number = 20
    GS_sum = []
    for i in range(Total_test_number):
        log_GS_lengths,GS_lengths,dir = read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,i+1,d)
        if i==0:
            GS_sum = copy.deepcopy(log_GS_lengths)
        else:
            for k in range(Tours):
                for j in range(d):
                    GS_sum[k][j] = GS_sum[k][j] + log_GS_lengths[k][j]
    
    average = GS_sum.copy()
    
    #Calculate the average GS values of multiple experiments
    for i in range(Tours):
        for j in range(d):
            average[i][j] = average[i][j]/Total_test_number
        
    log_GS_lengths = average
    
    
    #log_GS_lengths we get are 2*log||bi*||, so we do following pre-processingg:
    for i in range(len(Blocksizes)+1):
        for j in range(d):
            log_GS_lengths[i][j] = log_GS_lengths[i][j]/2
    
    mid_log_gs=[]
    for N in range(len(Blocksizes)):
        if N == 0:
            if Blocksizes[0]>44:
                sim_log_gs_lengths = copy.deepcopy(pnjBKZ_simulator(log_GS_lengths[0], Sim_Blocksizes[0], 1, d, jump))
                CN_11_sim = copy.deepcopy(CN11_simulator(log_GS_lengths[0], Blocksizes[0], 1, d))
            else:
                sim_log_gs_lengths = copy.deepcopy(CN11_simulator(log_GS_lengths[0], Sim_Blocksizes[0], 1, d))
                CN_11_sim = copy.deepcopy(CN11_simulator(log_GS_lengths[0], Blocksizes[0], 1, d))
            mid_log_gs = copy.deepcopy(sim_log_gs_lengths)
            CN_11_mid_log_gs = copy.deepcopy(CN_11_sim)

            square_error1 = compute_square_error(sim_log_gs_lengths,log_GS_lengths[N+1],1)
            square_error2 = compute_square_error(CN_11_sim,log_GS_lengths[N+1],1)
            print("\nOur square error: %f, CN11 square error: %f" %(square_error1,square_error2))
           
        else:
            if Blocksizes[N]>44:
                sim_log_gs_lengths = copy.deepcopy(pnjBKZ_simulator(mid_log_gs, Sim_Blocksizes[N], 1, d, jump))
                CN_11_sim = copy.deepcopy(CN11_simulator(CN_11_mid_log_gs, Blocksizes[N], 1, d))
            else:
                sim_log_gs_lengths = copy.deepcopy(CN11_simulator(mid_log_gs, Sim_Blocksizes[N], 1, d))
                CN_11_sim = copy.deepcopy(CN11_simulator(CN_11_mid_log_gs, Blocksizes[N], 1, d))
            mid_log_gs = copy.deepcopy(sim_log_gs_lengths)
            CN_11_mid_log_gs = copy.deepcopy(CN_11_sim)
            
            square_error1 = compute_square_error(sim_log_gs_lengths,log_GS_lengths[N+1],1)
            square_error2 = compute_square_error(CN_11_sim,log_GS_lengths[N+1],1)
            print("\nOur square error: %f, CN11 square error: %f" %(square_error1,square_error2))
            # show_gs_slope_figure(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
            #show_gs_slope_figure2(dir,log_GS_lengths[N+1],sim_log_gs_lengths,CN_11_mid_log_gs,n,d,alpha_,square_error1,square_error2,Blocksizes[N],jump,N,test_number)
            
           
        show_gs_slope_figure3(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
        
         

        show_gs_slope_figure4(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
        
        
        ratio_txt(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
        
        # if (N+1)%4 == 0:
        #     show_gs_slope_figure5(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
        


    Ratio_figure(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number) #多个不同Tours对应的Ratio图画在同一张图上，然后用颜色区分不同的Tours所对应的Ratio
