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


"""basic functions used in pro-pnj-bkz"""




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

def read_file(n,alpha_,k,fs,blocksize_start,blocksize_end):
    #n=50,alpha=0.01,k=1,f=0,beta=10-66
    dir = "pro-bkz-tests/gs-lengths-simulator/"+"n=%d,alpha=%s,k=%d," %(n,alpha_,k)
    if len(fs)==1:
        dir += "f="+str(fs[0])+","
    else:
        dir += "f="+str(fs[0])+"-"+str(fs[-1])+","
    dir += "beta="+str(blocksize_start)+"-"+str(blocksize_end)
    dir +='/'

    log_GS_lengths,GS_lengths = read_gs_lengths(dir+"rr_set.txt")
    blocksizes = read_blocksizes(dir+"blocksizes.txt" )
    times = read_times(dir+"times.txt")
    
    return log_GS_lengths,GS_lengths,blocksizes,times,dir



#calculate gh
def gaussian_heuristic(log_rr):
    """
    Return squared norm of shortest vector as predicted by the Gaussian heuristic.

    :param log_rr: vector of log(squared) Gram-Schmidt norms

    """
    n = len(list(log_rr))
    log_vol = sum([x for x in log_rr])
    log_gh =  1./n * (log_vol - 2 * ball_log_vol(n))
    return exp(log_gh)


def ball_log_vol(n):
    """
    Return volume of `n`-dimensional unit ball

    :param n: dimension

    """
    return (n/2.) * log(pi) - lgamma(n/2. + 1)

#calculate slope in G6K
def get_current_slope(log_rr, start_row=0, stop_row=-1):
    """
    Return current slope of log-gs-lengths

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param start_row: start row
    :param stopr_row: stop row

    """
    
    x = [ _ for _ in log_rr]
    n = stop_row - start_row
    i_mean = (n - 1) * 0.5 + start_row
    x_mean = sum(x)/n
    v1, v2 = 0.0, 0.0
    for i in range(start_row, stop_row):
        v1 += (i - i_mean) * (x[i] - x_mean)
        v2 += (i - i_mean) * (i - i_mean)
    return round(v1 / v2,6) 


def show_pump_gs_figure(dir,log_gs_length0,log_gs_length,llb,beta,f,t,gs0_name,gs_name):
    plt.figure(figsize=(15, 10), dpi=100)
    d  = len(log_gs_length)
    plt.scatter([_+1 for _ in range(d)],log_gs_length0,marker="*")
    plt.scatter([_+1 for _ in range(d)],log_gs_length,marker="*")
    plt.plot([llb,llb], [min(log_gs_length),max(log_gs_length)])
    plt.plot([llb+f,llb+f], [min(log_gs_length),max(log_gs_length)])
    function3 = gs0_name
    function4 = gs_name
    function1 = "llb=%d"%llb
    function2 = "llb+f=%d"%(llb+f)
    plt.legend([function1,function2,function3,function4])
    plt.title("gs-lengths")

    try:
        os.mkdir(dir)
    except FileExistsError:
        pass
    plt.savefig(dir+"pump(%d,%d,%d)-%d.png" %(llb,beta,f,t))
    plt.close()



#simulate gs and show it in figure
def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,blocksize,pnj,f,n,dimension,alpha_,blocksize_start, blocksize_end,square_error,t,N):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(4,17) #set range of y_ticks
    # plt.xlim(-5,210)
    t = 0
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    plt.title("gs-length change for %s , f = %d, n=%d, alpha = %s, dimension = %d, beta = %d, #Tours= %d,square error = %f" %(pnj,f,n,alpha_,dimension,blocksize,N,square_error))
    
    try:
        os.mkdir(dir+"gs-lengths-gh simulator/")
    except FileExistsError:
        pass
    plt.savefig(dir+"gs-lengths-gh simulator/f=%d,n=%d,alpha=%s,k=1,beta=%d,d=%d.png" %(f,n,alpha_,blocksize,dimension))
    plt.show()
    
    
def show_figure(log_gs_length,dimension):
    plt.figure(figsize=(15, 10), dpi=100)
    
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    #plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    #plt.title("gs-length change for %s , f = %d, n=%d, alpha = %s, dimension = %d, beta = %d, #Tours= %d,square error = %f" %(pnj,f,n,alpha_,dimension,blocksize,N,square_error))
    
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

    return square_error



def dim4free_wrapper(dim4free_fun, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize)
    return int(min((blocksize - 40)/2, dim4free))



def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)


def theo_dim4free_fun1(blocksize):
    """
    Theoretical Dimension-for-free function 1 without e in [Duc18]
    """

    return int(blocksize*log(4/3.)/log(blocksize/2./pi)) 


def theo_dim4free_fun2(blocksize):
    """
    Theoretical Dimension-for-free function 2 with e in [Duc18]
    """

    return int(blocksize*log(4/3.)/log(blocksize/2./pi/e)) 

#40, 0.035 , for pnj-BKZ, with pump/down = False
# def get_k1_k2_pnj_BKZ(beta):
#     if beta <= 40 and beta>=10:
#         k1 = 0.0589
#         k2 = 5.0574
#     elif beta<=75 and beta > 40:
#         k1 = 0.2032
#         k2 = - 2.564
#     else:
#         k1 = 0.2843
#         k2 = - 7.6978
#     return k1,k2



#g6k-cpu40, 0.035 , for pump, with pump/down = True, threads = 20, no gpus
# def get_k1_k2_pump(beta):
#     if beta <= 40 and beta>=0:
#         k1 = 0.0589
#         k2 = 5.0574
#     elif beta<=75 and beta > 40:
#         k1 = 0.0555
#         k2 = - 0.8001
#     elif beta > 75:
#         k1 = 0.2609
#         k2 = - 15.33
#     else:
#         k1 = 0
#         k2 = 0
#     # else: #30 > 80
#     #     k1 = 0.3642 
#     #     k2 = - 24.398 
#     return k1,k2



#g6k-gpu test svp 180, threads = 32, gpus = 2
# def get_k1_k2_pump(beta):
#     if beta <=0:
#         k1 = 0
#         k2 = 0
#     elif beta<= 50 and beta >0:
#         k1 = 0.035657
#         k2 = -2.317327
#     elif beta > 50 and beta < 120:
#         k1 = 0.063925
#         k2 = 0.292015
#     else: # >=91的部分需要重新测试一下，目前的结果较乐观
#         # k1 = 0.312
#         # k2 = - 17.04 
#         # k1 = 0.128224
#         # k2 = -5.518341
#         k1 = 1.02
#         k2 = -118.68
#     # else: #30 > 80
#     #     k1 = 0.3642 
#     #     k2 = - 24.398 
#     return k1,k2



# #g6k-cpu 90-0.005 , 
# def get_k1_k2_pump(beta):
#     if beta<=65 :
#         k1 = 0.1273
#         k2 = - 5.0914
#     elif beta > 65:
#         k1 = 0.312
#         k2 = - 17.04
#     else:
#         k1 = 0
#         k2 = 0
#     # else: #30 > 80
#     #     k1 = 0.3642 
#     #     k2 = - 24.398 
#     return k1,k2



# threads = 32, gpus = 2,  pnj-bkz cost
def get_k1_k2_pnj(beta,sieve):
    if beta >=0 and beta <10:
        k1 = 0 
        k2 = 0
    elif beta>=10 and beta<=42 and sieve == False:
        k1 = 0.03
        k2 = 5.188
    elif beta<=60 and sieve == False:
        k1 = 0.19
        k2 = -1.741
    elif beta <= 97:
        k1 = 0.056
        k2 = 7.85
    elif beta <= 118:
        k1 = 0.215
        k2 = -7.61
    elif beta <= 128:
        k1 = 0.314
        k2 = - 19.24
    else:
        k1 = 0.368
        k2 = -26.15
    return k1,k2


# threads = 32, gpus = 2, test pump
def get_k1_k2_pump(beta):
    if beta >=0 and beta <10:
        k1 = 0 
        k2 = 0
    elif beta>=10 and beta<=60:
        k1 = 0.035657
        k2 = -2.317327
    elif beta <= 96:
        k1 = 0.078794
        k2 = -0.039742
    elif beta <= 116:
        k1 = 0.231927
        k2 = -14.713430
    elif beta <= 128:
        k1 = 0.314
        k2 = -24.21
    else:
        k1 = 0.368
        k2 = -31.12
    # else: #30 > 80
    #     k1 = 0.3642 
    #     k2 = - 24.398 
    return k1,k2



#get pump time test in threads = 20
def get_pump_time(beta,d):
     #make sure not use the enum cost 
    k1, k2 = get_k1_k2_pump(beta) # threads = 20
    # k = (1/71.)*((1.33)**(beta/10.))
    T_pump = round((2 **(k1*(beta)+k2)),4)
    return T_pump  # n_expected = beta -f , beta = d-llb

    

#get pnj-BKZ time test in threads = 20
def get_pre_pnj_time(d,beta,f,jump):
    if beta <= 60:
        sieve = False
    else:
        sieve = True   
    k1,k2 = get_k1_k2_pnj(beta,sieve)
    c3, c4 = 0.018, -2.24
    T_pnj = 2**(k1*(beta-f)+k2)
    pre_pnj_time = T_pnj*(c3*d+c4)/jump

    return round(pre_pnj_time,4)
 



def delta_0f(k):
    """
    Auxiliary function giving root Hermite factors. Small values
    experimentally determined, otherwise from [Chen13]

    :param k: BKZ blocksize for which the root Hermite factor is required

    """
    small = (( 2, 1.02190),  # noqa
             ( 5, 1.01862),  # noqa
             (10, 1.01616),
             (15, 1.01485),
             (20, 1.01420),
             (25, 1.01342),
             (28, 1.01331),
             (40, 1.01295))

    k = float(k)
    if k <= 2:
        return (1.0219)
    elif k < 40:
        for i in range(1, len(small)):
            if small[i][0] > k:
                return (small[i-1][1])
    elif k == 40:
        return (small[-1][1])
    else:
        return (k/(2*pi*e) * (pi*k)**(1./k))**(1/(2*(k-1.)))


def predict_slope_with_beta(beta):
    # if beta < 50:
    #     raise TypeError("For slope prediction, beta must be larger than 50")
    delta = delta_0f(beta)
    slope = -4*log(delta)
    return slope



def dim4free_wrapper(dim4free_fun, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize)
    return int(min((blocksize - 40)/2, dim4free))





def get_beta_from_sieve_dim(sieve_dim,d,dim4free_fun):
    for beta in range(sieve_dim,d):
        f = dim4free_wrapper(dim4free_fun,beta)
        # print(beta,f,beta-f,sieve_dim)
        if beta - f >= sieve_dim:
            return beta

        