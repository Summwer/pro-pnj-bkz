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
#from .cost import gaussian_heuristic, get_pump_time, dim4free_wrapper, default_dim4free_fun, theo_dim4free_fun1, theo_dim4free_fun2
from scipy.special import chdtr 
from .cost import *


#pump dimension estimation in default G6K
def pump_estimation1(log_rr,q,alpha, goal_margin=1.5, ebeta = 50):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    #goal_margin = 1 #1.5
    d = len(log_rr)
    # file.write(d)
    #compute n_expected
    m = d - 1
    stddev = q*float(alpha)
    target_norm = goal_margin * (stddev**2) * m + 1 

    for n_expected in range(2, d-2):
        # file.write(n_expected,d)
        x = (target_norm/goal_margin) * n_expected/(1.*d)
        if 4./3 * gaussian_heuristic(log_rr[d-n_expected:]) > x:
            break
    
    llb = d - ebeta
    # print(d,llb,ebeta)
    while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
        llb -= 1
        if llb < 0:
            break
    llb = max(0, llb)
    f = max(d-llb-n_expected, 0)

    #pump time
    pump_time = get_pump_time(n_expected,d)  # n_expected = beta -f , beta = d-llb
 
    return pump_time,llb,d-llb,f


def sieve_estimation(log_rr,q, alpha, succ_prob = 0.999):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(log_rr)
    p = 0.
    rp = 1.
    for beta in range(50,d):
        GH = gaussian_heuristic(log_rr[d-beta:])
        
        length=(GH/(sigma**2))
        psvp = chdtr(beta, length) #integral_chi_squared distribution (n,x)    
        if(psvp > succ_prob):
            #print(beta)
            break

    llb = d - beta
    # while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
    #     llb -= 1
    #     if llb < 0:
    #         break
    llb = max(0, llb)
    # f = max(d-llb-beta, 0)
    f = dim4free_wrapper(default_dim4free_fun, beta)

    #pump time
    pump_time = get_pump_time(beta - f,d)

    return pump_time,llb,d-llb,f



#pump dimension estimation in [PV21] with chi-squared distribution， one probability
def pro_sieve_estimation(log_rr,q, alpha, succ_prob = 0.999):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(log_rr)
    p = 0.
    rp = 1.
    for beta in range(50,d):
        GH = gaussian_heuristic(log_rr[d-beta:])
        
        length=(GH/(sigma**2))
        psvp = chdtr(beta, length) #integral_chi_squared distribution (n,x)
        # print(p)
        p += rp * psvp
        rp = 1. - p
        if p > succ_prob:
        # if(psvp > succ_prob):
            #print(beta)
            break

    llb = d - beta
    # while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
    #     llb -= 1
    #     if llb < 0:
    #         break
    llb = max(0, llb)
    # f = max(d-llb-beta, 0)
    f = dim4free_wrapper(default_dim4free_fun, beta)

    #pump time
    pump_time = get_pump_time(beta - f,d)

    return pump_time,llb,d-llb,f


#pump estimation cum previous failure probability
def pro_sieve_estimation_20230414(log_rr,q, alpha, succ_prob = 0.999):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(log_rr)
    p = 0.
    rp = 1.
    pre_psvp = 0.
    cum_pump_time  = 0.
    for beta in range(50,d):
        GH = gaussian_heuristic(log_rr[d-beta:])
        length=(GH/(sigma**2))
        psvp = chdtr(beta, length) #integral_chi_squared distribution (n,x)
        
       
        # pump_time = get_pump_time(beta - f,d)
        # f = dim4free_wrapper(default_dim4free_fun, beta)
        f = dim4free_wrapper(theo_dim4free_fun2, beta)
        cum_pump_time += get_pump_time(beta-f, d) *(psvp-pre_psvp)
        if psvp > succ_prob:
            break
        pre_psvp = psvp

    llb = d - beta
    # while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
    #     llb -= 1
    #     if llb < 0:
    #         break
    llb = max(0, llb)
    # f = max(d-llb-beta, 0)
    # f = dim4free_wrapper(theo_dim4free_fun1, beta)
    f = dim4free_wrapper(theo_dim4free_fun2, beta)
    # f = dim4free_wrapper(default_dim4free_fun, beta)

    #pump time
    

    return cum_pump_time,llb,d-llb,f,GH



#pump estimation cum previous failure probability
def pro_sieve_estimation_20230609(log_rr,q, alpha, succ_prob = 0.999):
# def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
    """
    Return min pump time cost estimate according to progressive sieve following [Duc18]

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param beta: current bkz blocksize
    :param target_norm: target_norm = sigma^2*d following the LWE distribution

    """
    alpha = float(alpha)
    sigma = alpha * q
    d=len(log_rr)
    pre_psvp2 = 0.
    cum_pump_time  = 0.
    beta = 50
    flag1 = False
    flag2 = False
    while(beta <= d):
        GH = gaussian_heuristic(log_rr[d-beta:])
        length=(GH/(sigma**2))
        psvp1 = chdtr(beta, length)
        psvp2 = chdtr(beta, 4/3.*length)
        cum_pump_time += get_pump_time(beta,d) *(psvp2-pre_psvp2)
        if(psvp1 > succ_prob and not flag1):
            dsvp1 = beta
            flag1 = True
        if(psvp2 > succ_prob and not flag2):
            dsvp2 = beta
            flag2 = True
        if(flag1 and flag2):
            break
        pre_psvp2 = psvp2
        
        beta +=1

    llb = d - dsvp1
    f = dsvp1 - dsvp2
    llb = max(0, llb)

    return cum_pump_time,llb,d-llb,f,GH




#test whether the gs-lengths after pump satisfies the HKZ-reduced result.
# def HKZ_satification(rr_kappa, beta, q, alpha, succ_prob = 0.999):
# # def pump_estimation2(log_rr,q, alpha, succ_prob = 0.99, ebeta = 50, goal_margin=1.5):
#     """
#     Return min pump time cost estimate according to progressive sieve following [Duc18]

#     :param log_rr: vector of log(squared) Gram-Schmidt norms
#     :param beta: current bkz blocksize
#     :param target_norm: target_norm = sigma^2*d following the LWE distribution

#     """
#     sigma = q * alpha
#     length=(rr_kappa/(sigma**2))
#     psvp = chdtr(beta, length) #integral_chi_squared distribution (n,x)
#     print(rr_kappa,psvp)
#     if psvp > succ_prob:
#         return True
#     return False



def pump_estimation(log_rr,q,alpha,succ_prob=0.999):
    pump_cost1,llb1,blocksize1,f1 = pump_estimation1(log_rr,q,alpha)
    pump_cost2,llb2,blocksize2,f2 = pump_estimation2(log_rr,q,alpha,succ_prob=succ_prob)
    print(llb1,blocksize1,f1)
    print(llb2,blocksize2,f2)
    if pump_cost1 > pump_cost2:
        return pump_cost1,llb1,blocksize1,f1
    else:
        return pump_cost2,llb2,blocksize2,f2


# #pump dimension estimation in [PV21] with chi-squared distribution， cumulative probability
# def pump_estimation3(log_rr,q,alpha,target_norm, ebeta = 50, goal_margin=1.5):
#     """
#     Return min pump time cost estimate according to progressive sieve following [Duc18]

#     :param log_rr: vector of log(squared) Gram-Schmidt norms
#     :param beta: current bkz blocksize
#     :param target_norm: target_norm = sigma^2*d following the LWE distribution

#     """
#     alpha = float(alpha)

#     d=len(log_rr)
#     pt = 0
#     for beta in range(50,d):
#         sigma = alpha * q
#         GH = gaussian_heuristic(log_rr[d-beta:])

#         length=(GH/(sigma**2))
#         pn = chdtr(beta, length) #integral_chi_squared distribution (n,x)
#         p = (1 - pt) * pn
#         pt = pt+p
#         if pt > 0.9:
#             #print(beta)
#             break

#     llb = d - beta
#     while gaussian_heuristic(log_rr[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
#         llb -= 1
#         if llb < 0:
#             break
#     llb = max(0, llb)
#     f = max(d-llb-beta, 0)

#     #pump time
#     pump_time = get_pump_time(beta) 

#     return pump_time,llb,d-llb,f




# log_rr0 = [17.274924047614363, 17.242786842904202, 17.177568643430558, 17.21229237442227, 17.197899511104488, 17.183679994466193, 17.15644841751796, 17.099785039914813, 17.00155849799918, 16.953674513722028, 16.873452879577606, 16.885441590484827, 16.780150847372614, 16.771669348744936, 16.711204668323166, 16.67228188070577, 16.656582858074298, 16.660120045469117, 16.589710829324382, 16.469044056597806, 16.46903343164758, 16.390759571331415, 16.384405578268034, 16.414486000203233, 16.324312959856755, 16.244669566762404, 16.214044068910447, 16.17780966321049, 16.114686808449573, 16.119296978170752, 16.198263397902537, 16.027082683810168, 16.04977433434799, 15.922750291260286, 15.920138864871289, 15.768913584562355, 15.93736727569393, 15.78081581656566, 15.74966124813674, 15.724120484851758, 15.647637912799077, 15.510359764447088, 15.529221807178343, 15.483654803879563, 15.447915206717555, 15.418457267614711, 15.369748598257477, 15.363159297848627, 15.294101247359018, 15.152331276991072, 15.191633630882004, 15.201290234069504, 15.128455148961525, 15.147429586014825, 15.067182736472757, 15.013038172720444, 14.988744140070523, 14.953910965733858, 14.894145896468762, 14.862714520689769, 14.839523159616483, 14.710152273705486, 14.557505268344231, 14.564126395690575, 14.657907083053688, 14.629266207317979, 14.520935204853053, 14.475967771101857, 14.456435918588758, 14.414901680957623, 14.452269693277898, 14.339317167497235, 14.336758006445175, 14.23017689748226, 14.240042224967109, 14.188798423928635, 14.14008063018434, 14.080495405990323, 14.032786488745426, 14.01843395496149, 13.971927445400338, 13.951852272206436, 13.828335463222128, 13.791755561431408, 13.791917703517765, 13.740363749561277, 13.681596093208176, 13.66150822996532, 13.604533782388241, 13.566979198479684, 13.475086431907913, 13.481690874007505, 13.457156439554016, 13.42021659847496, 13.371754106108813, 13.333080643471455, 13.280281490664251, 13.281392260419153, 13.232417481756862, 13.141965122419627, 13.129481178631865, 13.137606022011633, 13.0435439168787, 12.89765839459612, 12.852959243024214, 12.787101105553928, 12.851630624625024, 12.645332448895399, 12.89147806981824, 12.748209993079207, 12.792182412033142, 12.746529335894394, 12.611095523860106, 12.534651538851968, 12.529008477591129, 12.45401272982848, 12.375081366471672, 12.422561922926196, 12.341445087649435, 12.351486341143739, 12.269881714430825, 12.268325010469102, 12.269219479983857, 12.14889419980636, 12.126317785778916, 12.025268325867474, 12.038514665678399, 11.988991005085294, 11.926329974420465, 11.848455543022629, 11.846203713905718, 11.843927075561519, 11.802118092709495, 11.757858270590159, 11.619336456970874, 11.6488275508746, 11.545246836673675, 11.580886492997283, 11.508074205920279, 11.429815724024516, 11.464972990588198, 11.445668159594405, 11.37628122669323, 11.220246176967084, 11.250101941353572, 11.21076461573791, 11.199835061430855, 11.177617716900233, 11.03747870737116, 11.016143184485404, 10.97187134715494, 10.934831678813051, 10.817796842115504, 10.885599454664913, 10.74140734906076, 10.740277625733048, 10.778064974299385, 10.70504509372303, 10.5337900880216, 10.48060336515011, 10.454811084209378, 10.603012781149612, 10.390664364234645, 10.349470173217322, 10.315379200592428, 10.27484957244597, 10.303335689998187, 10.265249702466447, 10.220858246560356, 10.147881980353016, 10.168498187713533, 10.081033022708002, 10.116309848382713, 9.9753533175682, 9.950599243860873, 9.964092081613128, 9.889674125242578, 9.859537170198996, 9.833382920927367, 9.674613535903244, 9.766423024045755, 9.710741737276209, 9.715184538761514, 9.682421514373559, 9.64926045340324, 9.539899547055965, 9.492725974669792, 9.460093426657586, 9.534040942844888, 9.413778480362087, 9.290640352447111, 9.358984850434382, 9.339505244183362, 9.28875655197216, 9.249047985590382, 9.227706973040355, 9.18528752963321, 9.13946092398042, 9.068024904839676, 9.035255184255817, 8.966017215863786, 8.96533230155463, 8.920913748402228, 8.920932052285899, 8.884673932013898, 8.874181971346006, 8.796834409903575, 8.740392567957807, 8.640074279060409, 8.584320028969078, 8.618675295272997, 8.603971174673632, 8.521989667822558, 8.641082592732051, 8.408402499840893, 8.372594286886322, 8.181014076357927, 8.31202526884671, 8.13849633476479, 8.196535694784346, 8.200399746529255, 8.07077970986425, 7.976528781213202, 8.113837043685162, 7.879261732661709, 7.77510092787272, 7.77816561167845, 7.680665284640871, 7.748034236909384, 7.661871981297552, 7.686596073643418, 7.552335317328355, 7.4939150036988105, 7.4078346288207975, 7.266234890937727, 7.151405749733018, 7.053947501552089, 6.983679815253555, 6.995559882750528, 7.02568800627537, 6.945275014590552, 7.019128778968396, 6.912269985450448, 6.724225296013704, 6.612922487255253, 6.497488542059103, 6.477305737940289, 6.605175354170545, 6.3640502452912795, 6.25085579747251, 6.183425403812565]

# n = 75
# alpha = 0.005
# q = 5639
# m = 250
# goal_margin = 1.5
# target_norm = goal_margin * (alpha*q)**2 * m + 1 #298110.259375



# # # beta = 40
# # pump_time1,llb1,beta1,f1 = pump_estimation1(log_rr0,q,alpha,goal_margin)
# pump_time2,llb2,beta2,f2 = pump_estimation2(log_rr0,q,alpha,goal_margin)
# print("-----")
# print(llb1,beta1,f1) #102, 149, 23
# print(llb2,beta2,f2)
# print(beta1 - f1 )
# print(2**pump_time1)
# print("...")
# print(beta2 - f2 )
# print(2**pump_time2)



# def d_svp_prediction(l, cost_model, progressive_sieve):
#     """
#     Dimension of sieve/progressive sieve chosen to find target vector.
#     Computes the probabilistic cumulated cost value for given gs-lengths.
#     :l: log(||b_i^*||), i = 0,...,d-1
#     :cumulated_proba: current scuccess cumulated probability of gs-lengths in reduction 
#     :cost_model: 1: gate model
#                  2: sec model with threads=32, gpus = 2 
#     :progressieve_sieve: True: progressieve sieve
#                          False: normal sieve

#     return value:
#         dsvp/avgdsvp: average dimension value to sieve
#         dsvp_r: the largest dimension value to sieve
#         G_sieve: the cumulated time cost for sieve
#         B_dsvp: the maximal memory cost for sieve

#     """

#     d = len(l)
#     l_ = [2*_*log(2.) for _ in l]
#     cumulated_proba = 0.
#     remaining_proba = 1. 
#     if(cumulated_proba>= 1.):
#         return (0,0,0,0)
#     if not progressive_sieve:
#         G_sieve, B_sieve = float("inf"), float("inf")
#         for dsvp in range(50, d):
#             psvp = 1.
#             #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
#             gh = gaussian_heuristic(l_[d-dsvp:])
#             psvp *= chisquared_table[dsvp].cum_distribution_function(gh)
            
#             p = cumulated_proba + remaining_proba * psvp
#             rp = 1 - p

#             if rp < 0.001:
#                 dsvp += dsvp * rp #Avoid too small of dsvp
#                 break
#         G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
        
#     else:           
#         #predict dimension of last sieve: progressive sieve
#         p = cumulated_proba
#         rp = 1. - p
#         avgdsvp = 0.
#         avgG_sieve,avgB_sieve = 0.,0.

#         for dsvp in range(50, d):
            
#             psvp = 1.
#             #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
#             gh = gaussian_heuristic(l_[d-dsvp:])
#             psvp *= chisquared_table[dsvp].cum_distribution_function(gh)
#             G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
        
#             avgdsvp += dsvp * rp * psvp
#             avgG_sieve = log(2**avgG_sieve+(2**G_sieve) * rp * psvp)
#             avgB_sieve = max(B_sieve,avgB_sieve)

#             p += rp * psvp
#             rp = 1. - p
#             #print(dsvp, gh)
        
#             if rp < 0.001:
#                 #raise ""
#                 #print(rp,avgdsvp,dsvp * rp)
#                 avgdsvp += dsvp * rp #Avoid too small of dsvp
                
#                 avgG_sieve = log(2**avgG_sieve + ((2**G_sieve) * rp))
#                 return (avgG_sieve,avgB_sieve,avgdsvp,dsvp)
                
#     return (G_sieve,B_sieve,dsvp,dsvp)
