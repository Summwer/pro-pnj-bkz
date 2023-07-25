#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 14:07:45 2022

@author: summer
"""
import time
import os
from .util import get_current_slope, predict_slope_with_beta, get_current_slope, default_dim4free_fun, dim4free_wrapper
from .pnjbkz_simulator import simulatepnjbkzloop
from .pump_estimation import pump_estimation
from math import log

"""v2: generate optimal strategy from pro-bkz based on an enumeration algorithm
quality is slope"""

#MAX_TIME = float("inf")



#Determine whether sim_slope1 can add into BS
def BS_add(BS,sslope,sim_log_rr1,sim_slope,local_time,beta,loop,MAX_TIME):
    slope_list = sorted(BS.keys())
    local_time0 = BS[sslope][1]
    
    strategy0 = BS[sslope][2]
    slope_changes0 = BS[sslope][3]
    strategy = strategy0 +  [beta for _ in range(loop)]
    slope_changes = slope_changes0 + ['%d * %d (%f) , cost = %.4f s' %(beta, loop, sim_slope, local_time)]
    local_time = local_time0 + local_time
    # print("---------------------")
    
    if sim_slope not in BS:
        slope_list.append(sim_slope)
        slope_list = sorted(slope_list)
        index = slope_list.index(sim_slope)
        #Compare to the next slope
        if local_time < MAX_TIME:
            if index == len(slope_list) -1:
                BS[sim_slope] = (sim_log_rr1,round(local_time,4),strategy,slope_changes)
                # print("Condition1")
                # print(sim_slope,local_time)
            elif index < len(slope_list) -1:
                slope2 = slope_list[index+1]
                T2= BS[slope2][1]
                if local_time < T2:
                    BS[sim_slope] = (sim_log_rr1,round(local_time,4),strategy,slope_changes)
                    # print("Conditon2")
                    # print(index,slope_list[index-1],sim_slope,slope2,BS[slope_list[index-1]][1],local_time,T2)
            #Compare to the previous slope  
            if index > 0:      
                it = 1
                slope1 = slope_list[index-it]
                T1 = BS[slope1][1]
                while(local_time <= T1):
                    del BS[slope1]
                    it += 1
                    if index - it < 0:
                        break
                    slope1 = slope_list[index-it]
                    T1 = BS[slope1][1] 
                    #return BS#,False
                    # print("Condition3")
                    # print(slope1,sim_slope,T1,local_time)
    else:
        index = slope_list.index(sim_slope)
        slope2 = slope_list[index]
        T2= BS[slope2][1]
        if local_time < T2:
            BS[sim_slope] = (sim_log_rr1,round(local_time,4),strategy,slope_changes)
            # print("Condition4")
            # print(sim_slope,local_time,T2)
    # print("=====================")
    return BS#,True

#Compute the maximum tours a pnj-bkz-beta can do
#ps: If pnj-bkz-beta implements exceeds maximum tours, the basis quality won't improve.
def max_tour_for_pnjbkz_beta(BS,sslope,beta,jump,MAX_TIME):
    log_rr0 = BS[sslope][0]
    d = len(log_rr0)
    loop = 0
    local_time = 0
    sim_slope0 = sslope
    
    sim_log_rr1,pre_pnj_time1 = simulatepnjbkzloop(log_rr0, beta, loop+1, d, jump)
    sim_slope1 = get_current_slope(sim_log_rr1,0,d)
    while sim_slope1 - sim_slope0 > 0.00001:
        loop += 1
        local_time += pre_pnj_time1
        sim_log_rr0 = [_ for _ in sim_log_rr1]
        sim_slope0 = sim_slope1
        sim_log_rr1,pre_pnj_time1 = simulatepnjbkzloop(sim_log_rr0, beta, 1, d, jump)
        sim_slope1 = get_current_slope(sim_log_rr1,0,d)
        
        #Determine whether sim_slope1 can add into BS
        BS = BS_add(BS,sslope,sim_log_rr0,sim_slope0,local_time,beta,loop,MAX_TIME) 
        
        #BS,EXIST_sslope = BS_add(BS,sslope,sim_log_rr1,sim_slope1,local_time,beta,loop) 
        #if EXIST_sslope == False:
        #    return BS,False
        
    return BS#,True

def gen_bkzstrategy_v2(file,BS,sslope, d , jump, MAX_TIME):
    slope_list = sorted([_ for _ in BS if _ >= sslope])
    file.write("Blocksize strategy selection process: %f --> %f, len_BS = %d" %(slope_list[0],slope_list[-1],len(BS)))
    file.write("\n")
    # file.write("===================")
    #file.write(len(slope_list))
    # for slope in slope_list:
    #     bs = BS[slope]
    #     file.write(slope,bs[1],bs[2])
    
    # file.write("-------end---------")
    
    strategy = BS[sslope][2]
    if strategy == []:
        sbeta = 50
    else:
        sbeta = strategy[-1]
    #Test all BKZs' max tours for the initial rr0, and store the required slope
    #Initialize BS
    for beta in range(max(50,sbeta+1),170): #170-f = 146
        BS= max_tour_for_pnjbkz_beta(BS,sslope,beta,jump,MAX_TIME) 
        
        #BS,EXIST_sslope = max_tour_for_pnjbkz_beta(BS,sslope,beta,jump) 
        #if EXIST_sslope == False:
        #    break
    
    sorted_BS = sorted(BS.items(), key = lambda kv:(kv[1], kv[0]))
    sorted_BS.reverse()
    sorted_BS = [bs for bs in sorted_BS if bs[0]>sslope]
    while(sorted_BS!=[]):
        bs = sorted_BS[0]
        sslope = bs[0]
        slope_list = sorted([_ for _ in BS if _ >= sslope])
        file.write("Blocksize strategy selection process: %f --> %f, len_BS = %d" %(slope_list[0],slope_list[-1],len(BS)))
        file.write("\n")

        strategy = BS[sslope][2]
        if strategy == []:
            sbeta = 50
        else:
            sbeta = strategy[-1]
        #Test all BKZs' max tours for the initial rr0, and store the required slope
        #Initialize BS
        for beta in range(max(50,sbeta+1),170):
            BS= max_tour_for_pnjbkz_beta(BS,sslope,beta,jump,MAX_TIME) 
        
        sorted_BS = sorted(BS.items(), key = lambda kv:(kv[1], kv[0]))
        sorted_BS.reverse()
        sorted_BS = [bs for bs in sorted_BS if bs[0]>sslope]


    # if sorted_BS !=[]:
    #     bs = sorted_BS[0]
    #     #file.write(bs[0],bs[1][1],bs[1][2])
    #     sslope = bs[0]
    #     BS = gen_bkzstrategy_v2(file,BS,sslope,target_norm, d , jump,MAX_TIME)
    
    return BS

#A strategy for all of the situations about diffection redcution.
def gen_bkzstrategy_min_cost_v2(file_path,n,alpha,log_rr0, d, jump, goal_margin,MAX_TIME,q, ebeta,succ_prob):
    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)
    file = open(file_path+"/%d-%s-%d-%s.log" %(n,alpha[2:],jump,str(goal_margin*10)[:2]),'w')
    
    file.write("New Strategy:\nn = %d, alpha = %s, q = %d, jump = %d, goal_margin = %.2f, MAX_TIME = %.2f h" %(n,alpha,q,jump,goal_margin,MAX_TIME/3600.))
    file.write("\n")
    
    file.write(str(log_rr0))
    file.write("\n")
    
    T0 = time.time()
    #initialize BS
    sslope = get_current_slope(log_rr0,0,d)
    BS = {sslope: (log_rr0,0,[],[],1)}#store the intial slope (slope) = (sim_logrr0,T,blocksize_strategy)  
    
    BS = gen_bkzstrategy_v2(file,BS,sslope,d, jump,MAX_TIME)
    
    
    T = time.time() - T0
    BS_list = sorted([(_,BS[_][1]) for _ in BS if _ >= sslope])
    #file.write()
    file.write("len(BS) = %d" %len(BS_list))
    file.write("\n")
    # file.write("===================\n")
    # for bs in BS_list:
    #     file.write(str(bs[0])+","+str(bs[1])+str(BS[bs[0]][2])+","+str(BS[bs[0]][3]))
    #     file.write('\n')
    # file.write("-------end---------\n")
    
    file.write("Time cost for generating Strategy: %f s" %T)
    file.write("\n")
    
    min_cost = float("inf")
    target_slope = predict_slope_with_beta(ebeta) 
    Min_Strategy = []
    Min_rr0 = []
    
    for bs in BS:
        log_rr = BS[bs][0]
        local_cost = BS[bs][1]
        pump_cost,llb,n_max,f = pump_estimation(log_rr,q,alpha,succ_prob=succ_prob)

        total_cost = local_cost + pump_cost
        
        # file.write(str(BS[bs][2]))
        # file.write("+ pump{%d,%d,%d}, T_bkz = %.3f s, T_pump = %.3f s, T_total = %.3f" %(llb,n_max,f, local_cost, pump_cost, total_cost))
        # file.write("\n")
        if BS[bs][2] !=[]:
            max_beta = BS[bs][2][-1]
            max_beta_f = dim4free_wrapper(default_dim4free_fun, max_beta)
            pnjbkz_dim = max_beta-max_beta_f
        else:
            pnjbkz_dim = 0
        if min_cost > total_cost:# and n_max-f < 146 and pnjbkz_dim < 146: #and abs(log_rr[0]-2*log(q))>0.1: #eliminate vector q
        #and  slope >= target_slope + (jump-1)*0.0003:
            llb_min = llb
            n_max_min = n_max
            f_min = f
            Min_Strategy = BS[bs][2]
            Min_slope_changes = BS[bs][3]
            min_cost = total_cost
            min_bkz_cost = local_cost
            Min_rr0 = [_ for _ in BS[bs][0]]
    str_Min_rr0 = '[' + ",".join([str(_) for _ in Min_rr0]) +']'
    file.write(str_Min_rr0)
    file.write("\n")
    file.write("MinStrategy is: %s + pump(%d,%d,%d) , %s + pump(%d,%d,%d)" %(str(Min_Strategy),llb_min,n_max_min,f_min,str( Min_slope_changes),llb_min,n_max_min,f_min))
    file.write("\n")
    str_Min_Strategy = '[' + ",".join([str(_) for _ in Min_Strategy]) +']'
    file.write(str(str_Min_Strategy))
    file.write("\n")
    file.write("Time cost for pnj-BKZ = %.2f h," %(min_bkz_cost/3600.))
    file.write("\n")
    file.write("Time cost for pump = %.2f h," %((min_cost - min_bkz_cost)/3600.))
    file.write("\n")
    file.write("Total time cost = %.2f h." %(min_cost/3600.))
    file.write("\n")
    file.close()
    #file.write(BS)
    
    return Min_Strategy,round(min_bkz_cost/3600.,2),round((min_cost - min_bkz_cost)/3600.,2),round(min_cost/3600.,2),round(T,2)




# if __name__ == "__main__":
    
#     #================================input==============================    
   
#     #45, 020
#     log_rr0 = [15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.182292534663446, 15.162373907298708, 15.138644751080975, 15.091995860973437, 14.96644254134134, 14.99680878236674, 14.916819500777095, 14.96554158380541, 14.929606200402734, 14.791515980863515, 14.750253657682144, 14.776438340453105, 14.685656956271723, 14.718991042040972, 14.628189578750614, 14.59529951093489, 14.376841008327638, 14.511010916017426, 14.523025404346434, 14.339179306578496, 14.349979237643218, 14.241839876747152, 14.196443565128352, 14.159957858581512, 14.15804275546448, 14.163075145354732, 14.058424367513343, 14.053950708723946, 14.015821448440162, 13.91212708183317, 13.908484725268655, 13.878417950173048, 13.800577609136264, 13.75878709605201, 13.74097365650203, 13.646854335781137, 13.652513805672886, 13.59686826908838, 13.538019507049594, 13.5083487099844, 13.422689567155986, 13.334854636269116, 13.465552372122108, 13.341482566147604, 13.243076746545487, 13.22214131050418, 13.191343524926452, 13.168316677732939, 13.00929390949677, 13.099724603241617, 13.051646827911208, 12.897246766430602, 12.937920658520257, 12.779079362592576, 12.754994039231123, 12.695760645463976, 12.779168393603555, 12.617030977564966, 12.751141710742438, 12.719186028829865, 12.632300648084737, 12.55445896433962, 12.451292588618525, 12.270931132147973, 12.432574974801677, 12.27199834733439, 12.332967950895801, 12.321950341955509, 12.232345402337048, 12.179036537643373, 12.163857649143878, 12.149465507189754, 11.991161407827473, 11.967442531594406, 11.982863662036072, 11.921394560244403, 11.847225628142283, 11.82865161032098, 11.788666629695605, 11.701690034840384, 11.711082072710221, 11.672688031947862, 11.557849322221147, 11.54944385913907, 11.467720169165093, 11.455318537736549, 11.355473894220106, 11.355029670178638, 11.316434032184496, 11.213516652314254, 11.267020030848478, 11.18504140480683, 11.035877878438376, 11.090036194226961, 11.02473258691957, 11.00423586406651, 10.93995961025409, 10.848323847319522, 10.874372816498825, 10.766039411595553, 10.706427152660826, 10.615313780862362, 10.644846204392383, 10.627185995324842, 10.598262210034052, 10.447158573354729, 10.461434882716404, 10.47227460253594, 10.441014909576266, 10.412982827923905, 10.385008044071409, 10.32524350072967, 10.184745329319535, 10.28682165120034, 10.225183057920658, 10.226140542933374, 10.129783542173232, 10.027205843065964, 10.057080296889746, 9.99424529276246, 9.98888825128024, 9.946943243314038, 9.894390673791495, 9.882186034207576, 9.854365305403741, 9.791219294543982, 9.754425727343701, 9.729454814847141, 9.69776865703124, 9.618344892588295, 9.577861954458522, 9.471779410397373, 9.529955781495124, 9.459646828604143, 9.40625345802174, 9.373490207825713, 9.260204604216636, 9.191898369894243, 9.254069550546768, 9.16657512877421, 9.11818134101556, 9.07617663545665, 9.080106490390477, 9.009551092865417, 8.953189287792851, 8.903192579347484, 8.863116366494006, 8.804155464906295, 8.745830194606834, 8.633542654835042, 8.641057198894435, 8.574905982086054, 8.48906020727908, 8.405924097405261, 8.360685355512762, 8.284673763465582, 8.169381534400637, 8.197018240710102, 8.247927830882139, 8.085118017889771, 8.051851088488325, 8.024306987218742, 8.031869783609238, 7.867171437106096, 7.857597851862546, 7.840169615485406, 7.699371099196372, 7.667285916378668, 7.481696438930203, 7.510337754843281, 7.429626878030258, 7.3337850263202915, 7.402704187369093, 7.1463900212133264, 7.277037746981621, 7.125007253385143, 7.013473461297095, 7.0772909550272045]

#     target_norm = 453604.6816
#     jump = 3
#     sbeta = 50
#     ebeta = 89
#     goal_margin = 1.5
    
#     d = len(log_rr0)
    
    
#     MAX_TIME = float("inf")
#     n = 45
#     alpha = "0.020"
#     print("New Strategy:\nn = %d, alpha = %s, b= %d, jump = %d, goal_margin = %.2f, MAX_TIME = %.2f h" %(n,alpha,ebeta,jump,goal_margin,MAX_TIME/3600.))
#     file_path = "test"
#     if os.path.exists(file_path):
#         pass
#     else:
#         os.mkdir(file_path)
#     blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen = gen_bkzstrategy_min_cost_v2(file_path+"/new_strategy",n,alpha,log_rr0, d, jump, goal_margin,ebeta,MAX_TIME)
        
#     #BS,ebeta = gen_bkzstrategy_min_cost(log_rr0,target_norm,sbeta,ebeta, d)
#     #BS_display(BS,sbeta,ebeta,target_norm)
#     #blocksizes = gen_bkzstrategy_from_slope2(log_rr0,target_norm,sbeta, ebeta, d)
    
  