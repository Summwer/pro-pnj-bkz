#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 14:07:45 2022

@author: summer
"""
import time
import os
from .util import get_current_slope, predict_slope_with_beta, get_current_slope
from .pnjbkz_simulator import simulatepnjbkzloop,PPNJBKZSimulate
from .pump_estimation import pump_estimation
from math import log
'''v1: previous strategy loop out from slope '''

#optimize the quality comparison from pump time to slope
def simtime_slope(log_rr0,target_norm,d,sbeta,ebeta,ubeta,jump): 

    eslope = predict_slope_with_beta(ebeta)

    loop = 1
    log_rr,pnj_time = simulatepnjbkzloop(log_rr0, ubeta, loop, d,jump)
    uslope = get_current_slope(log_rr,0,d)
    # file.write(uslope,eslope)
    if uslope >= eslope - 0.00001 or eslope >=0 :
        return pnj_time ,loop,log_rr
    
    slope_same_times = 0 #avoid the endless loop
    pnj_time1 = 0
    uslope1 = -1
    while( uslope1 < eslope - 0.00001 and abs(uslope1-uslope)>0.00001): #and slope_same_times < 5):
        # file.write(loop,slope_same_times,uslope,eslopes)
        uslope1 = get_current_slope(log_rr,0,d)
        log_rr1 = [_ for _ in log_rr]
        pnj_time1 +=pnj_time
        log_rr,pnj_time = simulatepnjbkzloop(log_rr1, ubeta,1, d, jump)
        uslope = get_current_slope(log_rr,0,d)
        #if abs(uslope1 - uslope) < 0.00001 :
        #    slope_same_times+=1
        loop+=1
        
    if uslope < eslope - 0.00001 :
        return -1,-1,log_rr0
    
    #loop = loop - slope_same_times
    #file.write(uslope1,eslope)
    return pnj_time1 ,loop-1,log_rr1
   


#optimize the quality comparison from pump time to slope
# def simtime_slope(log_rr0,target_norm,d,sbeta,ebeta,ubeta,jump): 

#     eslope = predict_slope_with_beta(ebeta)

#     loop = 1
#     log_rr,pnj_time = simulatepnjbkzloop(log_rr0, ubeta, loop, d,jump)
#     uslope = get_current_slope(log_rr,0,d)
#     # file.write(uslope,eslope)
#     if uslope >= eslope - 0.00001 or eslope >=0 :
#         return pnj_time ,loop,log_rr
    
#     slope_same_times = 0 #avoid the endless loop
#     while( uslope < eslope - 0.00001 and slope_same_times < 5):
#         # file.write(loop,slope_same_times,uslope,eslopes)
#         uslope1 = get_current_slope(log_rr,0,d)
#         log_rr,pnj_time = simulatepnjbkzloop(log_rr0, ubeta, loop, d, jump)
#         uslope = get_current_slope(log_rr,0,d)
#         if abs(uslope1 - uslope) < 0.00001 :
#             slope_same_times+=1
#         loop+=1
        
#     if slope_same_times == 5:
#         return -1,-1,log_rr0
    
#     loop = loop - slope_same_times

#     return pnj_time ,loop-1,log_rr
    
# simtimes:
# input: log_rr0,target_norm, beta , f
# return: T,log_rr,loop
'''
def simtime(log_rr0,target_norm,d,sbeta,ebeta,ubeta,goal_margin,q,succ_prob): 

    log_rr = [_ for _ in log_rr0]

    #initiate the gs-lengths after sufficient rounds of sbeta
    pump_cost,_,_,_= pump_estimation(log_rr,q,alpha,target_norm,goal_margin,succ_prob=succ_prob) #logfec(rr_sbeta)
    
    #target gs-lengths quality
    log_rr_ebeta,_,_ = PPNJBKZSimulate(log_rr0,ebeta,d) 
    eslope = get_current_slope(log_rr_ebeta,0,d) 
    pump_cost_ebeta,_,_,_ = pump_estimation(log_rr,q,alpha,target_norm,goal_margin,succ_prob=succ_prob) #logfec(rr_ebeta)

    loop = 1
    log_rr,pnj_time = simulatepnjbkzloop(log_rr0, ubeta, loop, d)
    uslope = get_current_slope(log_rr,0,d)
    # file.write(eslope,uslope)
    # pump_cost,_,_,_ = pump_estimation(log_rr,ubeta,target_norm)
    
    if pump_cost <= pump_cost_ebeta :
        return pnj_time ,loop,log_rr
    
    
    slope_same_times = 0
    while(pump_cost > pump_cost_ebeta and slope_same_times < 5 ):
        #pnj-BKZ*loop time 
        uslope1 = get_current_slope(log_rr,0,d)
        log_rr,pnj_time = simulatepnjbkzloop(log_rr0, ubeta, loop, d)
        pump_cost,_,_,_ = pump_estimation(log_rr,q,alpha,succ_prob=succ_prob)
        #logfec(rr_sbeta)
        uslope = get_current_slope(log_rr,0,d)
        if abs(uslope1 - uslope) < 0.0001 :
            slope_same_times+=1
        loop+=1
        
    if slope_same_times == 5 :
        return -1,-1,log_rr0
    
    loop = loop - slope_same_times

    return pnj_time ,loop-1,log_rr
'''

#return the ubeta that makes sbeta --> ebeta and the minimum cost.
def simtimeopt(log_rr0,target_norm,sbeta,ebeta,d, jump):
    mincost = -1
    minbeta = d
    minloop = -1
    minlog_rr = [_  for _ in log_rr0]
    
    for ubeta in range(ebeta+1,d):
        t_bkz_ubeta,loop,log_rr = simtime_slope(log_rr0,target_norm,d,sbeta,ebeta,ubeta, jump)
        
        if(loop!=-1 and mincost < 0 or t_bkz_ubeta < mincost):
            mincost = t_bkz_ubeta 
            minbeta = ubeta
            minloop = loop
            minlog_rr = [_ for _ in log_rr]
        # help to break the loop in advance.    
        if (ubeta >= minbeta+3): 
            break
        
    return mincost,minbeta,minloop,minlog_rr

        
#Optimize blocksizes strategy: BS={ebeta: (log_rr,sbeta,ebeta,ubeta,t_bkz,loop,total_cost)}
# total_cost = T(BKZ(sbeta))+t_bkz*loop
def gen_bkzstrategy(file,log_rr0,target_norm, BS, sbeta, ebeta, d, jump):
    
    # predict a beta_goal such that basis can be solved through this beta_goal, estimate 2016.
    
    for beta in range(sbeta+1,ebeta+1):
        file.write( "Simulating optimal strategy: " +str( sbeta ) + "->" + str( beta))
        file.write("\n")
        mincost = -1 #minimum total cost
        for ssbeta in range(sbeta,beta):
            if(ssbeta == sbeta): 
                localcost = 0
                log_rr = [_ for _ in log_rr0]
            if(ssbeta > sbeta):
                if(ssbeta not in BS): 
                    log_rr = [_ for _ in log_rr0]
                    BS= gen_bkzstrategy(log_rr, target_norm, BS, sbeta, ssbeta, d, jump )
                localcost = BS[ssbeta][-1]  
                log_rr = [_ for _ in BS[ssbeta][0]]
            t_bkz,ubeta,loop,sim_log_rr = simtimeopt(log_rr,target_norm,ssbeta,beta,d, jump)
            if (mincost < 0 or mincost > localcost + t_bkz) and loop!=-1:
                mincost = round(localcost + t_bkz,4)
                log_rr = [_ for _ in sim_log_rr]
                BS[beta] = (log_rr,get_current_slope(log_rr,0,len(log_rr)),ssbeta,beta,ubeta,round(t_bkz,4),loop,mincost)
            
            if loop == -1:
                file.write( "\rSimulating optimal strategy: " +str(sbeta)+ "->" +str( beta) + " finished")
                file.write("\n")
                #sys.stdout.flush()
                
                file.write("We cannot find an replaced blocksizes stragety for beta_goal larger than %d." %beta)
                file.write('\n')
                # for key in BS:
                    # file.write(str(BS[key][1:]))
                
                
                return BS
  
        file.write( "Simulating optimal strategy: " +str(sbeta)+ "->" +str( beta) + " finished")
        file.write("\n")
        #sys.stdout.flush()
   
    #for key in BS:
    #    file.write(BS[key][1:])
    
    
    return BS


#Optimize blocksizes strategy: BS={ebeta: (log_rr,sbeta,ebeta,ubeta,t_bkz,loop,total_cost)}
# total_cost = T(BKZ(sbeta))+t_bkz*loop
def BS_display(file,process,sbeta,ebeta,target_norm,goal_margin,alpha,q):
    Result = [ebeta]
    Ubeta = []
    bkz_totalcost = process[ebeta][-1]
    log_rr = process[ebeta][0]
    slope = process[ebeta][1]
    file.write(str(log_rr))
    file.write('\n')
    file.write("Predict slope = %f" %slope)
    file.write('\n')
    
    t_pump,llb,blocksize,df= pump_estimation(log_rr,q,alpha)
    total_cost = bkz_totalcost + t_pump
    while(ebeta in process):
        Result.append(process[ebeta][2])
        Ubeta.append((process[ebeta][4],process[ebeta][-2],process[ebeta][1],process[ebeta][5]))
        ebeta = process[ebeta][2]
    Result.reverse()
    Ubeta.reverse()
    if len(process)==0: return
    file.write(str(sbeta))
    blocksizes = []
    for i in range(len(Ubeta)):
        file.write( "->( %d * %d ( %f ), cost = %.4f s )-> %d"  %( Ubeta[i][0], Ubeta[i][1] , Ubeta[i][2],Ubeta[i][3], Result[i+1] ))
        for _ in range(Ubeta[i][1]):
            # if j > 5:
            #     break
            blocksizes.append(Ubeta[i][0])
    file.write("+pump(%d,%d,%d)" %(llb,blocksize,df)) 
    file.write('\n')
    str_blocksizes = '[' + ",".join([str(_) for _ in blocksizes]) +']'
    file.write(str_blocksizes)
    file.write('\n')
    file.write("Time cost for pnj-BKZ = %.2f h," %(bkz_totalcost/3600.))
    file.write('\n')
    file.write("Time cost for pump = %.2f h," %((total_cost - bkz_totalcost)/3600.))
    file.write('\n')
    file.write("Total time cost = %.2f h." %(total_cost/3600.))
    file.write('\n\n')
    file.close()

    return blocksizes,round(bkz_totalcost/3600.,2),round((total_cost - bkz_totalcost)/3600.,2),round(total_cost/3600.,2)



def gen_bkzstrategy_from_slope(file,log_rr0,target_norm,sbeta, target_beta, d):
    BS = {}
    
    if target_beta < 50:
        target_beta = 50
        
    target_slope = predict_slope_with_beta(target_beta)
    ebeta = max(target_beta,sbeta+1)
    file.write("Gen BKZ strategy from slope: theoretically predict slope(beta = %d) = %f" %(target_beta,target_slope))
                
    file.write("Find BKZ strategy: BKZ-" +str(sbeta) + " -> " + str(ebeta))
    t_gen = time.time()
    BS = gen_bkzstrategy(file,log_rr0,target_norm, BS, sbeta, ebeta, d)
    t_gen = time.time() - t_gen
    file.write("Time cost for generate blocksizes strategy: %f s" %t_gen)
    
    return BS,ebeta

def gen_bkzstrategy_min_cost(file,log_rr0,target_norm,sbeta,target_beta,d,jump,goal_margin,alpha,q,succ_prob):
    target_beta = max(target_beta,sbeta+1)
    BS = {}
    file.write("Find BKZ strategy: BKZ-" +str(sbeta) + " -> " + str(target_beta+25))
    file.write("\n")
    t_gen = time.time()
    
    BS = gen_bkzstrategy(file,log_rr0,target_norm, BS, sbeta,target_beta+25, d, jump)
    min_total_cost = float('inf')
    minbeta = sbeta
    target_slope = predict_slope_with_beta(target_beta) 
    for beta in range(min(target_beta,max(BS.keys())),min(target_beta+25,1+max(BS.keys()))):
        if beta in BS:
            bkz_total_cost = BS[beta][-1]
            log_rr = BS[beta][0]
            slope = get_current_slope(log_rr,0,d)
            t_pump,_,n_max, f= pump_estimation(log_rr,q,alpha,succ_prob=succ_prob)
            
            total_cost = bkz_total_cost + t_pump
            #file.write(total_cost,beta,slope,n_max)
            if total_cost < min_total_cost and n_max - f <= 200: #and abs(log_rr[0]-2*log(q))>0.1: #and  slope >= target_slope + (jump-1)*0.0003: eliminate vector q
                min_total_cost = total_cost
                minbeta = beta
    if min_total_cost == float('inf'):
        min_total_cost = total_cost
        minbeta = beta
    t_gen = time.time() - t_gen
    file.write("Time cost for generate blocksizes strategy: %f s\n" %t_gen)
    # file.write(min_total_cost, minbeta)
    return BS,minbeta,round(t_gen,2)



def gen_strategy_v1(file_path, n,alpha,log_rr0,target_norm,ebeta,d,jump,goal_margin,q,succ_prob):
    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)
    file = open(file_path+"/%d-%s-%d-%s.log" %(n,alpha[2:],jump,str(goal_margin*10)[:2]),'w')
    sbeta = 50
    file.write("Previous Strategy:\nn = %d, alpha = %s, b= %d, jump = %d, goal_margin = %.2f" %(n,alpha,ebeta,jump,goal_margin))

    file.write(str(log_rr0))
    file.write("\n")


    BS,minbeta,t_gen = gen_bkzstrategy_min_cost(file,log_rr0,target_norm,sbeta,ebeta, d, jump,goal_margin,alpha,q,succ_prob)
    blocksize_strategy, bkz_cost, pump_cost, total_cost = BS_display(file,BS,sbeta,minbeta,target_norm,goal_margin,alpha,q)
    
    return blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen

if __name__ == "__main__":
    
    #================================input==============================    
    #45, 020
    log_rr0 = [15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.228624292904, 15.182292534663446, 15.162373907298708, 15.138644751080975, 15.091995860973437, 14.96644254134134, 14.99680878236674, 14.916819500777095, 14.96554158380541, 14.929606200402734, 14.791515980863515, 14.750253657682144, 14.776438340453105, 14.685656956271723, 14.718991042040972, 14.628189578750614, 14.59529951093489, 14.376841008327638, 14.511010916017426, 14.523025404346434, 14.339179306578496, 14.349979237643218, 14.241839876747152, 14.196443565128352, 14.159957858581512, 14.15804275546448, 14.163075145354732, 14.058424367513343, 14.053950708723946, 14.015821448440162, 13.91212708183317, 13.908484725268655, 13.878417950173048, 13.800577609136264, 13.75878709605201, 13.74097365650203, 13.646854335781137, 13.652513805672886, 13.59686826908838, 13.538019507049594, 13.5083487099844, 13.422689567155986, 13.334854636269116, 13.465552372122108, 13.341482566147604, 13.243076746545487, 13.22214131050418, 13.191343524926452, 13.168316677732939, 13.00929390949677, 13.099724603241617, 13.051646827911208, 12.897246766430602, 12.937920658520257, 12.779079362592576, 12.754994039231123, 12.695760645463976, 12.779168393603555, 12.617030977564966, 12.751141710742438, 12.719186028829865, 12.632300648084737, 12.55445896433962, 12.451292588618525, 12.270931132147973, 12.432574974801677, 12.27199834733439, 12.332967950895801, 12.321950341955509, 12.232345402337048, 12.179036537643373, 12.163857649143878, 12.149465507189754, 11.991161407827473, 11.967442531594406, 11.982863662036072, 11.921394560244403, 11.847225628142283, 11.82865161032098, 11.788666629695605, 11.701690034840384, 11.711082072710221, 11.672688031947862, 11.557849322221147, 11.54944385913907, 11.467720169165093, 11.455318537736549, 11.355473894220106, 11.355029670178638, 11.316434032184496, 11.213516652314254, 11.267020030848478, 11.18504140480683, 11.035877878438376, 11.090036194226961, 11.02473258691957, 11.00423586406651, 10.93995961025409, 10.848323847319522, 10.874372816498825, 10.766039411595553, 10.706427152660826, 10.615313780862362, 10.644846204392383, 10.627185995324842, 10.598262210034052, 10.447158573354729, 10.461434882716404, 10.47227460253594, 10.441014909576266, 10.412982827923905, 10.385008044071409, 10.32524350072967, 10.184745329319535, 10.28682165120034, 10.225183057920658, 10.226140542933374, 10.129783542173232, 10.027205843065964, 10.057080296889746, 9.99424529276246, 9.98888825128024, 9.946943243314038, 9.894390673791495, 9.882186034207576, 9.854365305403741, 9.791219294543982, 9.754425727343701, 9.729454814847141, 9.69776865703124, 9.618344892588295, 9.577861954458522, 9.471779410397373, 9.529955781495124, 9.459646828604143, 9.40625345802174, 9.373490207825713, 9.260204604216636, 9.191898369894243, 9.254069550546768, 9.16657512877421, 9.11818134101556, 9.07617663545665, 9.080106490390477, 9.009551092865417, 8.953189287792851, 8.903192579347484, 8.863116366494006, 8.804155464906295, 8.745830194606834, 8.633542654835042, 8.641057198894435, 8.574905982086054, 8.48906020727908, 8.405924097405261, 8.360685355512762, 8.284673763465582, 8.169381534400637, 8.197018240710102, 8.247927830882139, 8.085118017889771, 8.051851088488325, 8.024306987218742, 8.031869783609238, 7.867171437106096, 7.857597851862546, 7.840169615485406, 7.699371099196372, 7.667285916378668, 7.481696438930203, 7.510337754843281, 7.429626878030258, 7.3337850263202915, 7.402704187369093, 7.1463900212133264, 7.277037746981621, 7.125007253385143, 7.013473461297095, 7.0772909550272045]

    target_norm = 453604.6816
    jump = 3
    sbeta = 50
    ebeta = 89
    goal_margin = 1.5
    
    d = len(log_rr0)
    
    
    MAX_TIME = float("inf")
    n = 45
    alpha = "0.020"
    
    print("Previous Strategy:\nn = %d, alpha = %s, b= %d, jump = %d, goal_margin = %.2f" %(n,alpha,ebeta,jump,goal_margin))
    file_path = "test"
    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)
    #blocksizes = gen_bkzstrategy_from_slope2(log_rr0,target_norm,sbeta, ebeta, d)
    blocksize_strategy, bkz_cost, pump_cost, total_cost,t_gen = gen_strategy_v1(file_path+"/previous_strategy",n,alpha,log_rr0,target_norm,ebeta,d,jump,goal_margin)