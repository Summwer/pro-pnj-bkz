#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 14:07:45 2022

@author: summer
"""
import time
import os
from .util import predict_slope_with_beta, get_current_slope
from .pnjbkz_simulator import simulatepnjbkzloop,PPNJBKZSimulate
from .pump_simulator import calculate_sim_log_gs_lengths_for_pump
from .util import gaussian_heuristic
from math import log
import matplotlib.pyplot as plt
from scipy.special import chdtr

'''Strategy in BKZ-only mode'''


def usvp_simulator_bkz_only_mode(file_path,n,alpha,log_rr0,target_norm,b,goal_margin,q,succ_prob):

    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)
    jump = 1
    file = open(file_path+"/%d-%s-%d-%s.log" %(n,alpha[2:],jump,str(goal_margin*10)[:2]),'w')   

    file.write(str(log_rr0))
    file.write("\n")

    
    local_cost = 0
    pump_cost = 0
    total_cost = 0

    d = len(log_rr0)

    log_rr = [_ for _ in log_rr0]
    
    # file.write(blocksizes)
    slopes = []
    bkz_time = []
    T0 = time.time()
    
    flag = False
    for beta in range(b,d):
        _,max_loop,_= PPNJBKZSimulate(log_rr0,beta,d,jump)
        # max_loop = 10
        for tours in range(max_loop):
            log_rr, pre_pnj_time = simulatepnjbkzloop(log_rr0, beta, tours, d,1)
            blocksize_strategy = [beta for _ in range(tours)]
            local_cost += pre_pnj_time
            total_cost = local_cost
            # print(log_rr[0],log(target_norm))
         
            # x = (target_norm/goal_margin) * beta/(1.*d)
            # if 4./3 * gaussian_heuristic(log_rr[d-beta:]) > x:
            #     flag = True
            #     break
            sigma = float(alpha) * q
            GH = gaussian_heuristic(log_rr[d-beta:])
            length=(GH/(sigma**2))
            p = chdtr(beta, length)
            if p >= succ_prob:
                flag = True
                break

        if flag == True:
            break

    T = time.time() - T0

    if(not flag):
        file.write("Compute Failed while using all blocksizes suggested!")
        file.write("\n")


    file.write("Time cost for generating Strategy: %f s" %T)
    file.write("\n")
    
    
    file.write(str(blocksize_strategy)+"\n")
    file.write("Time cost for pnj-BKZ = %.2f h," %(local_cost/3600.))
    file.write("\n")
    file.write("Total time cost = %.2f h." %(total_cost/3600.))
    file.write("\n")

    return blocksize_strategy, round(local_cost/3600.,2),round((total_cost - local_cost)/3600.,2),round(total_cost/3600.,2),round(T,2)



    
    
if __name__ == "__main__":
    
    #================================input==============================    
    # log_rr0 = [15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.650490582863549, 15.55933529703487, 15.61075772076923, 15.631086381240998, 15.528485580835149, 15.450593850852897, 15.218706772138988, 15.234187821037004, 15.072122289797365, 15.038837130002292, 14.975602769434236, 14.896543887541416, 14.718288778679362, 14.573677614962541, 14.581369961785299, 14.345071325866146, 14.421768054130096, 14.394105029408056, 14.260242997994858, 14.195327741225153, 14.176556853266774, 13.96799518075128, 13.951543928369318, 14.071471010282904, 13.84986373526558, 13.662433559293195, 13.436666392470515, 13.347087056170649, 13.169839091068972, 13.086767444197442, 13.181284540297348, 12.954145589923579, 12.980440380830165, 12.91072160780591, 12.83743674248719, 12.7817117606773, 12.708990386290974, 12.42801490862084, 12.603743794872914, 12.37531787449946, 12.293224605287456, 12.200994830636876, 12.015108698815048, 11.868418789842353, 12.001795838053667, 11.832364466676783, 12.098798184866487, 12.000219145781502, 11.73056521611099, 11.91043273872411, 11.694123819783027, 11.531791129239318, 11.557968691092144, 11.395052976617194, 11.348830359438624, 11.159010946366063, 11.102727272410892, 11.080042531052616, 10.893179990443617, 11.173654481447215, 11.146137382372885, 10.8808667743356, 10.738613333839876, 10.57617004270703, 10.427671140362362, 10.215022768867083, 9.986839899003382, 9.91345242741361, 9.719547755200528, 9.92117218169636, 9.751365981931396, 9.81922113939154, 9.621797857397633, 9.527697972807522, 9.379945405978164, 9.463532925168083, 9.34426679438853, 9.236597952425795, 9.113727245791885, 9.058851273571893, 8.828267601418364, 8.71881492254498, 8.606056776232887, 8.377059425847358, 8.245264573227866, 8.594052082751618, 8.427857315395059, 8.191023946853138, 8.174944978975265, 7.9218012289874515, 7.7531954230023885, 7.714569364377547, 7.772506100881505, 7.748488438231633, 7.672631930424717, 7.4164386849573445, 7.608269365239694, 7.507599235336425, 7.241646933985539, 7.091658274614096, 7.067241941770455, 6.9539568651048755, 7.06095133894324, 6.876611537052968, 6.788599683581194, 6.705610999440179, 6.671532222687652, 6.5435893128274225, 6.452307402357074, 6.511028488016522, 6.419174057781553, 6.226649182401291, 6.227623278227309, 6.172207196211846, 6.007080126001774, 6.2687896327143795, 6.070517315801322, 6.002667621348028, 5.775265190799222, 5.76236455942957, 5.5647916822222445, 5.5382132870191025, 5.48941811405098, 5.398038655061165, 5.19585577365003, 5.203978377093699, 5.160182195056955, 5.102283745605793, 5.0185231474663485, 5.037694981467047, 5.032766472413318, 4.898562807699726, 4.7199345117967635, 4.6934932157427, 4.558512173506229, 4.566143046823501, 4.540458599832846, 4.420582782201622, 4.358386687596291]
    # n = 50
    # alpha = "0.005"
    # q = 2501
    # target_norm = 34536.8621125
    # ebeta = 40
    # jump = 1
    # goal_margin = 1.5



    
    d = len(log_rr0)
    
    
    MAX_TIME = float("inf")

    

    print("usvp simulator in default G6K:\nn = %d, alpha = %s, b= %d, jump = %d, goal_margin = %.2f" %(n,alpha,ebeta,jump,goal_margin))
    file_path = "test"
    if os.path.exists(file_path):
        pass
    else:
        os.mkdir(file_path)
    

    usvp_simulator_default(file_path,n,alpha,log_rr0,target_norm,ebeta,goal_margin,q)


   