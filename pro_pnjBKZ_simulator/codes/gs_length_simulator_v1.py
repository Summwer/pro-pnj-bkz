# coding=utf-8
import matplotlib.pyplot as plt
from math import pi,exp,log,sqrt,lgamma
import os

'''v1: pnj-bkz simulator without jump'''

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


def read_times(file_name):
    f = open(file_name,'r')
    #f = open("50-0010/slopes-%s-50-0010.txt" %(pnj),'r')
    times = f.read()
    f.close()
    times = times.split(' ')
    times.remove('')
    times = [float(_) for _ in times]
    return times

def read_blocksizes(file_name):
    f = open(file_name,'r')
    #f = open("50-0010/slopes-%s-50-0010.txt" %(pnj),'r')
    blocksizes = f.read()
    f.close()
    blocksizes = blocksizes.split(' ')
    blocksizes.remove('')
    blocksizes = [int(_) for _ in blocksizes]
    return blocksizes

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

# Here,the GH is actually GH^2, since rr[i] is b_i^* ^2 .
def calculate_gh(dimension,beta,log_rr0,t):
    log_GH = []
    for i in range(dimension):
        beta_ = min(beta,dimension-i)
        log_gh = (2 * lgamma(beta_/2.0 + 1) + sum([log_rr0[k] for k in range(i,beta_+i)])) / beta_ - log(pi)
        log_GH.append(log_gh)
    
    if(dimension != beta):
        slope =  (log_GH[dimension - beta - 1 ]-log_GH[t+1])/ ( dimension  - t -beta )
    else:
        slope = 1
    # print(t,beta,dimension)
    return log_GH,slope



# Here,the GH is actually GH^2, since rr[i] is b_i^* ^2 .
def calculate_sim_log_gs_lengths(dimension,alpha,tau,beta,f,log_GH,org_log_GS_length,t):
    sim_log_gs_lengths =[]
    if(dimension != beta):
        # slope =  (log_GH[dimension - beta]-log_GH[t+1])/ ( dimension  - t -beta -1)
        slope = max(-(log_GH[dimension - beta ]-log_GH[dimension - beta + 5])/5 , (log_GH[dimension - beta ]-log_GH[dimension - beta - 5])/5 )
        # print(slope)
    else:
        slope = (log_GH[1]-log_GH[0])
    
    for i in range(dimension):
        if(i<=t):
            sim_gs_length = org_log_GS_length[0]
        elif(i>dimension-beta):
            sim_gs_length = sim_log_gs_lengths[dimension-beta]+slope*(i-dimension+beta)
        else:
            log_gh = log_GH[i]
            sim_gs_length = alpha * log_gh + tau
        sim_log_gs_lengths.append(sim_gs_length)

        
    return sim_log_gs_lengths,slope




#show slopes
def show_slope_figure(pnj,dimension,slopes,slope_orgs,blocksizes):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(-0.04, -0.08)
    plt.xticks(ticks = [0]+[_+1 for _ in range(len(blocksizes))], labels =[0]+blocksizes)
    plt.scatter([_ for _ in range(len(slopes))], slopes,marker="*",c = "green")
    #plt.scatter([blocksizes[_] for _ in range(len(slopes))],slopes,marker="*")
    plt.title("slope change for %s , n=50, alpha = 0.005, dimension = %d" %(pnj,dimension))
    plt.show()
    
    
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(-0.04, -0.08)
    plt.xticks(ticks = [0]+[_+1 for _ in range(len(blocksizes))], labels =[0]+blocksizes)
    
    plt.scatter([_ for _ in range(len(slope_orgs))], slope_orgs,marker="*",c = "blue")
    
    #plt.scatter([blocksizes[_] for _ in range(len(slopes))],slopes,marker="*")
    plt.title("slope change for %s , n=50, alpha = 0.005, dimension = %d" %(pnj,dimension))
    plt.show()


def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,slope,blocksize,pnj,f,n,dimension,alpha_,blocksize_start, blocksize_end,square_error,t,k):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(4,17) #set range of y_ticks
    # plt.xlim(-5,210)
    t = 0
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    plt.title("gs-length change for %s , f = %d, n=%d, alpha = %s, dimension = %d, beta = %d, slope = %f, k= %d,square error = %f" %(pnj,f,n,alpha_,dimension,blocksize,slope,k,square_error))
    
    try:
        os.mkdir(dir+"gs-lengths-gh simulator/")
    except FileExistsError:
        pass
    plt.savefig(dir+"gs-lengths-gh simulator/f=%d,n=%d,alpha=%s,k=1,beta=%d,d=%d.png" %(f,n,alpha_,blocksize,dimension))
    plt.show()
    
    
def o_show_gs_slope_figure(log_gs_length,sim_log_gs_lengths,dimension):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(4,17) #set range of y_ticks
    # plt.xlim(-5,210)
    t = 0
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    plt.show()
    
def show_times_figure(dir,times,pnj,f,n,dimension,alpha_,blocksizes):
    plt.figure(figsize=(15, 10), dpi=100)
    plt.scatter(blocksizes,[log(_,2) for _ in times],marker="*")#,c = color)
    plt.scatter(blocksizes,[log(d,2) + 0.24563*beta-14.2521 for beta in blocksizes],marker="*")#,c = color)
    
    plt.title("gs-length change for %s , f = %d, n=%d, alpha = %s, dimension = %d, beta = %d-%d" %(pnj,f,n,alpha_,dimension,blocksizes[0],blocksizes[-1]))
    
    try:
        os.mkdir(dir+"gs-lengths-gh simulator/")
    except FileExistsError:
        pass
    #plt.savefig(dir+"gs-lengths-gh simulator/f=%d,n=%d,alpha=%s,k=1,beta=%d,d=%d.png" %(f,n,alpha_,blocksize,dimension))
    #plt.show()
    
    
def time_cost_predict_beta(dir, times):
    return

#gs-lengths predict with beta change
def find_rule_of_beta(dir,f,n,d,alpha_,k_bound,blocksizes,log_GS_lengths,GS_lengths):
    Sum_square_error = 0
    slope_orgs = []
    slopes = []
    print(len(log_GS_lengths))
    
    #original gs
    show_gs_slope_figure(dir,log_GS_lengths[0],[],0,0,pnj,f,n,len(log_GS_lengths[0]),alpha_,blocksize_start, blocksize_end,0,0,0)
    t_org = 0
    for j in range(1,len(log_GS_lengths[0])):
        if log_GS_lengths[0][0]==log_GS_lengths[0][j]:
            t_org+=1
    print(t_org)
    k = 0
    for i in range(len(blocksizes)):
        t  = -1
        # t = max(round(t_org - 0.35*blocksizes[i]- 9*k**(0.01+0.016*blocksizes[i])),0)
        # if(k<21-blocksizes[i]):
        k+=1
        beta = blocksizes[i]
        dimension = len(log_GS_lengths[2*i+1])
        
        # t = 0
        # for j in range(1,len(log_GS_lengths[2*i+1])):
        #     if log_GS_lengths[2*i+1][0]==log_GS_lengths[2*i+1][j]:
        #         t+=1
        
        slope_org =  (log_GS_lengths[2*i+1][dimension - beta - 1 ]-log_GS_lengths[2*i+1][t+1])/ (dimension -beta-t)

        log_GH,GH_slope = calculate_gh(dimension, blocksizes[i],log_GS_lengths[2*i],t)
        # alpha1 = slope_org/GH_slope #prams[0] - prams[1]  *( dimension + beta- f)
        # alpha = 0.8 - 0.001*blocksizes[i]+ 0.15*k**(0.003*blocksizes[i]**1.2)
        alpha =1
        tau =0
        # print(alpha1,alpha1_)
        # tau =  -(alpha * log_GH[t]- log_GS_lengths[2*i+1][t]) #prams[2]  +prams[3] * ( dimension-beta+f)
        # tau = 0.965 + 0.028 * (blocksizes[i]-2*k)
        # tau = log_GS_lengths[2*i][0] - alpha * log_GH[t+1]
        # print(tau1,tau1_)
        # alpha2 = 0
        # tau2 = 0
        #alpha,tau1,tau2 = 0.8 -0.006 * i,  2.3 +0.1 * i, 1 - 0.03 *  i #f=10,n=50,alpha=0.005,k=1,beta=61-70
        
        sim_log_gs_lengths,slope = calculate_sim_log_gs_lengths(dimension,alpha,tau, blocksizes[i],f,log_GH,log_GS_lengths[2*i],t)
        print(slope_org,slope)
        
        slope_orgs.append(slope_org/GH_slope)
        # slopes.append(slope)

        
        # print(slope_org/GH_slope , slope_org/GH_slope*log_GH[1]- log_GS_lengths[2*i+1][1])

        # sim_log_gs_lengths = log_GH
        # slope = GH_slope
        #calculate the square error from simulator to log_GS_lengths.
        square_error = compute_square_error(sim_log_gs_lengths,log_GS_lengths[2*i+1],flag)
        Sum_square_error += square_error
        
        #show_gs_slope_figure(log_GS_lengths[i+1],log_GH,slope,blocksizes[i],pnj,f,n,dimension,alpha_,blocksize_start, blocksize_end)
        show_gs_slope_figure(dir,log_GS_lengths[2*i+1],sim_log_gs_lengths,slope,blocksizes[i],pnj,f,n,dimension,alpha_,blocksize_start, blocksize_end,square_error,t,k)
    
    print(round(Sum_square_error/(len(blocksizes)),5) )
    
    # show_slope_figure(pnj,dimension,slopes,slope_orgs,blocksizes)
    
# print(round(Sum_square_error,5) )


#dimension change
def find_rule_of_dimension(dir,f,n,alpha_,k,blocksize,pram_dict,d_start,d_end):
    dir = "pro-bkz-tests/gs-lengths-simulator/f=%d,n=%d,alpha=%s,k=%d,beta=%d,d=%d-%d/" %(f,n,alpha_,k,blocksize,d_start,d_end)
    blocksizes = read_blocksizes(dir+"blocksizes.txt" )
    log_GS_lengths,GS_lengths = read_gs_lengths(dir+"rr_set.txt" )
    Sum_square_error = 0
    for i in range(len(log_GS_lengths)//2):
        dimension = len(log_GS_lengths[2*i])
        prams = pram_dict['f=%d,n=%d,alpha=%s,k=%d,beta=%d,d=%d-%d' %(f,n,alpha_,k,blocksize,d_start,d_end )]
        alpha1 = prams[0] - prams[1]  *( dimension + blocksize- f)
        tau1 = prams[2]  +prams[3] * ( dimension-blocksize+f)
        alpha2 = 0
        tau2 = 0
        alpha1,tau1,alpha2,tau2 = prams[0] - prams[1]  * i,  prams[2]  +prams[3] * i, prams[4]  - prams[5]  *  i,prams[6]  - prams[7]  *  i #f=10,n=40,alpha=0.005,k=1,beta=61-70
        #alpha,tau1,tau2 = 0.8 -0.006 * i,  2.3 +0.1 * i, 1 - 0.03 *  i #f=10,n=50,alpha=0.005,k=1,beta=61-70
        log_GH,GH_slope = calculate_gh(dimension,alpha1,tau1,alpha2,tau2, blocksizes[0],f,GS_lengths[2*i])
        
        sim_log_gs_lengths,slope = calculate_sim_log_gs_lengths(dimension,alpha1,tau1,alpha2,tau2, blocksizes[0],f,log_GH)
        #calculate the square error from simulator to log_GS_lengths.
        square_error = compute_square_error(sim_log_gs_lengths,log_GS_lengths[2*i+1],flag)
        Sum_square_error += square_error
        
        #show_gs_slope_figure(log_GS_lengths[i+1],log_GH,slope,blocksizes[i],pnj,f,n,dimension,alpha_,blocksize_start, blocksize_end)
        show_gs_slope_figure(dir,log_GS_lengths[2*i+1],sim_log_gs_lengths,slope,blocksize,pnj,f,n,dimension,alpha_,blocksize, blocksize,square_error)
    

    print(round(Sum_square_error/(len(log_GS_lengths)//2),5) )
    
    

#=========================================input==============================================

if __name__ == "__main__":

    pnj = 'pnj-bkz'
    #pnj = 'bkz'

    flag = 1 #remove error in square error(1) or not (0)



    pram_dict = {}



    pram_dict['f=0,n=40,alpha=0.005,k=1,beta=51-60'] = (0.8, 0.001,  3.2 , 0.006, 0, 0,0, 0)
    pram_dict['f=0,n=40,alpha=0.01,k=1,beta=51-60'] = (0.9, 0.001,  3.2 , 0.006, 0, 0,0, 0)
    pram_dict['f=0,n=50,alpha=0.005,k=1,beta=51-60'] = (0.8, 0.000,  -2.8688638651511944 , 0, 0, 0,0, 0)
    pram_dict['f=0,n=50,alpha=0.01,k=1,beta=10-66'] = (0.8,- 100,  -2 , 1, 0, 0,0, 0)

    pram_dict['f=10,n=40,alpha=0.005,k=1,beta=61-70'] =(0.8, 0.001,  3.2 , 0.006, 0, 0,0, 0)
    pram_dict['f=10,n=40,alpha=0.01,k=1,beta=61-70'] =(0.9, 0.001,  3.2 , 0.006, 0, 0,0, 0)
    pram_dict['f=10,n=50,alpha=0.005,k=1,beta=61-70'] = (0.8, 0.000,  2.8688638651511944 , 0, 0, 0,0, 0)     #(0.8, 0.005,  2.3, 0.1, 1, 0.04)
    pram_dict['f=10,n=50,alpha=0.01,k=1,beta=61-70'] = (0.9, 0.001,  3.2 , 0.006, 0, 0,0, 0)


    # pram_dict['f=15,n=50,alpha=0.005,k=1,beta=66-75'] = (0.8, 0.005,  2.3, 0.08, 1.3, 0.04)
    pram_dict['f=15,n=50,alpha=0.005,k=1,beta=75-75'] = (0.9, 0.001,  3.2 , 0.006, 0, 0,0, 0)
    pram_dict['f=15,n=50,alpha=0.005,k=1,beta=75-75'] = (0.9, 0.001,  3.2 , 0.006, 0, 0,0, 0)

    pram_dict['f=15,n=50,alpha=0.005,k=1,beta=70'] = (0.75, 0,  3.5, 0, 0 , 0, 0, 0)
    pram_dict['f=10,n=50,alpha=0.015,k=1,beta=70,d=151'] = (0.75, 0,  3.5, 0, 0 , 0, 0, 0)
    pram_dict['f=10,n=50,alpha=0.015,k=1,beta=70,d=151-200'] = (0.75, 0.001,  3.5, 0.006, 0 , 0, 0, 0)



    # #-------------------input:f=0, n=50, d=184, alpha=0.010, k=1, beta=10-66-------------------------------
    # f, n, d, alpha_, k, blocksize_start, blocksize_end = 0,50,184,"0.01",1,10,66
    # fs = [0]
    # log_GS_lengths,GS_lengths,blocksizes,times,dir = read_file(n,alpha_,k,fs,blocksize_start,blocksize_end) #read files

    # find_rule_of_beta(dir,f,n,d,alpha_,k,blocksizes,log_GS_lengths,GS_lengths)

    # for i in range(len(blocksizes)):
    #     beta = blocksizes[i]
    #     print(beta,times[i], d * 2**(0.24563*beta-14.2521))

    # show_times_figure(dir,times,pnj,f,n,d,alpha_,blocksizes)
    
    
    
    
    # #-------------------input:f=0, n=50, d=184, alpha=0.010, k=1, beta=10-66-------------------------------
    # f, n, d, alpha_, k, blocksize_start, blocksize_end = 0,50,184,"0.01",50,10,10
    # fs = [0]
    # log_GS_lengths,GS_lengths,blocksizes,times,dir = read_file(n,alpha_,k,fs,blocksize_start,blocksize_end) #read files

    # find_rule_of_beta(dir,f,n,d,alpha_,k,blocksizes,log_GS_lengths,GS_lengths)

    # #-------------------input:f=0, n=50, d=184, alpha=0.010, k=1, beta=10-66-------------------------------
    # f, n, d, alpha_, k, blocksize_start, blocksize_end = 0,50,184,"0.01",20,20,20
    # fs = [0]
    # log_GS_lengths,GS_lengths,blocksizes,times,dir = read_file(n,alpha_,k,fs,blocksize_start,blocksize_end) #read files

    # find_rule_of_beta(dir,f,n,d,alpha_,k,blocksizes,log_GS_lengths,GS_lengths)
    
    
    
      # -------------------input:f=0, n=50, d=184, alpha=0.010, k=1, beta=10-66-------------------------------
    f, n, d, alpha_, k, blocksize_start, blocksize_end = 0,50,184,"0.01",10,50,50
    fs = [0]
    log_GS_lengths,GS_lengths,blocksizes,times,dir = read_file(n,alpha_,k,fs,blocksize_start,blocksize_end) #read files

    find_rule_of_beta(dir,f,n,d,alpha_,k,blocksizes,log_GS_lengths,GS_lengths)
