import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt
from scipy.interpolate import make_interp_spline
import numpy as np
from scipy import optimize
from adjustText import adjust_text
#from mpl_toolkits.mplot3d import Axes3D
#import plotly.graph_objects as go


def text_point(ax, x1, y1):
    texts = []
    i = 1
    for x,y in zip(x1,y1):
        texts.append(ax.text(x,y+1,'%.0f' % y,fontdict={'fontsize':10},  zorder = 10, ha='center',va='bottom'))
        i+=1
    return texts
    
def draw_costs( file_name, ns,alpha, default_g6k_cost,enumbs_cost_mjl1, enumbs_cost_mje1, enumbs_cost_mjl1_sj1, bssa_cost_mjl1, bssa_cost_mje1 , bssa_cost_mjl1_sj1, marker_type, ymax, ylabel_name, xlabel_name):
    plt.rc('font', size=24)
#    fig, ax = plt.subplots(figsize=(9,7),dpi=600)
    fig, ax = plt.subplots(figsize=(9,8),dpi=600)
#    ax.plot(ns,default_g6k_cost,  linewidth =2,markersize=8,marker=marker_type, zorder = 2, color= "black", label = r"default G6K")
#    ax.plot(ns,enumbs_cost_mjl1, linewidth =2,markersize=8, marker=marker_type, zorder = 5, color = 'red', label = r"ProPnJBKZ(EnumBS($J_{max}= d4f_{slope}(s)$))")
#    ax.plot(ns,enumbs_cost_mje1, linewidth =2,markersize=8, marker=marker_type, zorder = 4, color = 'blue', label = r"ProPnJBKZ(EnumBS($J_{max}=1$))")
#    ax.plot(ns,enumbs_cost_mjl1_sj1, linewidth =2,markersize=8, marker=marker_type, zorder = 3, color = 'green', label = r"ProPnJBKZ(EnumBS($J_{max}= d4f_{slope}(s)$)), set $J = 1$ in S")
#    ax.plot(ns,bssa_cost_mjl1, linewidth =2,markersize=8, marker=marker_type, zorder = 8, color = 'orange', label = r"ProPnJBKZ(BSSA($J_{max}= d4f_{slope}(s)$))")
#    ax.plot(ns,bssa_cost_mje1, linewidth =2,markersize=8, marker=marker_type, zorder = 7, color = 'brown', label = r"ProPnJBKZ(BSSA($J_{max}=1$))")
#    ax.plot(ns,enumbs_cost_mjl1_sj1, linewidth =2,markersize=8, marker=marker_type, zorder = 6, color = 'gray', label = r"ProPnJBKZ(BSSA($J_{max}= d4f_{slope}(s)$)), set $J = 1$ in S")
    xticks = np.arange(len(ns))
    
    bars = ax.bar(xticks ,enumbs_cost_mjl1, width = 0.13, zorder = 5, color = '#EA5B58',  label = r"ProPnJBKZ(EnumBS($J_{max}= d4f_{slope}(s)$))", alpha = 1,  edgecolor="k")
    bars += ax.bar(xticks +0.13*1,bssa_cost_mjl1, width =0.13, zorder = 8, color = '#EB8781', label = r"ProPnJBKZ(BSSA($J_{max}= d4f_{slope}(s)$))", alpha = 1, edgecolor="k")
    bars += ax.bar(xticks +0.13*2,enumbs_cost_mje1, width = 0.13, zorder = 10, color = '#719AAC', label = r"ProPnJBKZ(EnumBS($J_{max}=1$))",edgecolor="k")
    bars += ax.bar(xticks +0.13*3,bssa_cost_mje1, width = 0.13, zorder = 7, color = '#E29135', label = r"ProPnJBKZ(BSSA($J_{max}=1$))", alpha = 1,edgecolor="k")
    bars += ax.bar(xticks +0.13*4,enumbs_cost_mjl1_sj1,width = 0.13,zorder = 3, color = '#94C6CD', label = r"ProPnJBKZ(EnumBS($J_{max}= d4f_{slope}(s)$)), set $J = 1$ in S", alpha = 1,edgecolor="k")
    bars += ax.bar(xticks +0.13*5,bssa_cost_mjl1_sj1, width = 0.13, zorder = 6, color = '#909291', label = r"ProPnJBKZ(BSSA($J_{max}= d4f_{slope}(s)$)), set $J = 1$ in S", alpha = 1,edgecolor="k")
    bars += ax.bar(xticks +0.13*6,default_g6k_cost, width = 0.13, zorder = 2, color= "#F0D49C", label = r"default G6K", alpha = 1,edgecolor="k")
    
    
    
    
    texts =text_point(ax, xticks,enumbs_cost_mjl1)
    texts +=text_point(ax, xticks+0.13*1,bssa_cost_mjl1)
    texts +=text_point(ax, xticks+0.13*2,enumbs_cost_mje1)
    texts +=text_point(ax, xticks+0.13*3,bssa_cost_mje1)
    texts +=text_point(ax, xticks+0.13*4,enumbs_cost_mjl1_sj1)
    texts +=text_point(ax, xticks+0.13*5,bssa_cost_mjl1_sj1)
    texts += text_point(ax, xticks+0.13*6,default_g6k_cost)
    
#    adjust_text(texts, add_objects=bars, autoalign=False,
#                  only_move={'points': 'y', 'text': 'xy', 'objects': 'xy'},
#                 ha='center', va='bottom')
#    
    
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.legend(fontsize = 13.5, loc = "upper left")#bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    ax.set_xticks(xticks+0.13*3)
    ax.set_xticklabels(ns)
    plt.xlabel(xlabel_name)# with $\alpha$ = %.3f' %alpha, fontsize = 24)
    plt.ylabel(ylabel_name,fontsize = 24)
#    ax.autoscale(tight=False)
    plt.ylim(0,ymax)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    plt.grid(True)
    fig.savefig(file_name,bbox_inches='tight')
    
    
    

def draw_strategy_gen_costs( ns,alpha,default_g6k_cost,enumbs_cost,bssa_cost, marker_type):
    plt.rc('font', size=24)
    fig, ax = plt.subplots(figsize=(9,7),dpi=600)
        

    ax.plot(ns,enumbs_cost, linewidth =2,markersize=8, marker=marker_type, zorder = 5, color = 'red', label = r"EnumBS(threads = 10)")
    ax.plot(ns,bssa_cost,linewidth =2,markersize=8, marker=marker_type, zorder = 3, color = 'orange', label = r"BSSA" )
    
    
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.legend(fontsize = 20, loc = "upper left")#bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    plt.xlabel(r'$n$')# with $\alpha$ = %.3f' %alpha, fontsize = 24)
    plt.ylabel(r'T/sec',fontsize = 24)
    
    
    
#    ax.autoscale(tight=False)
#    plt.ylim(0,300)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    
    plt.grid(True)
    fig.savefig(r'CostComparison-%.3f.png' %alpha,bbox_inches='tight')
    
    
    
if __name__ == '__main__':

        
        alpha = 0.010
        ns = list(range(55,62))
        #Time cost comparison for lwe instance
        g6k_cost_t = [4092.76,3314.65,5409.96,7211.49,17251.79,15738.62,37882.41]
        enumbs_cost_mjg1_t =  [905.68,892.468,1417.344,1831.302,1811.794,2471.522,3859.664] #max jump>1
        enumbs_cost_mje1_t =  [1604.46,1805.396,3575.234,4078.176,5272.624,6617.398,8447.51] #max jump =1
        enumbs_cost_mjl1_sj1_t = [6329.52,7006.025,6922.92,12940.575,13685.205,18841.495,27976.12]#generate strategy in maxjump>1 but set jump = 1
        bssa_cost_mjg1_t  = [1122.252,1361.206,1788.354,2551.786,2671.528,3235.468,4388.066]
        bssa_cost_mje1_t  = [2904.162,3110.94,3770.016,5509.688,6623.534,7669.176,11115.456]
        bssa_cost_mjl1_sj1_t  = [5155.025,9465.05,9696.205,13563.87,15953.94,18303.285,23200.09]
        draw_costs(r'TimeCostComparison-dd-alpha-%.3f.png' %alpha,  ns,alpha, g6k_cost_t,enumbs_cost_mjg1_t, enumbs_cost_mje1_t, enumbs_cost_mjl1_sj1_t, bssa_cost_mjg1_t, bssa_cost_mje1_t , bssa_cost_mjl1_sj1_t, "s", 40000, r'T/sec',r'$n$')
        
        #Memory cost comparison for LWE instance
        g6k_cost_t = [5.9001,4.6466,6.5407,8.115,12.7987,11.2229,18.3033]
        enumbs_cost_mjg1_t = [5.7222,4.033,14.7288,14.7318,11.4436,12.4288,22.7856]  #max jump>1
        enumbs_cost_mje1_t = [12.2152,11.1204,22.9096,28.1028,10.8872,11.1362,23.7152] #max jump =1
        enumbs_cost_mjl1_sj1_t = [1.965,2.2325,1.98,4.41,5.7675,2.88,3.28] #generate strategy in maxjump>1 but set jump = 1
        bssa_cost_mjg1_t  = [11.2015,6.27514,8.2494,18.7454,5.581,9.7688,28.4878]
        bssa_cost_mje1_t  = [6.467,4.7158,8.4818,16.8776,4.821,7.297,27.617]
        bssa_cost_mjl1_sj1_t  = [7.5625,4.664,4.5115,8.5695,3.8455,3.554,3.769]
        draw_costs(r'MemoryCostComparison-dd-alpha-%.3f.png' %alpha,  ns,alpha, g6k_cost_t,enumbs_cost_mjg1_t, enumbs_cost_mje1_t, enumbs_cost_mjl1_sj1_t, bssa_cost_mjg1_t, bssa_cost_mje1_t , bssa_cost_mjl1_sj1_t, "s", 50, r'RAM/GB',r'$n$')
        

#        alpha = 0.015
#        ns = list(range(48,53))
#        enumbs_cost = [3.0925,4.605,6.3125,10.64,12.575]
#        bssa_cost = [123.24,146.54,166.8025,198.465,225.2425]
#        draw_strategy_gen_costs(ns, alpha, g6k_cost,enumbs_cost,bssa_cost,"s")
#
#
#        alpha = 0.010
#        ns = list(range(54,62))
#        enumbs_cost =  [3.06,4.0875,5.7675,7.4125,10.8525,13.43,20.45,20.131666667]
#        bssa_cost = [89.09,107.2075,129.4375,147.485,178.8725,199.4525,233.2133333,265.33]
#        draw_strategy_gen_costs(ns, alpha, g6k_cost,enumbs_cost,bssa_cost,"s")



        
        ns = ['(40,0.025)', '(50,0.015)', '(45,0.020)', '(40,0.030)']
        #Time cost comparison for lwe challenge
        g6k_cost_t = [6168.33,10333.16,12780.86,80973.86]
        enumbs_cost_mjg1_t =  [861.7525,1344.33,1134.83,4754.68] #max jump>1
        enumbs_cost_mje1_t =  [1575.306,2860.626,2749.842,9816.77] #max jump =1
        enumbs_cost_mjl1_sj1_t = [6110.89,8976.73,9224.72,23182.36]#generate strategy in maxjump>1 but set jump = 1
        bssa_cost_mjg1_t  = [1192.7075,1771.1975,1695.61,7921.53225]
        bssa_cost_mje1_t  = [2420.742,3427.76,3963.926,11755.666]
        bssa_cost_mjl1_sj1_t  = [7719.51,10625.45,11139.39,24477.8]
        draw_costs(r'TimeCostComparison-dd-lwechal.png',  ns,alpha, g6k_cost_t,enumbs_cost_mjg1_t, enumbs_cost_mje1_t, enumbs_cost_mjl1_sj1_t, bssa_cost_mjg1_t, bssa_cost_mje1_t , bssa_cost_mjl1_sj1_t, "s", 85000, r'T/sec',r'$(n,\alpha)$')
        
        #Memory cost comparison for LWE challenge
        g6k_cost_t = [8.2042,7.5408,5.5316,25.0064]
        enumbs_cost_mjg1_t =[9.59275,13.773,7.214,53.771]  #max jump>1
        enumbs_cost_mje1_t = [18.0726,14.1642,17.3322,77.602] #max jump =1
        enumbs_cost_mjl1_sj1_t = [12.795,4.634,4.616,60.317] #generate strategy in maxjump>1 but set jump = 1
        bssa_cost_mjg1_t  = [7.3625,12.7355,6.7705,65.93225]
        bssa_cost_mje1_t  = [10.8546,8.1398,8.9224,85.4864]
        bssa_cost_mjl1_sj1_t  = [10.04,5.48,5.978,35.1]
        draw_costs(r'MemoryCostComparison-dd-lwechal.png',  ns,alpha, g6k_cost_t,enumbs_cost_mjg1_t, enumbs_cost_mje1_t, enumbs_cost_mjl1_sj1_t, bssa_cost_mjg1_t, bssa_cost_mje1_t , bssa_cost_mjl1_sj1_t, "s", 100, r'RAM/GB',r'$(n,\alpha)$')
