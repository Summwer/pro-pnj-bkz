import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt,floor
from scipy import optimize
import numpy as np

def f_1(x, A, B):
    return A * x + B

def compute_R2(actual,predict):
    corr_matrix = np.corrcoef(actual, predict)
    corr = corr_matrix[0,1]  #相关系数
    R_sq = corr**2
    return R_sq
    
    
def draw_A_line(file_name, sieve_dims, log2_costs):
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,4),dpi=600)
        
        x = sieve_dims
        A, B = optimize.curve_fit(f_1, x, log2_costs)[0]
        y = [ A * _ + B for _ in sieve_dims]
        
        R = compute_R2(log2_costs,y)
        ax.scatter(sieve_dims, log2_costs,color = "purple", label=r"A of $T_{\rm mid}$ in dd,$R^2$ = %f" %(R), marker=".", zorder = 10)
        

        if(B<0):
            function =r"$A(n)$ = %.3f $n$ - %.3f, $%d\leq n\leq %d$" %(round(A,3),round(-B,3), sieve_dims[0],sieve_dims[-1])
        else:
            function =r"$A(n)$ = %.3f $n$ + %.3f, $%d\leq n\leq %d$" %(round(A,3),round(B,3),sieve_dims[0],sieve_dims[-1])
        ax.plot(x, y, "red", label=function, zorder = 4, linestyle=":")
        
#        plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        plt.legend(fontsize="9.5")
        ax.set(xlabel=r'$n$')
        ax.set(ylabel=r'$\log_2T (\log_2(sec))$')
        ax.autoscale(tight=False)
        fig.savefig(file_name)
        plt.close()

def draw_cost_model(file_name, sieve_dims, log2_costs,midnode1,midnode2,midnode3,y_name, func_name):
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,4),dpi=600)
        
        
        x1 = sieve_dims[1:(midnode1)]
        py1 = log2_costs[1:(midnode1)]
        A1, B1 = optimize.curve_fit(f_1, x1, py1)[0]
#        R1 = sum([(y1[i]-py[i])**2 for i in range(len(x))])/len(x)
        
        
        
        x2 = [sieve_dims[midnode1-1], sieve_dims[midnode2-1]]
#        py2 = log2_costs[midnode1: midnode2]
        A2 = (log2_costs[midnode2-1] - log2_costs[midnode1-1])/(sieve_dims[midnode2-1] - sieve_dims[midnode1-1])
        B2 = log2_costs[midnode1-1] - A2*sieve_dims[midnode1-1]
##        R2 = sum([(y2[i]-py2[i])**2 for i in range(len(x2))])/len(x2)
#
#
        x3 = sieve_dims[(midnode2-1):(midnode3)]
        py3 = log2_costs[(midnode2-1):(midnode3)]
        A3, B3 = optimize.curve_fit(f_1, x3, py3)[0]

        x4 = sieve_dims[(midnode3):]
        py4 = log2_costs[(midnode3):]
        A4, B4 = optimize.curve_fit(f_1, x4, py4)[0]
#        
#        x00 = (B2-B1)/(A1-A2)
#        x01 = (B3-B2)/(A2-A3)
        x02 = (B4-B3)/(A3-A4)
        x1 = [sieve_dims[0]] + x1
        y1 = [ A1 * _ + B1 for _ in x1]
        
        
        x2 = x2#+[x01]
        y2 = [ A2 * _ + B2 for _ in x2]
#        
        x3 = x3 +[x02]
        y3 = [ A3 * _ + B3 for _ in x3]
#        
        x4 = [x02] +x4
        y4 = [ A4 * _ + B4 for _ in x4]
        
        R = compute_R2(log2_costs,y1[:len(y1)-1]+y2+y3[1:len(y3)-1]+y4[1:])
        ax.scatter(sieve_dims, log2_costs, color = "purple", label=y_name+" in dd,$R^2$ = %f" %(R), marker=".", zorder = 10)
        
        if(B1<0):
            function1 =func_name+" = %.3f $n$ - %.3f, $%d\leq n\leq %d$" %(round(A1,3),round(-B1,3),sieve_dims[0],sieve_dims[midnode1-1])
        else:
            function1 =func_name+" = %.3f $n$ + %.3f, $%d\leq n\leq %d$" %(round(A1,3),round(B1,3),sieve_dims[0],sieve_dims[midnode1-1])
        ax.plot(x1, y1, "red", label=function1, zorder = 4,linestyle="-")
        if(B2<0):
            function2 =func_name+" = %.3f $n$ - %.3f, $%d< n\leq %d$" %(round(A2,3),round(-B2,3),sieve_dims[midnode1-1],sieve_dims[midnode2-1])
        else:
            function2 =func_name+" = %.3f $n$ + %.3f, $%d< n\leq %d$" %(round(A2,3),round(B2,3),sieve_dims[midnode1-1],sieve_dims[midnode2-1])
        ax.plot(x2, y2, "red", label=function2, zorder = 4,linestyle="--")
        if(B3<0):
            function3 =func_name+" = %.3f $n$ - %.3f, $%d< n\leq %d$" %(round(A3,3),round(-B3,3),sieve_dims[midnode2-1],sieve_dims[midnode3-1])
        else:
            function3 =func_name+" = %.3f $n$ + %.3f, $%d< n\leq %d$" %(round(A3,3),round(B3,3),sieve_dims[midnode2-1],sieve_dims[midnode3-1])
        ax.plot(x3, y3, "red", label=function3, zorder = 4,linestyle="-.")
#
        if(B4<0):
            function4 =func_name+" = %.3f $n$ - %.3f, $%d< n\leq %d$" %(round(A4,3),round(-B4,3),sieve_dims[midnode3-1],sieve_dims[-1])
        else:
            function4 =func_name+" = %.3f $n$ + %.3f, $%d< n\leq %d$" %(round(A4,3),round(B4,3),sieve_dims[midnode3-1],sieve_dims[-1])
        ax.plot(x4, y4, "red", label=function4, zorder = 4, linestyle=":")
#
#        ax.scatter(round(x00),A2 * x00 + B2, marker=".",color = "r", zorder = 5)
#        ax.text(round(x00),A2 * x00 + B2,(round(x00),round(A2 * x00 + B2,2)),ha='left', va='top', fontsize=7,zorder = 7)
#        
#        ax.scatter(round(x01),A2 * x01 + B2, marker=".",color = "r",zorder = 6)
#        ax.text(round(x01),A2 * x01 + B2,(round(x01),round(A2 * x01 + B2,2)),ha='left', va='top', fontsize=7,zorder = 6)
#        
#        ax.scatter(round(x02),A3 * x02 + B3, marker=".",color = "r",zorder = 6)
#        ax.text(round(x02),A3 * x02 + B3,(round(x02),round(A3 * x02 + B3,2)),ha='left', va='top', fontsize=7,zorder = 6)
        

#        plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        plt.legend(fontsize="9.5")
        ax.set(xlabel=r'$n$')
        ax.set(ylabel=r'$\log_2T (\log_2(sec))$')
        ax.autoscale(tight=False)
        fig.savefig(file_name)
        plt.close()
        
        


sieve_dims = [46,47,48,49,50,51,52,53,54,55,56,57,58,60,62,64,66,68,69,71,73,75,77,79,81,82,84,86,88,90,92,94,95,97,99]


A_of_AvgMidCosts = [0.00534,0.00764,0.00798,0.00911,0.00951,0.01014,0.01055,0.01104,0.01135,0.01132,0.01635,0.01503,0.01537,0.01542,0.01789,0.01731,0.0204,0.02304,0.02557,0.02061,0.02211,0.03018,0.0228,0.02443,0.03344,0.0312,0.02769,0.0358,0.02453,0.04554,0.02219,0.05102,-0.07373,0.17768,0.08192]


B_of_AvgMidCosts = [2.84042,2.84849,2.95271,2.7975,2.87286,2.96269,3.11718,3.18518,3.31825,3.32144,3.0762,3.37666,3.58225,3.91818,4.24551,4.81582,7.09192,7.85439,8.24816,9.54307,10.30168,10.91955,12.69289,13.86835,15.08675,16.44398,18.85646,22.22364,27.61019,33.01663,43.40242,53.0382,69.92405,77.53214,106.91253]


first_pump_costs = [4.07,2.73,2.79,2.58,2.85,2.9,2.8,3.07,3.06,3.27,3.31,3.17,3.63,3.91,4.38,4.71,7.29,8.78,9.52,10.85,12.06,13.18,15.04,17.79,21.43,21.73,27.44,35.95,44.84,59.26,75.66,97.34,118.09,153.94,211.73]


pnjbkz_pre_costs = [43.32,47.65,52.36,54.33,58.96,64.44,69.49,74.31,79.24,86.36,88.18,95.36,104.18,109.55,121.91,132.37,195.62,220.86,239.29,259.07,284.61,305.28,343.05,379.95,427.57,467.88,542.89,642.26,755.5,927.05,1161.6,1467.14,1822.28,2447.82,3059.83]

pnjbkz_later_costs = [51.15,60.95,66.41,67.83,74.7,80.49,89.16,92.68,104.98,106.37,111.5,129.64,131.89,142.3,155.48,158.97,236.04,251.94,273.5,304.16,327.15,360.35,389.14,412.24,472.97,516.08,577.36,681.03,774.54,964.41,1138.99,1500.69,1681.68,2314.65,3022.42]



log2_first_pump_costs = [log2(_) for _ in first_pump_costs]
log2_pnjbkz_pre_costs = [log2(_) for _ in pnjbkz_pre_costs]
log2_pnjbkz_later_costs = [log2(_) for _ in pnjbkz_later_costs]
log2_B_of_AvgMidCosts = [log2(_) for _ in B_of_AvgMidCosts]


draw_cost_model("pnjbkz_cost_model/log2_first_pump_costs.png", sieve_dims, log2_first_pump_costs,16,17,26,r"$T_{\rm first}$", r"$T_{\rm first}(n)$")

draw_cost_model("pnjbkz_cost_model/log2_pnjbkz_pre_costs.png", sieve_dims, log2_pnjbkz_pre_costs,16,17,26,r"$T_{\rm pre}$", r"$T_{\rm pre}(n)$")

draw_cost_model("pnjbkz_cost_model/log2_pnjbkz_later_costs.png", sieve_dims, log2_pnjbkz_later_costs,16,17,26,r"$T_{\rm later}$", r"$T_{\rm later}(n)$")

draw_cost_model("pnjbkz_cost_model/log2_B_of_AvgMidCosts.png", sieve_dims, log2_B_of_AvgMidCosts,16,17,26,r"B of $T_{\rm mid}$", "B(n)")

draw_A_line("pnjbkz_cost_model/A_of_AvgMidCosts.png", sieve_dims, A_of_AvgMidCosts)
