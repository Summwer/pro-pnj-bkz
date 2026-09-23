import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt
from scipy.interpolate import make_interp_spline
import numpy as np
from scipy import optimize


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


#theoretical f1
def theoretical_dim4free_wrapper1(blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0

    return int(blocksize*log(4/3.)/log(blocksize/(2*pi)))


#theoretical f2
def theoretical_dim4free_wrapper2(blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0

    return int(blocksize*log(4/3.)/log(blocksize/(2*pi*e)))



def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)



def f_1(x, A, B):
    return A * x + B
    
    
def compute_R2(actual,predict):
    
# actual = [1,2,3,4,5]
# predict = [1,2.5,3,4.9,4.9]

    corr_matrix = np.corrcoef(actual, predict)
    corr = corr_matrix[0,1]  #相关系数
    R_sq = corr**2
    return R_sq


def draw_practical_pump_cost():
#    qd_type = True
    #qd float_type
    if(True):
        x = list(range(75,155,5))
        x = [_-dim4free_wrapper(default_dim4free_fun,_) for _ in x]
        # py = [23.62,31.63,42.53,48.81,62.14,71.36,98.07,139.11,204.22,369.88,762.16,1487.49,3485.89,7359.20,21271.08,64823.46,299970.496]
        py = [23.62,31.63,42.53,48.81,62.14,71.36,98.07,139.11,204.22,369.88,762.16,1487.49,3485.89,7359.20,21271.08,64823.46]
        print(x) #[58, 63, 68, 72, 77, 81, 86, 91, 95, 100, 105, 109, 114, 118, 123, 128]
        py = [log(py[i],2) for i in range(len(py))]
        
#        RAM_x = list(range(120,155,5))
#        RAM_x = [_-theoretical_dim4free_wrapper1(_) for _ in RAM_x]
#        RAM_py = [4.38,6.77,10.44,19.09,32.31,61.35,124.67]
#        RAM_py = [log(RAM_py[i],2) for i in range(len(RAM_py))]
        
        thre1 = 8
        thre2 = 12
        
    
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,3),dpi=600)
        
        #ax.plot(x, y, label="pnj-BKZ Cost")
        
        
        A1, B1 = optimize.curve_fit(f_1, x[:thre1+1], py[:thre1+1])[0]
        A2, B2 = optimize.curve_fit(f_1, x[thre1+1:thre2+1], py[thre1+1:thre2+1])[0]
        A3, B3 = optimize.curve_fit(f_1, x[thre2+1:], py[thre2+1:])[0]
        # A3, B3 = optimize.curve_fit(f_1, x[thre2+1:], py[thre2+1:])[0]
        A4 = 0.368
        B4 = py[-1] - A4 * x[-1]
        
        A1 = round(A1,3)
        B1 = round(B1,2)
        A2 = round(A2,3)
        B2 = round(B2,2)
        A3 = round(A3,3)
        B3 = round(B3,2)
        A4 = round(A4,3)
        B4 = round(B4,2)
        
        
        x0 = (B2-B1)/(A1-A2)
        x01 = (B3-B2)/(A2-A3)
        
        x1 = x[:thre1+1]+[x0] #30和75要对应x0的两个端点，0.01为步长
        y1 = [ A1 * _ + B1 for _ in x1]
        
        # R = sum([(y1[i]-py[i])**2 for i in range(len(x1)-1)])/len(x)
        function1 =r"$\log_2T_{\rm Pump}$ = %.3f $n$  %.2f" %(A1,B1)
        
    
        
        x2 = [x0] + x[thre1+1:thre2+1]+[x01] #30和75要对应x0的两个端点，0.01为步长
        y2 = [ A2 * _ + B2 for _ in x2]
        
        # R = R + sum([(y2[i+1]-py[i+thre1+1])**2 for i in range(len(x2)-1)])/len(x)
        function2 =r"$\log_2T_{\rm Pump}$ = %.3f $n$ %.2f" %(A2,B2)
        
        ax.scatter(round(x0),A2 * x0 + B2, marker=".",color = "r")
        ax.text(round(x0),A2 * x0 + B2,(round(x0),round(A2 * x0 + B2,2)),ha='right', va='bottom', fontsize=7)
        
        x3 = [x01] + x[thre2+1:] #30和75要对应x0的两个端点，0.01为步长
        y3 = [ A3 * _ + B3 for _ in x3]
        
        
        # R = R + sum([(y3[i+1]-py[i+thre2+1])**2 for i in range(len(x3)-1)])/len(x)
        function3 =r"$\log_2T_{\rm Pump}$ = %.3f $n$ %.2f" %(A3,B3)
        
        x4 = [_-dim4free_wrapper(default_dim4free_fun,_) for _ in range(150,155,2)]
        y4 = [ A4 * _ + B4 for _ in x4]
        function4 =r"$\log_2T_{\rm Pump}$ = %.3f $n$ %.2f" %(A4,B4)
        
        ax.scatter(round(x01),A2 * x01 + B2, marker=".",color = "r")
        ax.text(round(x01),A2 * x01 + B2,(round(x01),round(A2 * x01 + B2,2)),ha='right', va='bottom', fontsize=7)
        
        ax.text(x3[-1],y3[-1],(round(x3[-1]),round(y3[-1],2)),ha='right', va='bottom', fontsize=7)
        # ax.text(x4[-1],y4[-1],(round(x4[-1]),round(y4[-1],2)),ha='left', va='top', fontsize=7)
        
        R = compute_R2(py,y1[:len(y1)-1]+y2[1:len(y2)-1]+y3[1:])
#        if draw_pump == "test7.down_sieve=False,f=g6k.png":
#            scatter = r"$\log_2T_{\rm Pump}$,$R^2$ = %f, down_sieve = False, $f$ in G6K" %(R)
#        elif draw_pump == "test7.down_sieve=True,f=theo1.png":
#            scatter = r"$\log_2T_{\rm Pump}$,$R^2$ = %f, down_sieve = True, $f = \frac{n \ln 4/3}{\ln(n/2\pi)}$" %(R)
#        elif draw_pump == "test7.down_sieve=True,f=g6k.png":
        scatter_label = r"$\log_2T_{\rm Pump}$ in qd,$R^2$ = %f" %(R)
        
        ax.scatter(x, py, label=scatter_label, marker=".",color='black')
        ax.plot(x1, y1, "#696969",linestyle = "-.", label=function1)
        ax.plot(x2, y2, "#696969",linestyle = "-.",label=function2)
        ax.plot(x3, y3, "#696969",linestyle = "-.", label=function3)
        ax.plot(x4, y4, "#696969",linestyle = "-.", label=function4)
    

    

        #dd float_type
        sieve_dim = list(range(51,133,2))
        log2_pump_cost = [3.304511042,3.304511042,3.339137385,3.415488271,3.547203025,3.642701572,3.753818443,4.110196178,4.275007047,4.475084883,4.552745937,4.618825953,4.723012396,4.797012978,4.885086225,5.032100843,5.163498732,5.364572432,5.699607154,5.976821852,6.333960351,6.722739236,7.047342073,7.451705665,7.848935856,8.248781488,8.665016072,9.115407855,9.529430554,9.983777435,10.45878595,10.8985107,11.40961899,11.9339451,12.44607526,13.02799236,13.55211309,14.15709504,14.75927237,15.50735209,16.30194036]
        
        
        midnode1 = 87
        midnode2 = 111
        midnode3 = 125
    
#    with plt.style.context(['seaborn-white',"science"]):
#        fig, ax = plt.subplots(figsize=(5,3),dpi=600)
        
        
        
        x = sieve_dim[:(midnode1-51)//2]
        py = log2_pump_cost[:(midnode1-51)//2]
        A1, B1 = optimize.curve_fit(f_1, x, py)[0]
#        R1 = sum([(y1[i]-py[i])**2 for i in range(len(x))])/len(x)

        x2 = sieve_dim[(midnode1-51)//2: (midnode2-51)//2]
        py2 = log2_pump_cost[(midnode1-51)//2: (midnode2-51)//2]
        A2, B2 = optimize.curve_fit(f_1, x2, py2)[0]
#        R2 = sum([(y2[i]-py2[i])**2 for i in range(len(x2))])/len(x2)

        x3 = sieve_dim[(midnode2-51)//2:(midnode3-51)//2]
        py3 = log2_pump_cost[(midnode2-51)//2:(midnode3-51)//2]
        A3, B3 = optimize.curve_fit(f_1, x3, py3)[0]
        
        x4 = sieve_dim[(midnode3-51)//2:]
        py4 = log2_pump_cost[(midnode3-51)//2:]
        A4, B4 = optimize.curve_fit(f_1, x4, py4)[0]
        
        x00 = (B2-B1)/(A1-A2)
        x01 = (B3-B2)/(A2-A3)
        x02 = (B4-B3)/(A3-A4)
        x1 = x+[x00]
        y1 = [ A1 * _ + B1 for _ in x1]
        
        
        x2 = [x00]+x2+[x01]
        y2 = [ A2 * _ + B2 for _ in x2]
        
        x3 = [x01] +x3 +[x02]
        y3 = [ A3 * _ + B3 for _ in x3]
        
        x4 = [x02] +x4
        y4 = [ A4 * _ + B4 for _ in x4]
        
        
#        print(x1)
#
#        print(x2)
        
        R = compute_R2(log2_pump_cost,y1[:len(y1)-1]+y2[1:len(y2)-1]+y3[1:len(y3)-1]+y4[1:])
        
        ax.scatter(sieve_dim, log2_pump_cost, label="log2(Pump Cost) in dd,$R^2$ = %f" %(R), marker=".", zorder = 3)
        
        if(B1<0):
            function1 =r"$T_{\rm Pump}$ = %f $n$ - %f" %(A1,-B1)
        else:
            function1 =r"$T_{\rm Pump}$ = %f $n$ + %f" %(A1,B1)
        ax.plot(x1, y1, "blue", label=function1, zorder = 4)
        if(B2<0):
            function2 =r"$T_{\rm Pump}$ = %f $n$ - %f" %(A2,-B2)
        else:
            function2 =r"$T_{\rm Pump}$ = %f $n$ + %f" %(A2,B2)
        ax.plot(x2, y2, "blue", label=function2, zorder = 4)
        if(B3<0):
            function3 =r"$T_{\rm Pump}$ = %f $n$ - %f" %(A3,-B3)
        else:
            function3 =r"$T_{\rm Pump}$ = %f $n$ + %f" %(A3,B3)
        ax.plot(x3, y3, "blue", label=function3, zorder = 4)
        
        if(B4<0):
            function4 =r"$T_{\rm Pump}$ = %f $n$ - %f" %(A4,-B4)
        else:
            function4 =r"$T_{\rm Pump}$ = %f $n$ + %f" %(A4,B4)
        ax.plot(x4, y4, "blue", label=function4, zorder = 4)
        
        ax.scatter(round(x00),A2 * x00 + B2, marker=".",color = "r", zorder = 5)
        ax.text(round(x00),A2 * x00 + B2,(round(x00),round(A2 * x00 + B2,2)),ha='left', va='top', fontsize=7,zorder = 7)
        
        ax.scatter(round(x01),A2 * x01 + B2, marker=".",color = "r",zorder = 6)
        ax.text(round(x01),A2 * x01 + B2,(round(x01),round(A2 * x01 + B2,2)),ha='left', va='top', fontsize=7,zorder = 6)
        
        ax.scatter(round(x02),A3 * x02 + B3, marker=".",color = "r",zorder = 6)
        ax.text(round(x02),A3 * x02 + B3,(round(x02),round(A3 * x02 + B3,2)),ha='left', va='top', fontsize=7,zorder = 6)
        
#
#
#        ax.legend([r"$T_{\rm Pump}(beta=80)$", function1,r"$T_{\rm Pump}(beta=90)$", function2],loc = (0,1))
#        ax.legend(frameon = True)
        plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        ax.set(xlabel=r'$n$')
        ax.set(ylabel=r'$\log_2T/\log_2(sec)$')
        ax.autoscale(tight=False)
        fig.savefig(r'PumpCost-dd.png')
        plt.close()
        
        
        
        

def draw_mulvar_pump_cost():
    x = list(range(75,155,5))
    x = [_-dim4free_wrapper(default_dim4free_fun,_) for _ in x]
    # py = [23.62,31.63,42.53,48.81,62.14,71.36,98.07,139.11,204.22,369.88,762.16,1487.49,3485.89,7359.20,21271.08,64823.46,299970.496]
    py = [23.62,31.63,42.53,48.81,62.14,71.36,98.07,139.11,204.22,369.88,762.16,1487.49,3485.89,7359.20,21271.08,64823.46]
    print(x) #[58, 63, 68, 72, 77, 81, 86, 91, 95, 100, 105, 109, 114, 118, 123, 128]
    
    
    fit_y = [0.05 * 2**(0.367*n-27.108)+0.1*n*2**(0.2075*n-16.187) for n in x]
    
    R = compute_R2(fit_y, py)
    
    fit_x = np.arange(x[0],x[-1]+0.1, 0.1)
    fit_y = 0.05 * 2**(0.367*fit_x-27.108)+0.1*fit_x*2**(0.2075*fit_x-16.187)
    
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(7,3),dpi=600)
        
        ax.scatter(x, py, label="Pump Cost in dd,$R^2$ = %f" %(R), marker="x", color="purple", zorder = 3)
        ax.plot(fit_x, fit_y, "red", linestyle = "--", label=r"SimTpump($n$) = $0.05 \times 2^{0.367n-27.108}+0.1n\times 2^{0.2075n-16.187}$", zorder = 2)
        plt.legend(borderaxespad=1)
        #plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y') #
        ax.set(xlabel=r'$n$')
        ax.set(ylabel=r'walltime/sec')
        ax.autoscale(tight=False)
        fig.savefig(r'SimTpump-dd.png')
        plt.close()
        
        
def draw_practical_pnjbkz_cost(blocksizes, jumps, sieve_dimes, costs,file_name,thre1,thre2,thre3):
#    indices_num = 2*(floor((f+extra_dim4free)/jump)+1) + (floor((d-beta-extra_dim4free)/jump)+1) + 1
    
    sieve_dims_jump1= []
    costs_jump1 = []
    for i in range(len(blocksizes)):
        blocksize,jump,sieve_dim,cost = blocksizes[i], jumps[i], sieve_dims[i], costs[i]
        if(jump==1):
            sieve_dims_jump1.append(sieve_dim)
            costs_jump1.append(cost)
            
    py = [log2(_) for _ in costs_jump1]
    x = sieve_dims_jump1
    
    print(py)
        
        
    
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,3),dpi=600)

        #ax.plot(x, y, label="pnj-BKZ Cost")


        A1, B1 = optimize.curve_fit(f_1, x[:thre1+1], py[:thre1+1])[0]
        A2, B2 = optimize.curve_fit(f_1, x[thre1+1:thre2+1], py[thre1+1:thre2+1])[0]
        
        
        A3, B3 = optimize.curve_fit(f_1, x[thre2+1:thre3+1], py[thre2+1:thre3+1])[0]
        A4, B4 = optimize.curve_fit(f_1, x[thre3+1:], py[thre3+1:])[0]
        
        A1 = round(A1,3)
        B1 = round(B1,2)
        A2 = round(A2,3)
        B2 = round(B2,2)
        A3 = round(A3,3)
        B3 = round(B3,2)
        A4 = round(A4,3)
        B4 = round(B4,2)


        x0 = (B2-B1)/(A1-A2)
        x01 = (B3-B2)/(A2-A3)
        x02 = (B4-B3)/(A3-A4)

        x1 = x[:thre1+1]+[x0] #30和75要对应x0的两个端点，0.01为步长
        y1 = [ A1 * _ + B1 for _ in x1]

        # R = sum([(y1[i]-py[i])**2 for i in range(len(x1)-1)])/len(x)
        function1 =r"$\log_2T_{\rm PnjBKZ}$ = %.3f $n$ + %.2f" %(A1,B1)



        x2 = [x0] + x[thre1+1:thre2+1] #+[x01] #30和75要对应x0的两个端点，0.01为步长
        y2 = [ A2 * _ + B2 for _ in x2]

        # R = R + sum([(y2[i+1]-py[i+thre1+1])**2 for i in range(len(x2)-1)])/len(x)
        function2 =r"$\log_2T_{\rm PnjBKZ}$ = %.3f $n$ + %.2f" %(A2,B2)

        ax.scatter(round(x0),A2 * x0 + B2, marker=".",color = "r")
        ax.text(round(x0),A2 * x0 + B2,(round(x0),round(A2 * x0 + B2,2)),ha='right', va='bottom', fontsize=7)
        
        
        ax.scatter(x2[-1],A2 * x2[-1] + B2, marker=".",color = "r")
        ax.text(x2[-1],A2 * x2[-1] + B2,(round(x2[-1]),round(A2 * x2[-1] + B2,2)),ha='right', va='bottom', fontsize=7)

        x3 = x[thre2+1:thre3+1] + [x02] #30和75要对应x0的两个端点，0.01为步长 #[x01] +
        y3 = [ A3 * _ + B3 for _ in x3]
        
        ax.scatter(x3[0],A3 * x3[0] + B3, marker=".",color = "r")
        ax.text(x3[0],A3 * x3[0] + B3,(x3[0],round(A3 * x3[0] + B3,2)),ha='right', va='bottom', fontsize=7)


        # R = R + sum([(y3[i+1]-py[i+thre2+1])**2 for i in range(len(x3)-1)])/len(x)
        function3 =r"$\log_2T_{\rm PnjBKZ}$ = %.3f $n$ + %.2f" %(A3,B3)

#        x4 = [_-dim4free_wrapper(default_dim4free_fun,_) for _ in range(150,155,2)]
#        y4 = [ A4 * _ + B4 for _ in x4]
        x4 = [x02]+x[thre3+1:] #30和75要对应x0的两个端点，0.01为步长 #[x01] +
        y4 = [ A4 * _ + B4 for _ in x4]
        function4 =r"$\log_2T_{\rm PnjBKZ}$ = %.3f $n$ %.2f" %(A4,B4)

        ax.scatter(round(x02),A3 * x02 + B3, marker=".",color = "r")
        ax.text(round(x02),A3 * x02 + B3,(round(x02),round(A3 * x02 + B3,2)),ha='right', va='bottom', fontsize=7)
        
        
        
        ax.scatter(round(x4[-1]),A4 * x4[-1] + B4, marker=".",color = "r")
        ax.text(round(x4[-1]),A4 * x4[-1] + B4,(round(x4[-1]),round(A4 * x4[-1] + B4,2)),ha='right', va='bottom', fontsize=7)
        
#        ax.text(x3[-1],y3[-1],(round(x3[-1]),round(y3[-1],2)),ha='right', va='bottom', fontsize=7)
        # ax.text(x4[-1],y4[-1],(round(x4[-1]),round(y4[-1],2)),ha='left', va='top', fontsize=7)

        R = compute_R2(py,y1[:len(y1)-1]+y2[0:len(y2)-1]+y3[1:len(y3)-1]+y4)
#        if draw_pump == "test7.down_sieve=False,f=g6k.png":
#            scatter = r"$\log_2T_{\rm Pump}$,$R^2$ = %f, down_sieve = False, $f$ in G6K" %(R)
#        elif draw_pump == "test7.down_sieve=True,f=theo1.png":
#            scatter = r"$\log_2T_{\rm Pump}$,$R^2$ = %f, down_sieve = True, $f = \frac{n \ln 4/3}{\ln(n/2\pi)}$" %(R)
#        elif draw_pump == "test7.down_sieve=True,f=g6k.png":
        scatter_label = r"$\log_2T_{\rm PnjBKZ}$ in dd,$R^2$ = %f" %(R)

        ax.scatter(x, py, label=scatter_label, marker=".",color='black')
        ax.plot(x1, y1, "#696969",linestyle = ":", label=function1)
        ax.plot(x2, y2, "#696969",linestyle = "-.",label=function2)
        ax.plot(x3, y3, "#696969",linestyle = "-", label=function3)
        ax.plot(x4, y4, "#696969",linestyle = "--", label=function4)



#
#
#        ax.legend([r"$T_{\rm Pump}(beta=80)$", function1,r"$T_{\rm Pump}(beta=90)$", function2],loc = (0,1))
#        ax.legend(frameon = True)
        plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
        ax.set(xlabel=r'$n$')
        ax.set(ylabel=r'$\log_2T/\log_2(sec)$')
        ax.autoscale(tight=False)
        fig.savefig(file_name)

        plt.close()

    

draw_practical_pump_cost()
draw_mulvar_pump_cost()

##Default PnjBKZ
#blocksizes = [51,53,53,55,55,57,57,57,59,59,59,61,61,61,61,63,63,63,63,65,65,65,65,65,67,67,67,67,67,69,69,69,69,69,71,71,71,71,71,73,73,73,73,73,75,75,75,75,75,77,77,77,77,77,79,79,79,79,79,81,81,81,81,81,83,83,83,83,83,85,85,85,85,85,87,87,87,87,87,89,89,89,89,89,91,91,91,91,91,93,93,93,93,93,95,95,95,95,95,97,97,97,97,97,99,99,99,99,99,101,101,101,101,101,103,103,103,103,103,105,105,105,105,105,107,107,107,107,107,109,109,109,109,109,111,111,111,111,111,113,113,113,113,113]
#jumps = [1,1,2,1,2,1,2,3,1,2,3,1,2,3,4,1,2,3,4,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]
#    
#sieve_dims = [46,47,47,48,48,49,49,49,50,50,50,51,51,51,51,52,52,52,52,53,53,53,53,53,54,54,54,54,54,55,55,55,55,55,56,56,56,56,56,57,57,57,57,57,58,58,58,58,58,60,60,60,60,60,62,62,62,62,62,64,64,64,64,64,66,66,66,66,66,68,68,68,68,68,69,69,69,69,69,71,71,71,71,71,73,73,73,73,73,75,75,75,75,75,77,77,77,77,77,79,79,79,79,79,81,81,81,81,81,82,82,82,82,82,84,84,84,84,84,86,86,86,86,86,88,88,88,88,88,90,90,90,90,90,92,92,92,92,92,94,94,94,94,94]
#    
#costs = [218.05,229.4,116.93,232.58,122.94,237.2,131.26,86.59,238.03,129.82,90.37,252.46,137.9,95.58,73.1,263.25,143.89,100.56,79.26,281.19,151.44,104.08,80.13,63.39,293.71,160.74,113.71,87.06,66.89,313.32,170.89,117.64,90.5,74.08,330.4,179.66,120.32,95.56,77.26,348.31,184.62,132.85,99.49,79.56,377.44,200.6,137.9,114.05,83.37,396.36,213.25,147.6,120.35,93.27,429.7,231.55,163.11,128.47,102.83,470.77,262.41,175.23,143.59,117.63,776.68,411.47,288.62,234.65,181.17,872.69,468.47,323.73,260.63,207.35,931.81,501.62,341.9,278.76,227.48,1048.23,548.54,385.23,301.74,247.07,1109.48,589.46,411.34,326.97,266.69,1202.82,636.19,438.65,352.46,284.38,1276.11,679.22,479.17,387.07,313.15,1403.7,742.55,525.36,415.59,342.83,1533.72,830.1,568.75,468.65,374.76,1649.36,894.69,627.54,484.21,422,1879.88,1011.11,706.32,568.65,473.59,2200.37,1187.36,842,655.73,578.22,2652.25,1405.08,1023.47,806.99,699.85,3159.65,1741.36,1234.88,977.26,836.27,4068.67,2248.07,1576.71,1277.9,1098.35,5177.49,2829.48,2029.69,1548.8,1458.53]
#    
#thre1 = 4
#thre2 = 15
#thre3 = 25
#draw_practical_pnjbkz_cost(blocksizes, jumps, sieve_dims, costs, 'PnjBKZCost-dd.png', thre1, thre2, thre3)
#
#
##PnjBKZ with pump_params["down_stop"] = blocksize - dim4free
#blocksizes = [51, 51, 51, 51, 51, 51, 53, 53, 53, 53, 53, 53, 55, 55, 55, 55, 55, 55, 57, 57, 57, 57, 57, 57, 59, 59, 59, 59, 59, 59, 61, 61, 61, 61, 61, 61, 63, 63, 63, 63, 63, 63, 65, 65, 65, 65, 65, 65, 67, 67, 67, 67, 67, 67, 69, 69, 69, 69, 69, 69, 69, 71, 71, 71, 71, 71, 71, 71, 73, 73, 73, 73, 73, 73, 73, 73, 75, 75, 75, 75, 75, 75, 75, 75, 77, 77, 77, 77, 77, 77, 77, 77, 79, 79, 79, 79, 79, 79, 79, 79, 81, 81, 81, 81, 81, 81, 81, 81, 83, 83, 83, 83, 83, 83, 83, 83, 85, 85, 85, 85, 85, 85, 85, 85, 87, 87, 87, 87, 87, 87, 87, 87, 87, 89, 89, 89, 89, 89, 89, 89, 89, 89, 91, 91, 91, 91, 91, 91, 91, 91, 91, 93, 93, 93, 93, 93, 93, 93, 93, 93, 95, 95, 95, 95, 95, 95, 95, 95, 95, 97, 97, 97, 97, 97, 97, 97, 97, 97, 99, 99, 99, 99, 99, 99, 99, 99, 99, 101, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 103, 103, 103, 103, 103, 103, 103, 105, 105, 105, 105, 105, 105, 105, 105, 107, 107, 107, 107, 107, 107, 109, 109, 109, 109, 109, 111, 111, 111, 111, 111, 111, 111, 113, 113, 113, 113, 115, 115, 115, 115]
#jumps = [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 7, 8, 9, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 8, 9, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 6, 7, 8, 9, 6, 7, 8, 9, 6, 7, 8, 9]
#    
#sieve_dims = [46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 60, 60, 60, 60, 60, 60, 60, 60, 62, 62, 62, 62, 62, 62, 62, 62, 64, 64, 64, 64, 64, 64, 64, 64, 66, 66, 66, 66, 66, 66, 66, 66, 68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, 69, 69, 71, 71, 71, 71, 71, 71, 71, 71, 71, 73, 73, 73, 73, 73, 73, 73, 73, 73, 75, 75, 75, 75, 75, 75, 75, 75, 75, 77, 77, 77, 77, 77, 77, 77, 77, 77, 79, 79, 79, 79, 79, 79, 79, 79, 79, 81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 82, 82, 82, 84, 84, 84, 84, 84, 84, 84, 84, 84, 86, 86, 86, 86, 86, 86, 86, 86, 88, 88, 88, 88, 88, 88, 90, 90, 90, 90, 90, 92, 92, 92, 92, 92, 92, 92, 94, 94, 94, 94, 95, 95, 95, 95]
#    
#costs = [437.54, 232.42, 157.29, 124.5, 101.61, 83.94, 447.2, 232.85, 159.66, 124.26, 101.48, 87.23, 465.87, 241.76, 167.44, 129.08, 103.89, 92.42, 467.78, 243.86, 167.56, 128.22, 106.09, 93.87, 482.34, 253.68, 172.99, 136.7, 110.93, 96.19, 492.44, 253.33, 180.05, 136.07, 113.67, 94.79, 506.83, 262.6, 180.15, 138.95, 113.56, 97.36, 527.27, 271.86, 186.07, 141.43, 117.58, 101.37, 530.93, 280.5, 195.39, 151.37, 121.32, 106.16, 548.03, 283.98, 196.14, 154.26, 127.13, 111.19, 95.51, 567.17, 295.87, 202.96, 157.48, 130.34, 113.31, 93.89, 591.71, 303.78, 213.83, 156.12, 132.3, 112.33, 95.77, 86.92, 604.53, 317.07, 215.51, 169.63, 133.11, 116.65, 108.66, 90.66, 620.93, 326.4, 226.01, 180.32, 146.43, 127.42, 110.9, 99.46, 674.55, 348.94, 241.01, 189.86, 150.02, 131.61, 122.32, 108.82, 722.83, 391.04, 267.83, 211.64, 171.65, 146.69, 135.11, 115.29, 1030.59, 546.49, 375.72, 296.59, 233.44, 209.79, 193.71, 162.99, 1133.09, 598.68, 400.47, 317.9, 261.85, 221.93, 205.76, 181.61, 1214.91, 638.48, 431.21, 346.09, 280.06, 237.54, 222.46, 197.49, 175.47, 1321.3, 685.68, 476.5, 370.04, 295.81, 258.02, 237.45, 200.81, 193.6, 271.97, 243.95, 222.69, 211.72, 1389.99, 730.94, 498.46, 398.86, 322.26, 1468.69, 758.17, 525.14, 417.88, 327.43, 289.12, 266.86, 237.17, 228.74, 1554.68, 814.18, 571.95, 448.6, 353.03, 321.05, 296.44, 258.69, 248.59, 1634.07, 878.93, 606.06, 474.95, 389.23, 333.13, 325.01, 267.95, 255.91, 1807.1, 942.15, 645.15, 524.01, 414.5, 374.53, 339.58, 298.07, 283.29, 1914.8, 1028.76, 722.87, 548.38, 469.42, 419.97, 362.69, 317.68, 301.44, 2158.23, 1125.46, 791.87, 613.97, 515.19, 456.47, 409.92, 362.04, 341.74, 2464.25, 1305.35, 903.28, 702.42, 589.2, 545.66, 405.99, 385.26, 2877.41, 1529.64, 1080.1, 840.88, 740.96, 634.38, 3521.55, 1861.95, 1315.44, 1018.26, 857.54, 4338.66, 2304.19, 1609.93, 911.13, 819.06, 714.21, 684.63, 1142.04, 966.92, 856.21, 840.46, 1228.53, 1080.08, 939.46, 896.63]
#
#thre1 = 13
#thre2 = 15
#thre3 = 26
#draw_practical_pnjbkz_cost(blocksizes, jumps, sieve_dims, costs, 'PnjBKZCost-dd-ds.png',thre1, thre2, thre3)
