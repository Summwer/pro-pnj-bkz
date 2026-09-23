import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt,floor
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


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)


def independent_pump_cost_fun(beta):
    if(beta <= 86):
        k1 = 0.064883
        k2 = - 0.172902
    elif(beta <= 110):
        k1 = 0.196298
        k2 = - 11.518570
    elif(beta <= 131):
        k1 = 0.277686
        k2 =  -20.507189
    return 2**(k1*beta+k2)

def read_each_pump_cost_file(d,blocksize,jump):
    filename = "each_pump_cost_in_pnjbkz/pnjbkz-(%d,%d,%d)-1.txt" %(d,blocksize,jump)
    try:
        fn = open(filename,"r")
    except FileNotFoundError:
        return None, None
    data = fn.readlines()
    data[0] = (data[0].split(":"))[1]
    data[0] = data[0].replace("\n","")
    indices = eval(data[0])
    data[1] = (data[1].split(":"))[1]
    data[1] = data[1].replace("\n","")
    pump_costs = eval(data[1])
    return indices,pump_costs
    
        


def f_1(x, A, B):
    return A * x + B
    
    
def compute_R2(actual,predict):
    corr_matrix = np.corrcoef(actual, predict)
    corr = corr_matrix[0,1]  #相关系数
    R_sq = corr**2
    return R_sq
    
    
#extract cost data
def draw_practical_pump_cost_in_each_pnjBKZ(d, blocksize, jump, indices, pump_costs,simulated_pre_pnjbkz_cost = None, simulated_later_pnjbkz_cost=None, simulated_first_pump_cost = None, A=None,B=None, xf=12):
    f = dim4free_wrapper(default_dim4free_fun, blocksize)
#    print(f+1,pre_indice_range, d-blocksize+f-12+1, later_indice_range, d-blocksize+2*f-xf+1, len(indices))
    
    
    if(jump>1):
        full_indices,_= read_each_pump_cost_file(d,blocksize,1)
        index_for_pump_cost = []
        for i in range(len(full_indices)):
            if(full_indices[i] in indices):
                index_for_pump_cost.append(i)
           
    else:
        index_for_pump_cost = list(range(len(pump_costs)))
    
    pre_indice_range = f+xf
    later_indice_range = d-blocksize+f+1
    
    #pre_indices_range = f+1
    #len(indices)-later_indice_range = f
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,3),dpi=600)
        ax.scatter(index_for_pump_cost,pump_costs, marker=".", s=10)
        ax.plot([pre_indice_range, pre_indice_range],[min(pump_costs),max(pump_costs[1:])],label="pre_indice_range = %d"%pre_indice_range,color = "r",linestyle="--")
        ax.plot([later_indice_range, later_indice_range],[min(pump_costs),max(pump_costs[1:])],label="later_indice_range = %d"%later_indice_range,color = "k",linestyle="--")
#        ax.plot([0,len(pump_costs)-1],[independent_pump_cost/2., independent_pump_cost/2.],label="T_pump/2 = %.2f sec"%(independent_pump_cost/2.),color = "g",linestyle="--")
#        ax.plot([0,len(pump_costs)-1],[independent_pump_cost, independent_pump_cost],label="T_pump = %.2f sec"%(independent_pump_cost),color = "purple",linestyle="--")
        
#        ax.plot([0,len(pump_costs)-1],[sum(pump_costs)/len(pump_costs), sum(pump_costs)/len(pump_costs)],label="AvgPumpCost = %.2f sec"%(sum(pump_costs)/len(pump_costs)),color = "g",linestyle="--")
        midAvgPumpCost = sum(pump_costs[pre_indice_range:later_indice_range])/(later_indice_range-pre_indice_range)
        if(A is None and B is None):
            A, B = optimize.curve_fit(f_1, list(range(pre_indice_range, later_indice_range)), pump_costs[pre_indice_range:later_indice_range])[0]
#        print(A,B)
#        ax.plot([pre_indice_range,later_indice_range],[midAvgPumpCost , midAvgPumpCost ],label="midAvgPumpCost = %.2f sec"%(midAvgPumpCost),color = "b",linestyle="--")
        ax.plot([pre_indice_range,later_indice_range],[A*pre_indice_range+B, A*later_indice_range+B],label=r"$T_{\rm mid}$ = %.3f i + %.3f"%(A, B),color = "b",linestyle="--")
        
#        MinPumpCost = min(pump_costs)
#        ax.plot([0,len(pump_costs)-1],[MinPumpCost , MinPumpCost ],label="MinPumpCost = %.2f sec"%(MinPumpCost),color = "brown",linestyle="--")
        plt.legend(fontsize="8")
        ax.autoscale(tight=False)
        ax.set_title("d = %d, blocksize = %d, jump = %d, f = %d, exf = 12, sieve_dim = %d"%(d,blocksize,jump,f,blocksize-f), fontsize="8")
        ax.set(xlabel=r'index: i')
        ax.set(ylabel=r'cost (sec)')
        fig.savefig(r'each_pump_cost_in_pnjbkz/auxiliary_fig/EachPumpCost-dd-in-PnJBKZ(%d,%d,%d).png'%(d,blocksize,jump))
        plt.close()
        
#        simulated_pnjbkz_cost = midAvgPumpCost*(d-blocksize+1 + f )
        if(simulated_pre_pnjbkz_cost is None):
            assert(jump == 1)
            simulated_pre_pnjbkz_cost = sum(pump_costs[1:f+xf])
        if(simulated_later_pnjbkz_cost is None):
            assert(jump == 1)
            simulated_later_pnjbkz_cost = sum(pump_costs[d-blocksize+f+1:])
        if(simulated_first_pump_cost is None):
            assert(jump == 1)
            simulated_first_pump_cost = pump_costs[0]
        simulated_mid_cost = ((A*(f+xf) + B) + (A*(d-blocksize+f)+B))/2. *(d-blocksize-xf+1)
#        simulated_pnjbkz_cost = (simulated_pre_pnjbkz_cost + simulated_later_pnjbkz_cost + simulated_mid_cost)/jump

        simulated_pnjbkz_cost = simulated_first_pump_cost + simulated_pre_pnjbkz_cost * (ceil((f+xf)/jump)-1)/(f+xf-1) + simulated_later_pnjbkz_cost * (ceil((f+xf)/jump))/(f+xf) + simulated_mid_cost *  (ceil((d-blocksize-xf)/jump)+1)/(d - blocksize - xf+1)
        
        pnjbkz_cost = sum(pump_costs)
        print("PnJBKZ-(%d,%d,%d), PnJBKZ Cost = %.2f sec, simulated PnJBKZ Cost = %.2f sec" %(d, blocksize, jump, sum(pump_costs),simulated_pnjbkz_cost))
        
        return round(A,5),round(B,5),round(pnjbkz_cost,2),round(simulated_pnjbkz_cost,2),round(simulated_pre_pnjbkz_cost,2),round(simulated_later_pnjbkz_cost,2),round(pump_costs[0],2)
        



def draw_practical_pump_cost_in_each_pnjBKZ_with_different_d(fig, ax, d, blocksize, jump, indices, pump_costs, pre_line, simulated_pre_pnjbkz_cost = None, simulated_later_pnjbkz_cost=None, simulated_first_pump_cost = None, A=None,B=None, xf=12):
    f = dim4free_wrapper(default_dim4free_fun, blocksize)

    if(jump>1):
        full_indices,_= read_each_pump_cost_file(d,blocksize,1)
        index_for_pump_cost = []
        for i in range(len(full_indices)):
            if(full_indices[i] in indices):
                index_for_pump_cost.append(i)
           
    else:
        index_for_pump_cost = list(range(len(pump_costs)))
    
    pre_indice_range = f+xf
    later_indice_range = d-blocksize+f+1
    
    with plt.style.context(['seaborn-white','science']):
        ax.scatter(index_for_pump_cost,pump_costs, marker=".", s=10, label="$d$ = %d" %d)
        if(not pre_line):
            ax.plot([pre_indice_range, pre_indice_range],[min(pump_costs),max(pump_costs[1:])],label=r"f+xf = %d"%pre_indice_range,color = "r",linestyle="--")
            pre_line = True
        ax.plot([later_indice_range, later_indice_range],[min(pump_costs),max(pump_costs[1:])],label=r"d-$\beta$+f+1 = %d"%(later_indice_range),linestyle="--") #color = "k",

    return pre_line
            
def draw_sim_pnjbkz_cost_prediction(d, blocksizes,jump, simulated_pnjbkz_costs, pnjbkz_costs):
    if(len(blocksizes)==0):
        return
    with plt.style.context(['seaborn-white','science']):
        fig, ax = plt.subplots(figsize=(5,3),dpi=600)
        R2 = compute_R2(pnjbkz_costs,simulated_pnjbkz_costs)
        ax.plot(blocksizes,simulated_pnjbkz_costs,label="Simulated PnJBKZ cost for J=%d, $R^2$ = %.6f"%(jump,R2),color = "r",marker='.', markersize='6', zorder = 5, alpha = 0.5)
        ax.plot(blocksizes,pnjbkz_costs,label="Actual PnJBKZ cost for J=%d" %jump,color = "b",marker='.', markersize='6', zorder = 4, alpha = 0.5)
        plt.legend(fontsize="10")
        ax.autoscale(tight=False)
        ax.set(xlabel=r'$\beta$')
        ax.set(ylabel=r'cost (sec)')
        fig.savefig("each_pump_cost_in_pnjbkz/PnJBKZ_Cost_Prediction/PnJBKZ_cost_prediction_in_svp-%d-with-jump-%d.png" %(d,jump))
        plt.close()
        
def draw_and_simulate_PnJBKZ_costs(d,jump,xf=12,blocksize_line=None,simulated_pre_pnjbkz_costs=None,simulated_later_pnjbkz_costs=None,simulated_first_pump_costs=None):
    compute_blocksize_line = False #Determine whether to compute the paramters for blocksize simulated cost
    if(blocksize_line is None):
        blocksize_line = {}
        compute_blocksize_line = True
        simulated_pre_pnjbkz_costs = {}
        simulated_later_pnjbkz_costs = {}
        simulated_first_pump_costs = {}
    pnjbkz_costs = []
    simulated_pnjbkz_costs = []
    blocksizes = []
    
    for blocksize in range(51,120,2):
        sieve_dim = blocksize-dim4free_wrapper(default_dim4free_fun, blocksize)
#        independent_pump_cost = independent_pump_cost_fun(sieve_dim)
        indices,pump_costs= read_each_pump_cost_file(d,blocksize,jump)
        if(indices is not None and  pump_costs is not None):
            if(compute_blocksize_line):
                #if there is no cost model, then generate one.
                A, B, pnjbkz_cost,simulated_pnjbkz_cost, simulated_pre_pnjbkz_cost, simulated_later_pnjbkz_cost, simulated_first_pump_cost= draw_practical_pump_cost_in_each_pnjBKZ(d, blocksize, jump, indices,pump_costs)
                blocksize_line[blocksize] = (A,B)
                simulated_pre_pnjbkz_costs[blocksize] = simulated_pre_pnjbkz_cost
                simulated_later_pnjbkz_costs[blocksize] = simulated_later_pnjbkz_cost
                simulated_first_pump_costs[blocksize] = simulated_first_pump_cost
            else:
                _,_, pnjbkz_cost,simulated_pnjbkz_cost,_,_,_ = draw_practical_pump_cost_in_each_pnjBKZ(d, blocksize,jump, indices,pump_costs,simulated_pre_pnjbkz_cost=simulated_pre_pnjbkz_costs[blocksize], simulated_later_pnjbkz_cost= simulated_later_pnjbkz_costs[blocksize],simulated_first_pump_cost = simulated_first_pump_costs[blocksize], A=blocksize_line[blocksize][0],B=blocksize_line[blocksize][1])
            pnjbkz_costs.append(pnjbkz_cost)
            simulated_pnjbkz_costs.append(simulated_pnjbkz_cost)
            blocksizes.append(blocksize)
        
    draw_sim_pnjbkz_cost_prediction(d, blocksizes, jump, simulated_pnjbkz_costs, pnjbkz_costs)

    if(compute_blocksize_line):
        return blocksize_line,simulated_pre_pnjbkz_costs,simulated_later_pnjbkz_costs,simulated_first_pump_costs
    else:
        return


#blocksize_line: store the mid simulated cost parameter of PnJBKZ
#Give a simulated cost model for both blocksize_line,simulated_pre_pnjbkz_costs,simulated_later_pnjbkz_costs,simulated_first_pump_costs
for d in range(180,150,-10):
    for jump in range(1, dim4free_wrapper(default_dim4free_fun, 180)+1):
        if(jump == 1 and d == 180):
            blocksize_line,simulated_pre_pnjbkz_costs,simulated_later_pnjbkz_costs,simulated_first_pump_costs = draw_and_simulate_PnJBKZ_costs(d,jump)
            
            f = open("pnjbkz_cost_model_data.txt", 'w')
            f.write("{0: <25} {1: <25} {2: <25} {3: <25} {4: <25} {5: <25} {6: <25}\n".format("blocksize","sieve_dim", "A_of_AvgMidCost","B_of_AvgMidCost","SimPnJBKZCost[0]","SimPnJBKZCost[1:f+xf]", "SimPnJBKZCost[d-β-xf+1:d]"))
            for blocksize in range(51,120,2):
                sieve_dim = blocksize-dim4free_wrapper(default_dim4free_fun, blocksize)
                f.write("{0: <25} {1: <25} {2: <25} {3: <25} {4: <25} {5: <25} {6: <25}\n".format(blocksize,sieve_dim, blocksize_line[blocksize][0],blocksize_line[blocksize][1],simulated_first_pump_costs[blocksize], simulated_pre_pnjbkz_costs[blocksize], simulated_later_pnjbkz_costs[blocksize]))
            f.close()
        else:
            draw_and_simulate_PnJBKZ_costs(d,jump,blocksize_line=blocksize_line,simulated_pre_pnjbkz_costs=simulated_pre_pnjbkz_costs,simulated_later_pnjbkz_costs=simulated_later_pnjbkz_costs,simulated_first_pump_costs=simulated_first_pump_costs )



#each pump cost with different d
for blocksize in range(51,120,2):
    f = dim4free_wrapper(default_dim4free_fun, blocksize)
    for jump in range(1, dim4free_wrapper(default_dim4free_fun, 180)+1):
        print("draw pump cost for blocksize = %d, jump = %d" %(blocksize,jump))
        fig, ax = plt.subplots(figsize=(5,4),dpi=600)
        flag = False
        pre_line = False
        for d in range(180,150,-10):
            indices,pump_costs= read_each_pump_cost_file(d,blocksize,jump)
            if(indices is not None and  pump_costs is not None):
                pre_line = draw_practical_pump_cost_in_each_pnjBKZ_with_different_d(fig, ax, d, blocksize, jump, indices,pump_costs, pre_line)
                flag = True
        if(flag):
            plt.legend(fontsize="8",ncol=2)
            ax.autoscale(tight=False)
            ax.set_title("blocksize = %d, jump = %d, f = %d, exf = 12, sieve_dim = %d"%(blocksize,jump,f,blocksize-f), fontsize="8")
            ax.set(xlabel=r'index: i')
            ax.set(ylabel=r'cost (sec)')
            plt.tight_layout()
            fig.savefig(r'each_pump_cost_in_pnjbkz/comparison_d/EachPumpCost-dd-in-PnJBKZ(%d,%d).png'%(blocksize,jump)) #, interpolation="nearest"
            
        plt.close()
            
            
