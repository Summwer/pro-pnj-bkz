from g6k.utils.util import gen_LWE_instance
from g6k.utils.lwe_estimation import gsa_params, primal_lattice_basis
from g6k.siever import Siever
import numpy as np
from fpylll import IntegerMatrix
from fpylll.util import gaussian_heuristic
from scipy.special import chdtr 
from math import sqrt
import matplotlib.pyplot as plt

def chi_square_estimate(rr,sigma, succ_prob =  0.999):
    d = len(rr)
    for beta in range(30,d):
        GH = gaussian_heuristic(rr[d-beta:])
        length=(GH/(sigma**2))
        if(chdtr(beta, length) >= succ_prob):
            return beta


def expected_value_estimate(rr,sigma):
    d = len(rr)
    for beta in range(30,d):
        if(gaussian_heuristic(rr[d-beta:]) >= (sigma**2 * (d-1) + 1) * beta/(1.*d)):
            return beta

  
def compute_projected_norm(target_vector,g6k,dsvp):      
    yl = g6k.M.from_canonical(target_vector)
    d = len(target_vector)
    return np.linalg.norm(np.array(g6k.M.to_canonical([0]*(d-dsvp) + list(yl)[d-dsvp:])))**2


#True(False): the condition is (not) satisfied.
def determine_correctness_condition(target_vector, g6k, dsvp,rr):
    projected_norm1 = compute_projected_norm(target_vector,g6k,dsvp)
    d = len(rr)
    if(projected_norm1 <= gaussian_heuristic(rr[d-dsvp:])):
        return 1
    else:
        return 0

def draw_length_figure(n,alpha, sigma, d, projected_target_lengths_for_LWE):
    plt.rc('font', size=24)
    fig, ax = plt.subplots(figsize=(9,7),dpi=600)
    
    dims = [_ for _ in range(d,50,-1)]
    for i in range(len(projected_target_lengths_for_LWE)-1):
        ax.scatter(dims,projected_target_lengths_for_LWE[i], marker=".", zorder = 2, color= "black") #linewidth =2,
    ax.scatter(dims,projected_target_lengths_for_LWE[-1], marker=".", zorder = 2, color= "black", label = r"actual projected target length")
        
    ax.scatter(dims,[sigma*sqrt(d-i) for i in range(d-50)], marker=".", zorder = 3, color= "blue", label = r"$\sigma\sqrt{d_{\rm svp}}$")
    # ax.plot(ns,enumbs_cost, linewidth =2,markersize=8, marker=marker_type, zorder = 5, color = 'red', label = r"ProPnjBKZ(EnumBS)")
    # ax.plot(ns,bssa_cost,linewidth =2,markersize=8, marker=marker_type, zorder = 3, color = 'orange', label = r"ProPnjBKZ(BSSA)" )
    
    
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.legend(fontsize = 20, loc = "upper left")#bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    plt.xlabel(r'$d_{\rm svp}$')# with $\alpha$ = %.3f' %alpha, fontsize = 24)
    plt.ylabel(r'length',fontsize = 24)
#    ax.autoscale(tight=False)
    # plt.ylim(0,40000)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    
    plt.grid(True)
    fig.savefig(r'projected-target-length-%d-%.3f.png' %(n,alpha),bbox_inches='tight')
    
    
    
def draw_fail_prob_figure(ns, alpha, fail_probs1, fail_probs2):
    plt.rc('font', size=24)
    fig, ax = plt.subplots(figsize=(9,7),dpi=600)
        
    ax.plot(ns,fail_probs1,  linewidth =2,markersize=8,marker=".", zorder = 2, color= "black", label = "default G6K")
    ax.plot(ns,fail_probs2,  linewidth =2,markersize=8,marker=".", zorder = 2, color= "red", label = "our work")
    
    
    # plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.legend(fontsize = 20, loc = "upper left")#bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    plt.xlabel(r'$n$')# with $\alpha$ = %.3f' %alpha, fontsize = 24)
    plt.ylabel(r"$Pr(\| \pi_{d-d_{\rm svp}}({\bf{t}}) \| > {\rm GH}({\bf{B}}_{\pi[d-d_{\rm svp}]})$)",fontsize = 24)
#    ax.autoscale(tight=False)
    plt.ylim(-0.1,1.1)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    
    plt.grid(True)
    fig.savefig(r'fail-prob-%.3f.png' %(alpha),bbox_inches='tight')
    


def determine_satisfied_projected_target_norm(n,alpha, m = None):
    A,c,e,s,q = gen_LWE_instance(n,alpha,store_file = False)

    if m is None:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True, samples = A.nrows)
            (b, _, m) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    else:
        try:
            min_cost_param = gsa_params(n=A.ncols, alpha=alpha, q=q,
                                        decouple=True)
            (b, _, _) = min_cost_param
        except TypeError:
            raise TypeError("No winning parameters.")
    # print("Chose %d samples." % (m))
    # print()
    
    sigma= alpha*q

    B = primal_lattice_basis(A, c, q, m=m) #debug
    
    g6k = Siever(B)
    d = g6k.full_n
    target_vector = list(e[:m])+[1]
    rr = [g6k.M.get_r(i,i) for i in range(d)]

    
    # projectd_target_norms =[] #sqaure_norm of projected target_vecotor t[i:]
    
    # for i in range(d-50):
    #     projectd_target_norms.append(np.linalg.norm(np.array(g6k.M.to_canonical([0]*i + list(yl)[i:]))))
    
  
    dsvp1 = expected_value_estimate(rr,sigma)
    dsvp2 = chi_square_estimate(rr,sigma)

    return determine_correctness_condition(target_vector, g6k, dsvp1,rr), determine_correctness_condition(target_vector, g6k, dsvp2,rr)
    
    
    
    
    # print(projectd_target_norms)
    # return projectd_target_norms,q,d
    # print([sigma*sqrt(d-i) for i in range(d - 50)])
    
    
    
#compute failure probability of the estimate (1) GH(rr[d-dsvp:])<=sigma*sqrt(dsvp)
#                                            (2) ChiProb(x<=GH(rr[d-dsvp:])) >= 0.999
def failure_probability_for_two_estimates(n,alpha, m = None, samples = 100):
    total_succ_amount1 = 0
    total_succ_amount2 = 0
    for i in range(samples):
        print("\r%d/%d"%(i+1,samples),end="")
        succ1, succ2 = determine_satisfied_projected_target_norm(n,alpha)
        total_succ_amount1+=succ1
        total_succ_amount2+=succ2
    
    print()
    # print(total_succ_amount1,total_succ_amount2)
    total_fail_probability1 = 1. - total_succ_amount1/samples
    total_fail_probability2 = 1. -  total_succ_amount2/ samples
    # draw_prob_figure(n,alpha,sigma, d, projected_target_lengths_for_LWE)
    return total_fail_probability1, total_fail_probability2

def est_prob_comparison(ns, alpha):
    fail_probs1 = []
    fail_probs2 = []
    for n in ns:
        print("n = %d, alpha = %.3f" %(n,alpha))
        fail_prob1, fail_prob2 = failure_probability_for_two_estimates(n,alpha)
        fail_probs1.append(fail_prob1)
        fail_probs2.append(fail_prob2)
        # print(n,alpha, "prob(estimate1) = ", )

    draw_fail_prob_figure(ns, alpha, fail_probs1, fail_probs2)

        
ns = list(range(48,53))
alpha = 0.015   
est_prob_comparison(ns, alpha)

ns = list(range(54,62))
alpha = 0.010
est_prob_comparison(ns, alpha)