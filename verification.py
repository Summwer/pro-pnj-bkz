from hnp_gpu import load_matrix,Euclidean_Norm,Inf_Norm,recover_alpha,verification,get_k_list
from readfile import challenge



def try_to_find_solution(db,Coeff,w,K_list,q,s,EquaNum,bitsize):
    vec = db[0]
    print(float(Euclidean_Norm(vec,q)))
    print(float(Inf_Norm(vec,q)))
    for index in range(EquaNum):
        if Coeff[index] %2  == 1:
            if(vec[-1] == w or vec[-1] == -w):
                # print(check_len_bin(vec,bitsize,s),end = ' ')
                # print(vec)
                alpha = recover_alpha((w//vec[-1])*vec[-2]+w+K_list[index],index,Coeff[index],q)
                if(verification(Coeff,K_list,q,s,EquaNum,bitsize,alpha)):
                    print("Solution = %d" %alpha)
                    return True
    return False


A = load_matrix("hnp-009-midmat-127.txt")
# A = load_matrix("hnp-008-after-pump2.txt")
q,s,EquaNum,Coeff,KnownNonce,bitsize,K_list,solution = challenge(9)

w = 1 << (bitsize - s  - 1)

for i in range(A.nrows):
    # print(float(Euclidean_Norm(A[i],q)))#, float(w))
    print(float(Inf_Norm(A[i],q)), float(w))
# print(float(Inf_Norm(A[0],q)))
# print(float(w))


try_to_find_solution(A,Coeff,w,K_list,q,s,EquaNum,bitsize)


# k_list = get_k_list(Coeff,K_list,q,s,EquaNum,bitsize,solution)
# print(k_list)


