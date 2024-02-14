from random import randint
from sympy import nextprime
import numpy as np
from fpylll import IntegerMatrix


def store_lwe_instance(n,alpha,m,q,A,b):
    alpha_ = int(alpha*1000)
    filename = 'lwe_instances/%03d-%03d.txt' % (n, alpha_)
    fn = open(filename, "w")
    fn.write(str(n)+'\n')
    fn.write(str(m)+'\n')
    fn.write(str(q)+'\n')
    fn.write(str(alpha)+'\n')
    fn.write(str(b)+'\n')
    fn.write(str(A)+'\n')
    # fn.write('[')
    # for i in range(A.nrows):
    #     fn.write('[')
    #     for j in range(A.ncols):
    #         fn.write(str(g6k.M.B[i][j]))
    #         if j<g6k.M.B.ncols-1:
    #             fn.write(' ')
    #     if i < g6k.M.B.nrows-1:
    #         fn.write(']\n')
    # fn.write(']]')
    fn.close()

def gen_LWE_instance(n,alpha):
    m = n**2
    q = nextprime(m)
    sigma = alpha * q
    print("Generate LWE instance: n = %d, m = %d, alpha = %f, q = %d,  sigma = %f" %(n,m,alpha,q, sigma))


    A = IntegerMatrix.from_matrix([[np.random.randint(0,q-1) for _ in range(n)] for _ in range(m)])

    s = IntegerMatrix.from_matrix([[np.random.randint(0,q-1) for _ in range(n)]])

    e = IntegerMatrix.from_matrix([[round(_) for _ in np.random.normal(loc = 0, scale = sigma, size = m)]])

    print("secret vector s =", s)

    print("noise vector e = ", e)

    # print("random matrix A = ", A)

    b = (A*s.transpose()).transpose()  #+ e.transpose()
    b = IntegerMatrix.from_matrix([[(b[0][i]+e[0][i])%q for i in range(m)]])

    # print("b = ", b)
    
    #store b and A into the lwechallenge folder.
    
    

    store_lwe_instance(n,alpha,m,q,A,b)
    return A,b


# n = 60
# alpha = 0.006
# gen_LWE_instance(n,alpha)

