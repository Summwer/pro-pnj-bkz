'''python version of the reference kernel--> simhash.inl'''

from random import shuffle 


def Euclidean_Norm(vector):
    Euclidean_norm = sum([abs(_)**2 for _ in vector])
    return Euclidean_norm

class SimHashes(object):
    def __init__(self):
        self.XPC_BIT_LEN = 256
        self.XPCes_WORD_LEN = 4
    
    

    '''initiate a indices array'''
    def reset_compress_pos(self,d):

        # create random permutation of 0..n-1:
        permut = [i for i in range(d)]
        shuffle(permut)

        #read the indices reference file v
        file_name = "/home/cryptothesis/g6k/g6k/spherical_coding/sc_"+str(d)+"_256.def"
        f = open(file_name,'r')
        v_list = ' '.join(f.read().split('\n')).split(' ')
        v_list.remove('')
        v_list = [int(_) for _ in v_list]

        #initiate compress_pos for dimension > 30
        compress_pos = []
        for y in range(self.XPC_BIT_LEN):
            row = []
            for x in range(6):
                #pull out the element in v_list in order
                index = (x+6*y)%len(v_list) 
                v = v_list[index]
                row.append(permut[v])
            compress_pos.append(row)

        return compress_pos


    def compress(self,vector,compress_pos):
        c = [] #compressed vector
        for j in range(self.XPCes_WORD_LEN):
            c_tmp = int(0)
            a = 0
            for i in range(64):
                k = 64 * j + i
                a   = vector[compress_pos[k][0]]
                a  += vector[compress_pos[k][1]]
                a  += vector[compress_pos[k][2]]
                a  -= vector[compress_pos[k][3]]
                a  -= vector[compress_pos[k][4]]
                a  -= vector[compress_pos[k][5]]

                c_tmp = c_tmp << 1
                c_tmp |= int(a>0)
            #bin_c_tmp = bin(c_tmp)[2:]
            #c.append(('0'*(64-len(bin_c_tmp ))+bin_c_tmp ))
            c.append(c_tmp)
        return c

    def hamming_distance(self,c1,c2):       
        hamming_distance = 0
        for i in range(self.XPCes_WORD_LEN):
            hamming_distance += bin(c1[i]^c2[i])[2:].count('1')
        return hamming_distance
    
    


''' test
simhash = SimHashes()
compress_pos = simhash.reset_compress_pos(50)


vector1 = (-13, 1, -146, 277, -107, -180, 673, -311, -167, 47, 200, 395, 167, -25, -136, -392, 117, -165, 147, -515, 185, 637, 343, 8, 247, 44, -220, -146, 52, 135, -347, -369, -332, -102, 469, -285, 1, 167, 397, 84, -97, -138, -135, 218, 567, 141, 72, 21, 312, -41)

vector2 = (0, -124, -146, 277, -107, -180, 673, -311, -167, 47, 200, 395, 167, -25, -136, -392, 117, -165, 147, -515, 185, 637, 343, 8, 247, 44, -220, -146, 52, 135, -347, -369, -332, -102, 469, -285, 1, 167, 397, 84, -97, -138, -135, 218, 567, 141, 72, 21, 312, -41)

vector3 = (0, 0, -146, 0, -107, -180, 673, -311, -167, 47, 200, 395, 167, -25, -136, -392, 117, -165, 147, -515, 185, 637, 343, 8, 247, 44, -220, -146, 52, 135, -347, -369, -332, -102, 469, -285, 1, 167, 397, 84, -97, -138, -135, 218, 567, 141, 72, 21, 312, -41)

c1 = simhash.compress(vector1,compress_pos)
c2 = simhash.compress(vector2,compress_pos)
c3 = simhash.compress(vector3,compress_pos)



print(simhash.hamming_distance(c1,c2))
print(simhash.hamming_distance(c1,c3))
print(simhash.hamming_distance(c2,c3))

'''
