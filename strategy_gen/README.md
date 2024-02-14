# Guidance

This project implements the EnumBS algorithm, and the security estimation of NIST-round3 schemes and lwe challenges(https://www.latticechallenge.org/lwe_challenge/challenge.php) in EnumBS and BSSA. 

All of the estimators are implemented according to the  (nist)-bkz estimator implemented in leaky-lwe-estimator regarding the norm of shortest vector (e,1) for LWE instance in primal attack as a chi-squared distribution, and computing a success probability in it. Let `worst_case = True` and `cost_model = 1`, then it is the form for computing success probability, we use this estimator to estimate the NIST schemes, corresponding to the values of Two-step column shown in Table 7 in https://eprint.iacr.org/2022/1343.pdf. Set `worst_case = False` and `cost_model = 2`, then it is an estimation used for estimating the solvability TU LWE challenge and it returns a blocksize strategy for solving TU LWE challenge by EnumBS or BSSA, corresponding to the Table 3 in https://eprint.iacr.org/2022/1343.pdf.


### Dependencies

#### Required

* GNU MP 4.2.0 or higher http://gmplib.org/ or MPIR 1.0.0 or higher [http://mpir.org](http://mpir.org/)
* MPFR 2.3.0 or higher, COMPLETE INSTALLATION http://www.mpfr.org/
* autotools 2.61 or higher
* g++ 4.9.3 or higher
* boost 1.66 or higher (but lower than a boost using c++14) : https://www.boost.org/

In unbuntu, you can install the above dependences by implementing the following code directly:

```bash
./install-dependency.sh
```

Or install the dependencies on yourserlf: 

1. C++11

2. mpfr: https://www.mpfr.org/

```bash
apt-get update
apt-get install libmpfr-dev #ubuntu
```

3. gmp: https://gmplib.org/

```bash
apt-get install m4 #should install m4 first
apt-get install libgmp-dev
```

4. boost: 

```bash
#apt-get install libboost-all-dev

#apt-get install libboost1.67-all-dev 
wget https://boostorg.jfrog.io/artifactory/main/release/1.66.0/source/boost_1_66_0.tar.gz
tar -xzvf boost_1_66_0.tar.gz
cd boost_1_66_0
sudo ./bootstrap.sh
sudo ./b2 install

sudo gedit /etc/profile
export CPLUS_INCLUDE_PATH=/usr/local/include/boost:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

5. autoconf

```
apt-get install autoconf
```


6. download and install fplll library(https://github.com/fplll/fplll) in the folder `cpp`: 

```bash
git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
make install
```


### Organization of the code
There are 3 code folders in `cpp`.  `framework` contains core code for EnumBS(`enumbs.cpp`), BSSA(`bssa.cpp`), pnj-BKZ simulator(`bkz_with_jump_simulator.cpp`) and cost model(`cost.cpp`) implementations. If developers modify the codes in it, please run `./rebuild.sh` in the main direction. 

There are two files in the folder `lwe-est`: 
- `lwechal-bssa-est.cpp` is an executable file for generating blocksize strategy by BSSA algorithm and give a cost estimation for lwe instances provided in TU LWE challenge(https://www.latticechallenge.org/lwe_challenge/challenge.php). One can run it by running code 

- `lwechal-enumbs-parallel-est.cpp` is an executable file for generating blocksize strategy by EnumBS algorithm and give a cost estimation for lwe instances provided in TU LWE challenge(https://www.latticechallenge.org/lwe_challenge/challenge.php). 

One can test the blocksize strategy generation method in EnumBS and BSSA by running code 
```
./lwechal-est.sh
```
in the main directory. It will return the blocksize strategy for some of LWE challenges in the folder `lwechal-est-result`.


We can use the strategy generated in `lwechal-est-result` to solve the corresponding LWE problem by the progressive pnj-BKZ solver(https://github.com/Summwer/pro-pnj-bkz.git).


The code file `NIST-round3-est-gate.cpp` in `NIST-round3-est` is used for estimating all the LWE-based NIST schemes(Kyber, Dilithium and Frodo) by EnumBS estimator. One can test it by running the code 
```
./implement_all_NIST_schemes.sh
```
in the main directory. It will return the blocksize strategy for some of LWE challenges in the folder `lwechal-est-result`, corresponding to the Table 3 in https://eprint.iacr.org/2022/1343.


To simulate a LWE instance in detail and return the simulated basis quality(we regard "slope" as a basis quality representation) and simulated cost(in practical cost model with threads = 32 and gpus = 2 * 3039Ti GPUS (sec)) after each pnj-BKZ$(\beta,J, {\rm tours})$ reduction, please refer to the file `strategy_simulation.cpp` and one can run it by the following command
```
./strategy_simulation.sh
```
The result will store in the file `lwechal-est-result/simulation.log`, it obtains the 'Simulation' column for LWE challenge ($40,0.030$) in Table 5 in https://eprint.iacr.org/2022/1343.


### clean the output *.so file of enumbs

```
cd framework
make clean
cd ..
```


### compile strategy generation method in python

```
cd ..
source ./activate #in main directory
cd strategy_gen
python setup.py build_ext --inplace

```

If one need change the paramter setting in struct Params in utils.h, one should firstly delete the `build` folder and `libest.so` in the folder `framework`, then recompile the library.





