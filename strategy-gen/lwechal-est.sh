export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir "lwechal-est-result"

# cd "lwechal-est-result"
# mkdir "enumbs(32+2gpus)-d4f-default-g6k"
# mkdir "enumbs(32+2gpus)-d4f-theo2"
# mkdir "enumbs(32+2gpus)-d4f-theo1"
# cd ..

# params input in main function in enumbs-est
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: threads number
# argv[6]: practical_pump_d4f
# argv[7]: start beta value.
# argv[8]: maximal RAM memory.
g++ -O3 -funroll-loops -o lwechal-enumbs-est ./lwe-est/lwechal-enumbs-parallel-est.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-enumbs-est 100 5 2 300 2 3 50 0 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=5.log" 

# ./lwechal-enumbs-est 100 5 2 300 2 3 50 1 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=5-limmaxRAM.log" 

# ./lwechal-enumbs-est 100 5 2 300 2 2 50 0 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-theo2-J=100-maxloop=5.log" 

# ./lwechal-enumbs-est 100 5 2 300 2 1 50 0 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-theo1-J=100-maxloop=5.log"

# ./lwechal-enumbs-est 100 5 1 300 2 1 50 0 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-theo1-J=100-maxloop=5-costmodel1.log"


# ./lwechal-enumbs-est 100 5 2 300 2 3 50 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=5-betastart-55.log" #{{80, 0.005}}

# ./lwechal-enumbs-est 8 5 2 300 2 3 50 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=8-maxloop=5.log" #{{80, 0.005}}

# params input in main function in bssa-est
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: enumbs_min_G: 0-- false, find minimal RAM strategy; 1--true, find minimal time cost strategy
# argv[6]: max_RAM
# argv[7]: practical_pump_d4f
# argv[8]: 1:tradional bssa; 2:improved bssa
# argv[9]: start beta value.
g++ -O3 -funroll-loops -o lwechal-bssa-est ./lwe-est/lwechal-bssa-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-bssa-est 100 5 2 300 1 0 3 1 50 | tee lwechal-est-result/"bssa(32+2gpus)-d4f-default-g6k=J=100-maxloop=5.log"

# ./lwechal-bssa-est 8 5 2 300 1 0 3 1 50 | tee lwechal-est-result/"bssa(32+2gpus)-d4f-default-g6k=J=8-maxloop=5.log"





