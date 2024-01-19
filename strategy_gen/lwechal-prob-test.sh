export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir "lwechal-prob-test"

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
# argv[8]: list decoding. 1: agps20 2:matzov22
g++ -O3 -funroll-loops -o lwechal-enumbs-est ./lwe-est/lwechal-enumbs-parallel-est.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-enumbs-est 100 5 1 300 2 3 50 1 | tee "lwechal-prob-test"/"EnumBS-Strategy(cumprob+prob)+list decoding in [AGPS20](gate min).log"

./lwechal-enumbs-est 100 5 1 300 2 3 50 2 | tee "lwechal-prob-test"/"EnumBS-Strategy(cumprob+prob)+list decoding in [MATZOV22](gate min).log" 
