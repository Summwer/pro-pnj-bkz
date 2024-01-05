export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
# g++ -O3 -funroll-loops -o strategy_simulation strategy_simulation.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest
# ./strategy_simulation | tee lwechal-est-result/simulation.log

mkdir "lwechal-prob-test"

g++ -O3 -funroll-loops -o strategy_simulation_for_cum_prob strategy_simulation_for_cum_prob.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest
./strategy_simulation_for_cum_prob | tee lwechal-prob-test/"EnumBS(cumprob+prob)-simulation.log"