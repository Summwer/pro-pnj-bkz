cd framework
export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/local/bin/:/usr/include/:/usr/lib"
# g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++11 -o ../libest.so utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp  bssa.cpp cost.cpp attack.cpp est.cpp load_lwechal.cpp lwe_estimation.cpp -L. -pthread -lgmp -lfplll -lmpfr  -fPIC -shared -lalglib -lcurl   #-Wall -Wextra


g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++11 -o ../libest.so utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp  bssa.cpp cost.cpp attack.cpp est.cpp load_lwechal.cpp lwe_estimation.cpp -shared -pthread -L/home/cryptothesis/summer/pro-pnj-bkz/strategy_gen/framework -Wl,-rpath=/home/cryptothesis/summer/pro-pnj-bkz/strategy_gen/framework -lpthread -lalglib -lcurl -lfplll -lmpfr -lgmp


cd ..
  
