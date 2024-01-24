target=libest.so
cc=g++
# srcs=enumbs.cpp
# includes=-I./include/
OBJS= test.cpp utils.cpp cost.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp  enumbs.cpp bssa.cpp

CXXFLAGS=-fPIC -Ofast -march=native -ftree-vectorize -funroll-loops -std=c++11 -pthread -Wall -Wextra -O3 -march=native -Wp,-U_FORTIFY_SOURCE $(EXTRAFLAGS)  #-Wp,-D_FORTIFY_SOURCE=1
LDFLAGS=-shared -pthread -L$(realpath .) -Wl,-rpath=$(realpath .) #
LIBADD=-lcurl  -lalglib  -lmpfr -lgmp -lpthread  -lfplll #
#KERNELHEADERS = utils.h bkz_with_jump_simulator.h d_svp_prediction.h  cost.h  enumbs.h



all: $(target) 
$(target): $(OBJS) #$(KERNELHEADERS) 
	$(cc) $(CXXFLAGS)  $(OBJS) -o $@ $(LDFLAGS) $(LIBADD)
	
# %.o:%.cpp
# 	$(cc) -c -pthread $(INCLUDES) $(CFLAGS) $^ -o $@


	
# for python setup.py to obtain final values of variables
printvariables:
	$(foreach v,$(.VARIABLES),$(info $(v) = $($(v))))


clean:
	-rm  *.so