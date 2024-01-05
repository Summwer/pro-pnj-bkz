target=libest.so
cc=g++
# srcs=enumbs.cpp
includes=-I./include/
OBJS=utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp cost.cpp 

CXXFLAGS+=-fPIC -Ofast -march=native -ftree-vectorize -funroll-loops -std=c++14 -pthread -Wall -Wextra -O3 -march=native -Wp,-U_FORTIFY_SOURCE $(EXTRAFLAGS)
LDFLAGS+=-shared -pthread -L$(realpath .) -Wl,-rpath=$(realpath .)
LIBADD+=-lpthread -lalglib -lcurl -lfplll -lmpfr -lgmp
KERNELHEADERS = utils.h bkz_with_jump_simulator.h d_svp_prediction.h enumbs.h cost.h 


all: $(target) 
$(target): $(OBJS) $(KERNELHEADERS) 
	$(cc) $(CXXFLAGS)  $(OBJS) -o $@ $(LDFLAGS) $(LIBADD)
	
%.o:%.cpp
	$(cc) -c -pthread $(INCLUDES) $(CFLAGS) $^ -o $@

clean:
	-rm *.o *.so ../*.so
	
# for python setup.py to obtain final values of variables
printvariables:
	$(foreach v,$(.VARIABLES),$(info $(v) = $($(v))))

