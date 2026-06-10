CC       = gcc
CXX      = g++

CFLAGS   = -std=c99 -O3 -march=native -mavx2 -mfma -g
CXXFLAGS = -std=c++11 -O3 -march=native -mavx2 -mfma -g

LDFLAGS  = -fopenmp
OMPFLAGS = -fopenmp

CPP_SRC  = main.cpp src/geometry.cpp src/model.cpp \
           src/wavelet.cpp src/seismogram.cpp \
           src/propagation.cpp

C_SRC    = src/utils.c

INIT ?= tests/parallel.c

CPP_OBJ  = $(CPP_SRC:.cpp=.o)
C_OBJ    = $(C_SRC:.c=.o)
INIT_OBJ = $(INIT:.c=.o)

run: $(CPP_OBJ) $(C_OBJ) $(INIT_OBJ)
	$(CXX) $^ $(LDFLAGS) -o run.out
	./run.out

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c $< -o $@

parallel: 
	$(MAKE) run INIT=tests/parallel.c

marmousi: 
	$(MAKE) run INIT=tests/marmousi.c

clean: 
	rm -f src/*.o tests/*.o run.out
