
CC          = gcc
NVCC        = nvcc
AR          = ar

CFLAGS      = -std=gnu11 -O3 -g -march=native -fopenmp -mavx2 -mfma
CUDAFLAGS   = -Xcompiler -fopenmp --std=c++11 -O3

CFLAGS += -ffp-contract=off
CUDAFLAGS += --fmad=false

INCLUDE     = -Iinclude -Isrc -Itests

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE = $(shell python3 -c 'import numpy; print(numpy.get_include())')
PYTHON_LIBS = $(shell python3-config --embed --ldflags)
PYTHON_LIBDIR = $(shell python3 -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')

LIBS        = -lm -lfftw3f

TARGET      = run.out

C_SRC       = $(shell find src -name "*.c") tests/moving_rickers.c tests/euclidian.c config/config.c main.c
CU_SRC      = $(shell find src -name "*.cu")

C_OBJ       = $(patsubst %.c,build/%.o,$(C_SRC))
CU_OBJ      = $(patsubst %.cu,build/%.cu.o,$(CU_SRC))

.PHONY: run clean

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(C_OBJ) $(CU_OBJ)
	$(NVCC) -ccbin $(CC) \
	$(C_OBJ) $(CU_OBJ) \
	$(LIBS) \
	$(PYTHON_LIBS) \
	-Xcompiler=-fopenmp \
	-Xlinker=-rpath,$(PYTHON_LIBDIR) \
	-o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(INCLUDE) \
	$(PYTHON_INCLUDE) \
	-I$(NUMPY_INCLUDE) \
	-c $< -o $@

build/%.cu.o: %.cu
	@mkdir -p $(dir $@)
	@$(NVCC) $(CUDAFLAGS) $(INCLUDE) \
	$(PYTHON_INCLUDE) \
	-I$(NUMPY_INCLUDE) \
	-c $< -o $@

clean:
	rm -rf build $(TARGET)

