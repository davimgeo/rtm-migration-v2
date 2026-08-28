CC       = gcc
NVCC     = nvcc

CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma

NVFLAGS  = -O3 -g -std=c++17
NVFLAGS += -Xcompiler=-march=native
NVFLAGS += -Xcompiler=-fopenmp

INCLUDE  = -Iinclude -Isrc

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE  = $(shell python3 -c 'import numpy; print(numpy.get_include())')
PYTHON_LIBS    = $(shell python3-config --embed --ldflags)

LIBS     = -lm -fopenmp

TARGET   = run.out

C_SRC    = $(shell find src -name "*.c") config/config.c
CU_SRC   = $(shell find src -name "*.cu")

C_OBJ    = $(patsubst %.c,build/%.o,$(C_SRC))
CU_OBJ   = $(patsubst %.cu,build/%.o,$(CU_SRC))

OBJ      = $(C_OBJ) $(CU_OBJ)
MAIN_OBJ = build/main.o

all: $(TARGET)
	@./$(TARGET)

$(TARGET): $(OBJ) $(MAIN_OBJ)
	@$(NVCC) $(OBJ) $(MAIN_OBJ) $(LIBS) $(PYTHON_LIBS) -o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(INCLUDE) \
		$(PYTHON_INCLUDE) \
		-I$(NUMPY_INCLUDE) \
		-c $< -o $@

build/%.o: %.cu
	@mkdir -p $(dir $@)
	@$(NVCC) $(NVFLAGS) $(INCLUDE) \
		$(PYTHON_INCLUDE) \
		-I$(NUMPY_INCLUDE) \
		-c $< -o $@

clean:
	rm -rf build $(TARGET)
