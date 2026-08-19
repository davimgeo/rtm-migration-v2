CC       = gcc
CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma
INCLUDE  = -Iinclude -Isrc

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE  = $(shell python3 -c 'import numpy; print(numpy.get_include())')
PYTHON_LIBS    = $(shell python3-config --embed --ldflags)

LIBS     = -lm

TARGET   = run.out

SRC      = $(shell find src -name "*.c") config/config.c
OBJ      = $(SRC:%.c=build/%.o)
MAIN_OBJ = build/main.o

all: $(TARGET)
	./$(TARGET)

$(TARGET): $(OBJ) $(MAIN_OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(MAIN_OBJ) $(LIBS) $(PYTHON_LIBS) -o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDE) $(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -c $< -o $@

clean:
	rm -rf build $(TARGET)
