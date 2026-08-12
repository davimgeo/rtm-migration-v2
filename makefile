CC       = gcc
CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma

INCLUDE  = -Iinclude -Isrc

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE  = $(shell python3 -c 'import numpy; print(numpy.get_include())')

LIBS    = -lm
PYTHON_LIBS = $(shell python3-config --embed --ldflags)

SRC = $(shell find . -name "*.c")

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(PYTHON_INCLUDE) \
	-I$(NUMPY_INCLUDE) $(SRC) $(LIBS) $(PYTHON_LIBS) -o run.out

	./run.out

run: all
	./run.out

clean:
	rm -f run.out
