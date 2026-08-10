CC       = gcc
CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma
INCLUDE  = -Iinclude -Isrc
LIBS     = -lm

SRC = $(shell find . -name "*.c")

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC) $(LIBS) -o run.out

	./run.out

clean:
	rm -f run.out