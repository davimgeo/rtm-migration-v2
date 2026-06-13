CC      = gcc
CFLAGS  = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma
INCLUDE = -I/home/lum/rtm-migration-v2/include

C_SRC := $(shell find . -name "*.c")

run:
	$(CC) $(CFLAGS) $(INCLUDE) $(C_SRC) -lm -o run.out
	./run.out
	$(MAKE) clean

clean:
	rm -f run.out
