CC      = gcc
CFLAGS  = -std=gnu99 -g -lm

# Model
C_SRC := $(wildcard *.c) $(shell find ../src/model/ -type f -name '*.c') \
         ../src/io.c ../src/utils.c
# All
#C_SRC := $(wildcard *.c) $(shell find ../src -type f -name '*.c')

INCLUDE = -I/home/lum/coding/rtm-migration-v2/rtm-migration-c/include

C_OBJ   = $(C_SRC:.c=.o)

run: $(C_OBJ)
	$(CC) $^ -o run.out
	./run.out
	$(MAKE) clean

%.o: %.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(C_OBJ) run.out
