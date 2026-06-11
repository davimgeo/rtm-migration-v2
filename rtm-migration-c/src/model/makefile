CC       = gcc

CFLAGS   = -std=gnu99 -g

C_SRC    = ../utils.c main.c extent.c ../io.c create.c

INCLUDE =  -I/home/lum/coding/rtm-migration-v2/rtm-migration-c/include

C_OBJ    = $(C_SRC:.c=.o)

run: $(C_OBJ)
	$(CC) $^ -o run.out
	./run.out

%.o: %.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(C_OBJ) run.out

