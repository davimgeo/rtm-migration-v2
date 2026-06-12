CC       = gcc

CFLAGS   = -std=gnu99 -g

C_SRC    = src/utils.c main.c src/model/extent.c src/io.c src/model/create.c

INCLUDE =  -I/home/lum/coding/rtm-migration-v2/rtm-migration-c/include

C_OBJ    = $(C_SRC:.c=.o)

run: $(C_OBJ)
	$(CC) $^ -o run.out
	./run.out
	$(MAKE) clean

%.o: %.c
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(C_OBJ) run.out


