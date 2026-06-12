CC      = gcc
CFLAGS  = -std=gnu99 -g
INCLUDE = -I/home/lum/rtm-migration-v2/include

# Sources
C_SRC := $(shell find . -name "*.c")
C_OBJ := $(patsubst ./%.c,obj/%.o,$(C_SRC))

run: $(C_OBJ)
	$(CC) $^ -lm -o run.out
	./run.out
	$(MAKE) clean

obj/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean:
	rm -rf obj run.out

