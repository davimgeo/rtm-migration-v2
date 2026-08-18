CC       = gcc
CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma
INCLUDE  = -Iinclude -Isrc
LIBS     = -lm

LIB      = librtm.a
TARGET   = run.out

SRC      = $(shell find src -name "*.c") $(shell find config -name "*.c")
OBJ      = $(SRC:.c=.o)

all: $(TARGET)
	./$(TARGET)

$(LIB): $(OBJ)
	ar rcs $@ $^

$(TARGET): main.o $(LIB)
	$(CC) $(CFLAGS) main.o -L. -lrtm $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm -f $(OBJ) main.o $(LIB) $(TARGET)
