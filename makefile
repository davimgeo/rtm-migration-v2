CC       = gcc
AR       = ar

CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma

INCLUDE  = -Iinclude -Isrc

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE  = $(shell python3 -c 'import numpy; print(numpy.get_include())')
PYTHON_LIBS    = $(shell python3-config --embed --ldflags)

LIBS     = -lm -fopenmp -lfftw3f

TARGET   = run.out
STATIC_LIB = libprop.a

C_SRC    = $(shell find src -name "*.c") config/config.c
C_OBJ    = $(patsubst %.c,build/%.o,$(C_SRC))

MAIN_OBJ = build/main.o


all: $(STATIC_LIB) $(TARGET)
	@./$(TARGET)


# Static library
$(STATIC_LIB): $(C_OBJ)
	@$(AR) rcs $@ $^


# Executable linked against libfwi.a
$(TARGET): $(STATIC_LIB) $(MAIN_OBJ)
	@$(CC) $(MAIN_OBJ) \
		-L. -lprop \
		$(LIBS) \
		$(PYTHON_LIBS) \
		-o $@


# Compile source files
build/%.o: %.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(INCLUDE) \
		$(PYTHON_INCLUDE) \
		-I$(NUMPY_INCLUDE) \
		-c $< -o $@


clean:
	rm -rf build $(TARGET) $(STATIC_LIB)
