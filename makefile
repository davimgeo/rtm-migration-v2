CC          = gcc
AR          = ar

CFLAGS      = -std=gnu11 -O3 -g -march=native -fopenmp -mavx2 -mfma

INCLUDE     = -Iinclude -Isrc

PYTHON_INCLUDE = $(shell python3-config --includes)
NUMPY_INCLUDE = $(shell python3 -c 'import numpy; print(numpy.get_include())')
PYTHON_LIBS = $(shell python3-config --embed --ldflags)
PYTHON_LIBDIR = $(shell python3 -c 'import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')

LIBS        = -lm -fopenmp -lfftw3f

TARGET      = run.out

C_SRC       = $(shell find src -name "*.c") config/config.c main.c
C_OBJ       = $(patsubst %.c,build/%.o,$(C_SRC))

.PHONY: run clean

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(C_OBJ)
	$(CC) $(CFLAGS) \
	$(C_OBJ) \
	$(LIBS) \
	$(PYTHON_LIBS) \
	-Wl,-rpath,$(PYTHON_LIBDIR) \
	-o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(INCLUDE) \
	$(PYTHON_INCLUDE) \
	-I$(NUMPY_INCLUDE) \
	-c $< -o $@

clean:
	rm -rf build $(TARGET)
