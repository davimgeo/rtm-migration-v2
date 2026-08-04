CC       = gcc
CFLAGS   = -std=gnu99 -O3 -g -march=native -fopenmp -mavx2 -mfma
INCLUDE  = -Iinclude -Isrc
LIBS     = -lm

COMMON_SRC = \
	main.c \
	config/config.c \
	src/internal.c \
	src/geometry/create.c \
	src/geometry/read.c \
	src/model/create.c \
	src/model/extent.c \
	src/propagation/damp.c \
	src/seismogram/load.c \
	src/seismogram/rm_direct_wave.c \
	src/wavelet/create.c

RTM_SRC = \
	src/rtm/forward.c \
	src/rtm/backward.c \
	src/rtm/rtm.c

ACOUSTIC_SRC = \
	src/propagation/acoustic/forward.c

ELASTIC_SRC = \
	src/propagation/elastic/forward.c

ACOUSTIC = $(COMMON_SRC) $(RTM_SRC) $(ACOUSTIC_SRC)
ELASTIC  = $(COMMON_SRC) $(RTM_SRC) $(ELASTIC_SRC)

.PHONY: acoustic elastic run clean

acoustic:
	$(CC) $(CFLAGS) $(INCLUDE) $(ACOUSTIC) $(LIBS) -o run.out

elastic:
	$(CC) $(CFLAGS) $(INCLUDE) $(ELASTIC) $(LIBS) -o run.out

run: acoustic
	./run.out

clean:
	rm -f run.out
