CFLAGS = -std=c++11 -Wall

SRC = main.cpp src/geometry.cpp

all:
	g++ $(CFLAGS) $(SRC) -o run.out

run: all
	./run.out

clean:
	rm -f run.out
