CFLAGS = -std=c++11

SRC = main.cpp src/geometry.cpp src/model.cpp

TESTS = tests/parallel.cpp

run: 
	g++ $(CFLAGS) $(SRC) $(TESTS) -o run.out
	./run.out

clean:
	rm -f run.out
