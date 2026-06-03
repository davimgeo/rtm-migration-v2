CFLAGS = -std=c++11 -Wall

SRC = main.cpp src/geometry.cpp 

run: 
	g++ $(CFLAGS) $(SRC) -o run.out
	./run.out

clean:
	rm -f run.out
