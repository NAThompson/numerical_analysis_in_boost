all:
	g++ --std=c++14 -O3 -I../boost/ -o examples.x main.cpp -lfftw3

clean:
	rm -f *.x
