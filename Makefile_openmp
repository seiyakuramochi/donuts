CXX := g++
CXXFLAGS := -Wall -c -std=c++11 -O
main: main.o torus.o
	$(CXX) -o main main.o torus.o -fopenmp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -fopenmp
torus.o: torus.cpp
	$(CXX) $(CXXFLAGS) -o torus.o torus.cpp
