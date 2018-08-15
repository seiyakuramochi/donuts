CXX := g++
CXXFLAGS := -Wall -c -std=c++11 -g -ggdb
main: main.o torus.o
	$(CXX) -o main main.o torus.o
main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp
torus.o: torus.cpp
	$(CXX) $(CXXFLAGS) -o torus.o torus.cpp