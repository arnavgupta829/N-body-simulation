CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -I include -fopenmp
LDFLAGS = -lm -fopenmp
 
SRCS = src/main.cpp src/QuadTree.cpp src/Simulation.cpp src/SequentialSimulation.cpp src/ParallelSimulation.cpp src/Morton.cpp
TARGET = nbody
 
.PHONY: all clean profile fast
 
all: $(TARGET)
 
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)
 
clean:
	rm -f $(TARGET)
 
profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: clean all
 
fast: CXXFLAGS = -std=c++17 -O3 -march=native -Wall -I include -fopenmp
fast: clean all