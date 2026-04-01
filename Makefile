# Makefile for the 'simplify' polygon simplification tool
# Builds a C++17 executable from all sources under src/ with include/ on the path

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude
TARGET   = simplify

# Collect all implementation files from the src directory
SRCS = src/main.cpp \
       src/geometry.cpp \
       src/io_csv.cpp \
       src/simplify.cpp \
       src/validation.cpp

OBJS = $(SRCS:.cpp=.o)

# Default rule: build the simplify executable in this directory
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Pattern rule: compile each .cpp to a .o object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
