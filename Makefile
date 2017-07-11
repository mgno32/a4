TARGET := test

CXX := g++
CFLAGS := -g -std=c++11
LDFLAGS :=  -pthread -lX11

SRCS := $(wildcard ./tutorial.cpp)
OBJS := $(patsubst %cpp,%o,$(SRCS))

all: $(OBJS) 
	$(CXX) $^ $(CFLAGS) $(LDFLAGS) -o $(TARGET)
clean:
	rm $(TARGET) $(OBJS)
%.o:%.cpp  
	$(CXX) $< $(CFLAGS) $(LDFLAGS) -c -o $@ 
