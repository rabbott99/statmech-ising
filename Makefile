NAME?=ising
CXX=g++
LDFLAGS=
CPPFLAGS= -g -Wall -std=c++11


CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)

all: $(NAME)

$(NAME): $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp *.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

.PHONY: clean

clean:
	rm -f *.o $(NAME)
