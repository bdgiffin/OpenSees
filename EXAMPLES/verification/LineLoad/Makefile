# Simple Makefile used to compile a collection of C++ API functions
# contained in the specified $(LIBNAME).cpp as a "shared object" library.

LIBNAME=line_load_example

CC=g++
FLAGS=-O2 -Wall -pedantic -shared -fPIC -lstdc++ -std=c++11

.PHONY: clean $(LIBNAME).so

all: $(LIBNAME).so

$(LIBNAME).so: $(LIBNAME).cpp
	$(CC) $(FLAGS) -o $(LIBNAME).so $(LIBNAME).cpp

clean:
	rm -f $(LIBNAME).so
