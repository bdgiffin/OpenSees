# Simple Makefile used to compile a collection of C++ API functions
# contained in the specified $(LIBNAME).cpp as a "shared object" library.

LIBNAME=ParticleDynamicsAPI

CC=g++
FLAGS=-O2 -Wall -pedantic -shared -fPIC -lstdc++ -std=c++11

.PHONY: clean $(LIBNAME).so debug

all: $(LIBNAME).so

debug: $(LIBNAME).cpp
	$(CC) $(FLAGS) -DDEBUG_ENABLE -o $(LIBNAME).so $(LIBNAME).cpp

$(LIBNAME).so: $(LIBNAME).cpp
	$(CC) $(FLAGS) -o $(LIBNAME).so $(LIBNAME).cpp

clean:
	rm -f $(LIBNAME).so
