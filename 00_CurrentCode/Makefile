# Makefile for building Genome Streamlining Model

CC = g++
CFLAGS = -Wall -O2 -L/usr/local/lib
LDFLAGS = -lgsl -lgslcblas -lm 

COMPILE = $(CC) $(CFLAGS) -c

OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: TheModel

TheModel: $(OBJFILES)
	$(CC) -o Model_with_HGT $(OBJFILES) $(LDFLAGS)

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	rm -f *~ *.o 
