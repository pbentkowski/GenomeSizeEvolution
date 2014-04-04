#!/bin/bash

g++ -O2 -Wall -o TheModel -L/usr/local/lib TheModel.cpp \
genotype.cpp gene.cpp cell.cpp ecosystem.cpp rngEngine.cpp \
./tagging_system.cpp -lgsl -lgslcblas -lm
