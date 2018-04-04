#! /bin/bash

cd ../utils
make
cd ../sticks
make
g++  -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0 -g3 -Wall -fmessage-length=0   -D'_PRODART_VERSION_="0.1"' -I/Users/michaelsadowski/me/prodart2/src  -I/opt/local/include/ -O2 -ffast-math -msse3 -funroll-loops -pipe -o stick teststickmatch.o  ../pose/*.o stickmatch.o  ../utils/line_fit.o  -L/opt/local/lib -lgsl -lgslcblas -lm -lboost_system-mt-d

mv stick ../apps
