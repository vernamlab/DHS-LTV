CC = g++
LD = g++
AR = ar

CFLAGS = -O3 -std=c++11 -I/usr/local/include -I..
GMP=-lgmp 
LDLIBS = -I/usr/local/include -L/usr/local/lib -lntl $(GMP) -lm

all: main

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

main: main.cpp fft_mult.cpp general.cpp ltv.cpp
	$(CC) $(CFLAGS) -o $@ $< fft_mult.cpp general.cpp ltv.cpp $(LDLIBS)

