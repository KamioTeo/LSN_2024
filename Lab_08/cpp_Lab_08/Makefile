randomdir:= ../RandomNumberGenerator

CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I$(randomdir)

EXECUTABLES = delta_finder.exe SA.exe

all: 
		make $(EXECUTABLES)

%.exe : %.o random.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

random.o : $(randomdir)/random.cpp $(randomdir)/random.h
	$(CC) -c $(randomdir)/random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o delta_finder.exe SA.exe seed.out