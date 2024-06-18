randomdir:= ../RandomNumberGenerator

CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I$(randomdir)

main : main.o random.o
	$(CC) random.o main.o -o main

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : $(randomdir)/random.cpp $(randomdir)/random.h
	$(CC) -c $(randomdir)/random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main seed.out