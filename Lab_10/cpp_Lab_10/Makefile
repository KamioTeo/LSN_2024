N_PROCESSES = 4

randomdir:= ../RandomNumberGenerator
CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -I$(randomdir)
EXECUTABLES = main.exe

all: 
	make $(EXECUTABLES)

%.exe : %.o random.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)

random.o : $(randomdir)/random.cpp $(randomdir)/random.h
	$(CC) -c $(randomdir)/random.cpp -o random.o $(CFLAGS)

run: $(EXECUTABLES)
	mpiexec -np $(N_PROCESSES) ./main.exe

clean :
	rm *.o main.exe seed.out