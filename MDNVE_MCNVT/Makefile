CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o NVE_NVT.exe

clean_dat:
	rm *.dat

mv_solid_es4:
	mv *.dat ../Lab_04/Solido
mv_liquid_es4:
	mv *.dat ../Lab_04/Liquido
mv_gas_es4:
	mv *.dat ../Lab_04/Gas

mv_solid_es7:
	mv *.dat ../Lab_07/Solido
mv_liquid_es7:
	mv *.dat ../Lab_07/Liquido
mv_gas_es7:
	mv *.dat ../Lab_07/Gas