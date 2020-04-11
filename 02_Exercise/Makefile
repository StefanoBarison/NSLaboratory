CC = g++
CFLAGS = -Wall -O3 --std=c++11

Integrals.exe : Integrals.o random.o
	$(CC) random.o Integrals.o -o Integrals.exe
Random_walks.exe : Random_walks.o random.o
		$(CC) random.o Random_walks.o -o Random_walks.exe
Integrals.o : Integrals.cpp
	$(CC) -c Integrals.cpp -o Integrals.o $(CFLAGS)
Random_walks.o : Random_walks.cpp
	$(CC) -c Random_walks.cpp -o Random_walks.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Integrals.exe seed.out
