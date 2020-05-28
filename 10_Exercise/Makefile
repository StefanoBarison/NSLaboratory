CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o tsp.o Random.o
	$(CC) main.o tsp.o Random.o -o main.exe
main.o : main.cpp 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
Random.o : Random.cpp Random.h
	$(CC) -c Random.cpp -o Random.o $(CFLAGS)
clean :
	rm *.o main.exe 
