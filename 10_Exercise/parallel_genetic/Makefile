CFLAGS = -Wall -O3 --std=c++11

main_p.exe : main_p.o tsp.o Random.o
	mpicxx main_p.o tsp.o Random.o -o main_p.exe
main_p.o : main_p.cpp
	mpicxx -c main_p.cpp -o main_p.o $(CFLAGS)
tsp.o: tsp.cpp tsp.h
	mpicxx -c tsp.cpp -o tsp.o $(CFLAGS)
Random.o : Random.cpp Random.h 
	mpicxx -c Random.cpp -o Random.o $(CFLAGS)
clean :
	rm *.o main_p.exe seed.out