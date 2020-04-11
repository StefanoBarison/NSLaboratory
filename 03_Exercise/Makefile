CC = g++
CFLAGS = -Wall -O3 --std=c++11

Vanilla_prices.exe : Vanilla_prices.o random.o
	$(CC) random.o Vanilla_prices.o -o Vanilla_prices.exe
Vanilla_prices.o : Vanilla_prices.cpp
	$(CC) -c Vanilla_prices.cpp -o Vanilla_prices.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Vanilla_prices.exe seed.out
