CC = g++
CFLAGS = -Wall -O3 --std=c++11

1_0_Exercise.exe : 1_0_Exercise.o random.o
	$(CC) random.o 1_0_Exercise.o -o 1_0_Exercise.exe
Buffon.exe : Buffon.o random.o
		$(CC) random.o Buffon.o -o Buffon.exe
1-0_Exercise.o : 1_0_Exercise.cpp
	$(CC) -c 1_0_Exercise.cpp -o 1_0_Exercise.o $(CFLAGS)
Buffon.o : Buffon.cpp
	$(CC) -c Buffon.cpp -o Buffon.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o 1_0_Exercise.exe seed.out
