CC = g++
CFLAGS = -Wall -O3 --std=c++11

mod_NVE.exe : mod_NVE.o 
	$(CC)  mod_NVE.o -o mod_NVE.exe
mod_NVE.o : mod_NVE.cpp
	$(CC) -c mod_NVE.cpp -o mod_NVE.o $(CFLAGS)
clean :
	rm *.o mod_NVE.exe
