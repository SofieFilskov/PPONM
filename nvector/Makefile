CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

.PHONEY: all clean

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o nvector.o

main.o nvector.o: nvector.h

clean:
	rm out.txt main main.o nvector.o test
