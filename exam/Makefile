CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.out.txt
A.out.txt: main
	./$< > $@
main: main.c print_matrix.c jacobi.c 



clean:
	rm -f A B *.svg *.txt
