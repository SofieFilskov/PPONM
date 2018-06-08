CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.out.txt B.out.txt

A.out.txt: A
	./$< > $@
A: A.c gramschmidt.c print_matrix.c
B.out.txt: B
	./$< > $@
B: B.c gramschmidt.c print_matrix.c

clean:
	rm -f A B *out.txt
