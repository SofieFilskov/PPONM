CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.out.txt B.svg

A.out.txt: A
	./$< > $@
A: A.c gramschmidt.c newton.c

B.svg: B.gpi B.out.txt
	gnuplot $<
B.out.txt: B
	./$< > $@
B: B.c gramschmidt.c qnewton.c print_matrix.c


clean:
	rm -f A B *out.txt
