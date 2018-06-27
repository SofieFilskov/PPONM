CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: main times.svg

times.svg: times.gpi times.txt
	gnuplot $<
times.txt: Makefile A
	rm -f times.txt
	for n in `seq 50 2 100`; do \time -f "$$n %U" -ao times.txt ./A $$n 1>/dev/null; done

main: A.out.txt B.out.txt
A.out.txt: A
	./$< > $@
A: A.c print_matrix.c jacobi.c
B.out.txt: B
	./$< > $@
B: B.c print_matrix.c jacobi_B.c jacobi.c




clean:
	rm -f A B *.svg *.txt
