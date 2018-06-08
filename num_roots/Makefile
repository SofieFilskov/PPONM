CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.out.txt gsl.out.txt

A.out.txt: A
	./$< > $@
A: A.c gramschmidt.c print_matrix.c newton.c

gsl.out.txt: gsl_first gsl_rosen gsl_himmel
	./$<
	./gsl_rosen
	./gsl_himmel
gsl_first: gsl_first.c
gsl_rosen: gsl_rosen.c
gsl_himmel: gsl_himmel.c


clean:
	rm -f A B *out.txt
