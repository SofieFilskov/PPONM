CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: sinus.svg function.svg

sinus.svg function.svg: A.gpi A.out.txt
	gnuplot $<

A.out.txt: A
	./$< > $@
A: A.c ode.c


clean:
	rm -f A *out.txt *.svg
