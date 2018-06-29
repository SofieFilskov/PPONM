CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.svg

A.svg: A.gpi A.out.txt
	gnuplot $<

A.out.txt: A
	./$< > $@
A: A.c ann.c



clean:
	rm -f A *.svg *out.txt
