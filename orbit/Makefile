CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: odiff.svg orbit.svg

odiff.svg: odiff.gpi odiff.out.txt
	gnuplot $<

odiff.out.txt: odiff
	./$< > $@

odiff: odiff.c

orbit.svg: orbit.gpi orbit.out.txt
	gnuplot $<

orbit.out.txt: orbit
	./$< > $@

orbit: orbit.c

clean:
	rm -f odiff orbit *.svg *out.txt
