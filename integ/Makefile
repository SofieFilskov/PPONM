CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: integ hamilton.svg
	./integ

hamilton.svg: hamilton.gpi hamilton.out.txt
	gnuplot $<

hamilton.out.txt: hamilton
	./$< > $@

integ: integ.c

hamilton: hamilton.c

clean:
	rm -f integ hamilton *.txt *.svg
