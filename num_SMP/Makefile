CFLAGS = -Wall -std=c99 $$(gsl-config --cflags) -fopenmp
LDLIBS = $$(gsl-config --libs) -lgomp

.PHONEY: all clean

all: A.out.txt

A.out.txt: A
	./$< > $@

A: A.c integration.c


clean:
	rm -f A *out.txt
