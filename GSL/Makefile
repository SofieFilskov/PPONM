CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean
all: airy.svg linear
	./linear
airy.svg: airy.gpi airy_data.txt
	gnuplot $<
airy_data.txt: airy
	./$< > $@
airy: airy.c
linear: linear.c

clean:
	$(RM) airy linear *.txt *.svg