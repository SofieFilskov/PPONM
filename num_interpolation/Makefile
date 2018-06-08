CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: linear.svg linearinteg.svg qua.svg quainteg.svg quaderiv.svg cubic.svg

linear.svg linearinteg.svg: linear.gpi linear.out.txt
	gnuplot $<
linear.out.txt: linear
	./$< > $@
linear: linear.c

qua.svg quainteg.svg quaderiv.svg: qua.gpi qua.out.txt
	gnuplot $<
qua.out.txt: qua
	./$< > $@
qua: qua.c quaspline.c

cubic.svg: cubic.gpi cubic.out.txt
	gnuplot $<
cubic.out.txt: cubic
	./$< > $@
cubic: cubic.c cspline.c

clean:
	rm -f linear cubic qua quaB2 *.svg *out.txt
