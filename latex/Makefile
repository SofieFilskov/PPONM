CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
ABDX = -3 3 0.2

.PHONEY: all clean

all: errfunc.pdf

errfunc.pdf: errfunc.tex plot.tex
	latex errfunc
	latex errfunc
	dvipdf errfunc

plot.tex: errfunc.gpi data.out.txt
	gnuplot $<

data.out.txt: errfunc
	./$< $(ABDX) > $@

errfunc: errfunc.c



clean:
	rm -f errfunc *.svg *out.txt *.aux *.dvi *.log plot.tex
