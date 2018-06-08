CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean
all: report.pdf

report.pdf: report.tex plot.tex
	latex report
	latex report
	dvipdf report
plot.tex: plot.gpi out.txt
	gnuplot $<

out.txt: main
	./$< 100 > $@
main: main.c

clean:
	rm -f main plot.tex *.pdf *out.txt *.aux *.dvi *.log
