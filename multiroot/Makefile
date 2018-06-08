CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: rosen.out.txt hydro.svg

rosen.out.txt: rosen
	./$< > $@
rosen: rosen.c

hydro.svg: hydro.gpi hydro.out.txt
	gnuplot $<
hydro.out.txt: hydro
	./$< > $@
hydro: hydro.c

clean:
	$(RM) rosen *.txt hydro *.svg 
