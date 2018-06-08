CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: rosen.svg radio.svg
rosen.svg: rosen.gpi rosen.out.txt
	gnuplot $<
rosen.out.txt: rosen
	./$< > $@
rosen: rosen.c

radio.svg: radio.gpi radio.out.txt
	gnuplot $<
radio.out.txt: radio
	./$< > $@
radio:radio.c

clean:
	$(RM) rosen radio *.txt *.svg
