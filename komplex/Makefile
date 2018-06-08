CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: default clean
default: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.o komplex.o

main.o komplex.o: komplex.h

clean:
	rm main main.o komplex.o out.txt
