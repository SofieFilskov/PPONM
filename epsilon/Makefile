CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: default clean
default: out.txt
	cat out.txt

out.txt: epsilon
	./epsilon > out.txt

epsilon: epsilon.c

clean:
	rm -f epsilon *out.txt
