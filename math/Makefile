CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: default clean
default: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math: math.c

clean:
	rm -f math out.txt
