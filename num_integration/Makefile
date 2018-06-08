CFLAGS = -Wall -std=c99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: A.out.txt B.out.txt C.out.txt

A.out.txt B.out.txt C.out.txt: A B C
	./$< > A.out.txt
	./B > B.out.txt
	./C > C.out.txt

A: A.c integration.c

B: B.c integration.c

C: C.c integration.c

clean:
	rm -f A B C *out.txt
