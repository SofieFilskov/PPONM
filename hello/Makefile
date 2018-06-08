CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean
all: hello
	./$<

hello: hello.c

clean:
	rm -f hello
