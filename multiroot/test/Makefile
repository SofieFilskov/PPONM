CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

all: test myecho
test: test.c
myecho: myecho.c

clean:
	$(RM) test myecho
