CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

.PHONEY: all clean 
	
all: test.io.out.txt test.cmd.out.txt io.out.txt cmd.out.txt plot.svg
	cat io.out.txt cmd.out.txt	

test.io.out.txt: main-stdio
	echo 1 2 3 4 5 | ./main-stdio > test.io.out.txt

test.cmd.out.txt: main-cmdline
	./main-cmdline 1 2 3 4 5 6 > test.cmd.out.txt

io.out.txt: main-stdio input.data.txt
	cat input.data.txt | ./main-stdio > io.out.txt

cmd.out.txt: main-cmdline input.data.txt
	./main-cmdline `cat input.data.txt` > cmd.out.txt

input.data.txt: ;
	seq 0 0.2 5 > input.data.txt
	
main-stdio: main-stdio.c

main-cmdline: main-cmdline.c

plot.svg: plot.gpi io.out.txt cmd.out.txt
	gnuplot plot.gpi

clean:
	rm -f main-cmdline main-stdio *.svg *.txt
