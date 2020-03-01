all: sph sep

sph:
	gcc *.c -o sph -lm -Wall

sep:
	gcc sep_plot/sep_plot.c -o sep -lm -Wall

clean:
	-rm -f *~ *.o
	-rm -f /sep_plot/*~ *.o

purge: clean
	-rm -f sph
	-rm -f sep