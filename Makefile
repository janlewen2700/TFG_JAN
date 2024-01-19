#suffixes used

.SUFFIXES: .o .h .c 

#compilation options

COMPOPS = -g -Wall

#linking options

LINKOPS = -lgsl -lgslcblas -lfftw3 -lm

#list of object files

objects = main.o modelB.o

#list of header files

headers = stdio.h stdlib.h math.h complex.h gsl_rng.h gsl_math.h fftw3.h headers.h

#list of source codes

sources = main.c modelB.c

#directory paths for the source and header files

vpath %.c ./source_code/
vpath %.h ./header_file/
vpath %.h /opt/homebrew/Cellar/gcc/13.2.0/include/c++/13/tr1/ 							#gcc compiler directory. To be changed accordingly based on the machine.
vpath %.h /opt/homebrew/Cellar/fftw/3.3.10_1/include/                       #directory path for fftw3 header function. To be changed accordingly based on the directory path on the machine.
vpath %.h /opt/homebrew/Cellar/gsl/2.7.1/lib/               #lgsl and lgslcblas file directory. To be changed accordingly based on the machine.
vpath %.h /opt/homebrew/Cellar/gsl/2.7.1/include/gsl/                     #gsl header functions directory. To be changed accordingly based on the machine.

#actions on the source files

simrun: $(objects) $(headers)
	gcc -o simrun.out -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -L /opt/homebrew/Cellar/fftw/3.3.10_1/lib/ $(objects) $(LINKOPS)
main.o: $(sources) $(headers)
	gcc -o $@ -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -I /opt/homebrew/Cellar/fftw/3.3.10_1/include/  -c ./source_code/main.c $(COMPOPS)
modelB.o: modelB.c stdio.h stdlib.h math.h gsl_math.h complex.h fftw3.h headers.h
	gcc -I /opt/homebrew/Cellar/fftw/3.3.10_1/include/ -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -c $(COMPOPS) $<

.PHONY : clean CLEAN
clean:
	rm -rf *.o
CLEAN:
	rm -rf CH_variable_mobility.out


#end of the Makefile
