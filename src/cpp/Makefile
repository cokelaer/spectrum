#CFL=-g
#F77=f77 -v
#.f.o: ; $(F77) $(CFL) -c $*.f


.c.o: ; ${CC} ${CFLAGS}  -c $*.c
CFLAGS          =  ${INCLUDE}

CC = gcc -g -Wall -O2

#


all: 
	$(CC) -c test.c -o test.o 
	$(CC) -c -shared *.o -o libtest.so 

clean:    
	rm -f *.o *.so
