CC = gcc

default: libDispInHalfSpace.a

libexamples.a: DispInHalfSpace.o
	ar rcs $@ $^
    
examples.o: DispInHalfSpace.c
	$(CC) -c $<

clean:
	rm *.o *.a
