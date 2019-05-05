CC = gcc

default: libexamples.a

libexamples.a: examples.o
	ar rcs $@ $^
    
examples.o: examples.c examples.h
	$(CC) -fPIC -c $<

clean:
	rm *.o *.a
