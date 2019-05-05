CC = gcc

default: libwrapper.a

libwrapper.a: wrapper.o
	ar rcs $@ $^
    
wrapper.o: wrapper.c wrapper.h
	$(CC) -c $<

clean:
	rm *.o *.a
