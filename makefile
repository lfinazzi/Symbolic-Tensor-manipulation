CC = gcc
CFLAGS = -Wall -O3
LIBS = -lm
GLIBS =
GLIBS +=
OBJECTS = main.o Einsum.o
HEADERS =
DEBUG_FLAG = -g

ALL : main.exe
	@echo File has been successfully compiled $(NEWLINE)

main.exe : $(OBJECTS)
	$(CC) $(OBJECTS) $(DEBUG_FLAG) -o main.exe $(LIBS) $(GLIBS) $(CFLAGS)

main.o : main.c $(HEADERS)
	$(CC) -c main.c $(DEBUG_FLAG) -o main.o $(CFLAGS)
	
Einsum.o : main.c $(HEADERS)
	$(CC) -c Einsum.c $(DEBUG_FLAG) -o Einsum.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
