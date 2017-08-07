CC = gcc
FILE = kchi.c kchi.h rand.c
CFLAGS = -Wall -Wextra -g -lm #-pg
all: Makefile main rand
main: kchi.h func.h algo.h kchi.o func.o algo.o
	$(CC) kchi.o func.o algo.o -o kchi $(CFLAGS)	
rand: rand.c func.o
	$(CC) rand.c func.o -o rand $(CFLAGS)
clean:
	rm -f kchi func algo rand *.o *~

