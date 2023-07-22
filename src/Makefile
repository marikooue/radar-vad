CC = gcc

SRCS = vad2.c marumi.c 

OBJS = ${SRCS:.c=.o}

HEAD = marumi.h

EXE = VAD2

CFLAGS = -O4 -Wall -g

all : ${EXE}

${OBJS}: ${HEAD}

${EXE} : ${OBJS} 
	${CC} ${OBJS} -L/usr/local/lib -lm  -lnetcdf -o $@

${OBJS} : ${SRCS}

clean : 
	\rm -f *.c~ *.h~ *.o Makefile~ core

