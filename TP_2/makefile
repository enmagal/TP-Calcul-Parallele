CC = gcc
Exec = prog
Src = $(wildcard *.c)
Obj = $(Src:.c=.o)
CFLAGS = -Wall -Wextra
LDFLAGS = -lm

all : $(Exec)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< -g -fopenmp

$(Exec) : $(Obj)
	$(CC) -o $@ $^ $(LDFLAGS) -fopenmp

clean :
	rm -rf *.o
	
mrproper : clean
	rm -rf prog