CC = mpicc
Exec = prog
Src = $(wildcard *.c)
Obj = $(Src:.c=.o)
CFLAGS = -Wall -Wextra
LDFLAGS = -lm

all : $(Exec)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< -g

$(Exec) : $(Obj)
	$(CC) -o $@ $^ $(LDFLAGS)

clean :
	rm -rf *.o
	
mrproper : clean
	rm -rf prog