CC = g++
CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address
#CFLAGS = -g -std=c17 -ffp-exception-behavior=strict

OBJECTS = main.o flux.o main_loop.o setup.o

.PHONY: clean

tonberry1d: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o tonberry1d -lm

clean:
	-rm *.o
	-rm tonberry1d

main.o: flux.o  main_loop.o setup.o

input.o: ini.o


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

