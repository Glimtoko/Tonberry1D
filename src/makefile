CC = mpic++
#CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address
CFLAGS = -O3 -std=c++17

OBJECTS = main.o flux.o main_loop.o setup.o comms.o

BOOST_DIR=/prod/boost
BOOST_INC=$(BOOST_DIR)/include
BOOST_LIB=$(BOOST_DIR)/lib

INC = -I$(BOOST_INC)
LIB = -L$(BOOST_LIB) -lboost_program_options -lboost_mpi

.PHONY: clean

tonberry1d: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o tonberry1d -lm $(LIB)

clean:
	-rm *.o
	-rm tonberry1d

setup.o: comms.o

main_loop.o: comms.o

main.o: flux.o  main_loop.o setup.o

input.o: ini.o


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

