CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

# Specify the target executable and the source files needed to build it
spkmeans: spkmeans.o
	$(CC) -o spkmeans spkmeans.o $(CFLAGS)
# Specify the object files that are generated from the corresponding source files
spkmeans.o: spkmeans.c #check for header file
	$(CC) -c spkmeans.c $(CFLAGS)

clean:
	rm -f spkmeans spkmeans.o