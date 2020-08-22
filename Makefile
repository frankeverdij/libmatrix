# libmatrix Makefile

# gcc/g++ compiler
CC = gcc
CXX = g++

# gcc compiler flags
CFLAGS = -O2 -g -fPIC -DNDEBUG -Wall -Wno-unused-but-set-variable \
         -I/usr/include/suitesparse -I../SuiteSparse/include -I.

# Archiver
AR = ar

# Archiver flags
ARFLAGS = cr

# linker
LD = ld

# linker flags
LDFLAGS = -L../SuiteSparse/lib -lumfpack -lamd -lsuitesparseconfig -llapack -lblas

# name of archive
LIB = libmatrix.a

# list of all objects
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)

# Makefile
all: $(OBJECTS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJECTS)
shared: $(OBJECTS)
	$(LD) -shared -soname,libmatrix.so -o libmatrix.so.1.2 \
		$(OBJECTS) $(LDFLAGS)
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Cleaning everything
clean:
	rm -f $(LIB) $(OBJECTS)
# End
