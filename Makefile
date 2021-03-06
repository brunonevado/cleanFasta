
PROGRAM=cleanFasta


CC=g++
LIBS=-lbpp-core -lbpp-seq -lbpp-phyl


CFLAGS=-Wall -O3 -std=c++11

# use below if bio++ libs/includes are not in default location 
# IDIR =/home/Bruno/local/include/bpp
# LDIR =/home/Bruno/local/lib/bpp/
# CFLAGS=-I$(IDIR) -L$(LDIR) -Wall -O3 -std=c++11


$PROGRAM: source/main.cpp source/args.h source/args.cpp 
	$(CC) $(CFLAGS) -o $(PROGRAM) source/main.cpp source/args.h source/args.cpp  $(LIBS)
	

