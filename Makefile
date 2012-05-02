#for quick tests, no optimization...to enable debug messages
#which may be printed by the Qubit classes, etc delete -DNODEBUG
CC			= g++
CFLAGS	= -O2 -g -DNODEBUG	
LNKOPT	= -L. -lOpenQubit
PERCEPS	= templates/perceps
PEROPT	= -h -a -b -e -m -r -t templates/ 

all: qubit

clean:
	rm -f *.o *.a

docs:
	$(PERCEPS) $(PEROPT) -d doc/

libOpenQubit.a: utility.o qstate.o qop.o iomanip.o
	ar rc libOpenQubit.a utility.o qstate.o qop.o iomanip.o
	ranlib libOpenQubit.a

utility.o: utility.cc utility.h
	$(CC) $(CFLAGS) -c utility.cc

qstate.o: qstate.cc qstate.h
	$(CC) $(CFLAGS) -c qstate.cc

iomanip.o: iomanip.cc
	$(CC) $(CFLAGS) -c iomanip.cc

qop.o: utility.o qop.cc qop.h
	$(CC) $(CFLAGS) -c qop.cc

qubit: main.cc libOpenQubit.a
	$(CC) $(CFLAGS) main.cc -o shor $(LNKOPT)
