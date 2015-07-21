

#CFLAGS = -I${GMXLDLIB}/../include -g -O3 -fomit-frame-pointer -finline-functions -Wall -msse2 -Wno-unused -fno-exceptions
#CXXFLAGS = -I${GMXLDLIB}/../include -g -O3 -fomit-frame-pointer -finline-functions -Wall -msse2 -Wno-unused -fno-exceptions

CFLAGS = -I${GMXLDLIB}/../include -g -Wall -Wno-unused -fno-exceptions -O3
CXXFLAGS = -I${GMXLDLIB}/../include -g  -Wall -Wno-unused -fno-exceptions -O3

#CFLAGS = -I${GMXLDLIB}/../include -g -Wall -msse2 -funroll-all-loops -Wno-unused -fno-exceptions
#CXXFLAGS = -I${GMXLDLIB}/../include -g -Wall -msse2 -funroll-all-loops -Wno-unused -fno-exceptions

LDFLAGS = -L${GMXLDLIB}  -g
LIBS = -lgmx -lmd -lpthread


timeproc: timeproc.o proccom.o event.o hist.o grid.o
	$(CXX) $^ ${LDFLAGS} ${LIBS} -o $@


clean:
	rm *.o

depend: 
	makedepend -Y. -f Makefile *.cpp

# DO NOT DELETE

event.o: hist.h event.h grid.h proccom.h
grid.o: grid.h
hist.o: hist.h
proccom.o: hist.h event.h grid.h proccom.h
timeproc.o: hist.h event.h grid.h proccom.h
