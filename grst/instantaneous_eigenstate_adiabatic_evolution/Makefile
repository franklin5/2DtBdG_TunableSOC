CC=icpc
eigenbin = /home/ld7/bin	
CFLAGS=-c -Wall -I${eigenbin}
LFLAGS=-limf -lm
all: tBdG

tBdG: main.o grst.o lgwt.o gapnumberintegral.o tBdG.o
	$(CC) main.o grst.o lgwt.o gapnumberintegral.o tBdG.o -o tBdG $(LFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

gapnumberintegral.o: gapnumberintegral.cpp
	$(CC) $(CFLAGS) gapnumberintegral.cpp

lgwt.o: lgwt.cpp
	$(CC) $(CFLAGS) lgwt.cpp

grst.o: grst.cpp
	$(CC) $(CFLAGS) grst.cpp

tBdG.o: tBdG.cpp
	$(CC) $(CFLAGS) tBdG.cpp

touch: 
	touch *.cpp *.h
clean:
	rm *.o chern *~ *#