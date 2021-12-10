CFLAGS = -std=c++17 -O3 -Wall 
OBJS = main.o analysis.o

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o onco_model
	
analysis.o: ../Simulation/analysis.cpp
	g++ ${CFLAGS} -c ../Simulation/analysis.cpp	

main.o: main.cpp
	g++ ${CFLAGS} -c main.cpp

clean:
	rm -f onco_model ${OBJS}
