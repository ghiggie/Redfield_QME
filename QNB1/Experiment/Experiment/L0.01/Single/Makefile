FC=gfortran

OPTS=-O3
LIBS=-L/usr/lib64 -llapack -lblas

all: constants.o parameters.o helper.o bath.o main.o BornMarkov1B

constants.o: constants.f95
	${FC} ${OPTS} -c constants.f95

parameters.o: parameters.f95
	${FC} ${OPTS} -c parameters.f95

helper.o: helper.f95
	${FC} ${OPTS} -c helper.f95

bath.o: bath.f95
	${FC} ${OPTS} -c bath.f95

main.o: main.f95
	${FC} ${OPTS} -c main.f95

BornMarkov1B: main.o constants.o parameters.o helper.o bath.o
	${FC} ${OPTS} ${LIBS} -o BornMarkov1B main.o constants.o parameters.o helper.o bath.o

clean:
	rm -rf *.o *.mod BornMarkov1B
