FC=gfortran
#FC=pgfortran

OPTS=-O3
#OPTS=-O3 -Minfo

LIBS=-L/usr/lib64 -llapack -lblas

all: constants.o variables.o helper.o bath.o main.o Redfield1B

constants.o: constants.f95
	${FC} ${OPTS} -c constants.f95

variables.o: variables.f95
	${FC} ${OPTS} -c variables.f95

helper.o: helper.f95
	${FC} ${OPTS} -c helper.f95

bath.o: bath.f95
	${FC} ${OPTS} -c bath.f95

main.o: main.f95
	${FC} ${OPTS} -c main.f95

Redfield1B: main.o constants.o variables.o helper.o bath.o
	${FC} ${OPTS} -o Redfield1B main.o constants.o variables.o helper.o bath.o ${LIBS}

clean:
	rm -rf *.o *.mod Redfield1B *.dat *.out
