FC=gfortran

OPTS=-O3
LIBS=-L/usr/lib64 -llapack -lblas

all: constants.o parameters.o helper.o integrator.o main.o redfield1b

constants.o: constants.f95
	${FC} ${OPTS} -c constants.f95

parameters.o: parameters.f95
	${FC} ${OPTS} -c parameters.f95

helper.o: helper.f95
	${FC} ${OPTS} -c helper.f95

integrator.o: integrator.f95
	${FC} ${OPTS} -c integrator.f95

main.o: main.f95
	${FC} -c main.f95

redfield1b: main.o integrator.o helper.o parameters.o constants.o
	${FC} ${OPTS} ${LIBS} -o redfield1b main.o integrator.o helper.o parameters.o constants.o

clean:
	rm -f *.mod *.o refield1b *.dat
