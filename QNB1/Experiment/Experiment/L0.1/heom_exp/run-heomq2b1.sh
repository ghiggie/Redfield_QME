#!/bin/bash

### PGI and other modules on cheaha
#module add PGI/18.5-GCC-6.4.0-2.28
#module add GCC/6.4.0-2.28
#module add lapack/gcc/64/3.8.0
#module add openblas/dynamic/0.2.20
#module add cuda80/toolkit/8.0.61

### PGI modules for desktop
#module load /opt/pgi/modules/2019



### setting for ACC
### specify ACC mode:  "core" or "gpu"
#ACC="cpu"
ACC="gpu"

### for multicores, specify the number of cores if necessary
#export ACC_NUM_CORES=24

### Specify the location of the executable.
### Example: On puma  it should be "~/bin/puma"
BIN="~/bin/puma"

### Name of the executable
EXE="heomq2b1"

### Tamplate of the params
PAR="params-${EXE}.template"

# OPT
#      (NULL) = use the initial rho defined in params.cfg
#      -G     = use Gibbs state as initial rho
#      -N     = normalized inituial rho
OPT=""

T1=2000
T2=500

KI=0.2
KF=5.0
DK=0.2

l=0.00

for k in 0.01 0.05 0.1 `seq $KI $DK $KF`;
do
   K=`printf "%3.2f" $k`
   L=`printf "%3.2f" $l`
   if [ $(echo "if (${k} < 0.10 || ${l} < 0.10) 1 else 0" | bc) -eq 1 ]; then
      sed -e "s/KAPPA/$K/g" -e "s/LS/$L/g" -e "s/TIME/$T1/g" $PAR > params.cfg;
      ${EXE}-${ACC} ${OPT} params.cfg
   else
      sed -e "s/KAPPA/$K/g" -e "s/LS/$L/g" -e "s/TIME/$T2/g" $PAR > params.cfg;
      ${EXE}-${ACC} ${OPT} params.cfg
   fi

   IDX="-K$K-L$L"
   mv ${EXE}.out ${EXE}${IDX}.out
   mv rho.dat rho${IDX}.dat
   mv eta.dat eta${IDX}.dat
   mv heat.dat heat${IDX}.dat
   mv thermo.dat thermo${IDX}.dat
   mv sprod.dat sprod${IDX}.dat
   mv params.cfg params${IDX}.cfg

done

