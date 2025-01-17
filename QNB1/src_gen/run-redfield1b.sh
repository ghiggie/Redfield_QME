#!/bin/bash

# Name of the executable
EXE="Redfield1B"
# Size of the Hilbert space
N="4"
# Name of the parameter file
PAR="params${N}.template"
# Execution options
OPT=""
# Execution time
T=2000
# Temperate and gamma
temp=1.0
gamma=`echo "0.1*$temp" | bc`
# Number of Matsubara terms
mat=5
# Coupling strength between system and environment
LI=0.01
LF=0.01
DL=0.01
# Coupling strength between atoms
KI=0.00
KF=0.00
DK=0.01

mkdir Data

for l in `seq $LI $DL $LF`;
do
   for k in `seq $KI $DK $KF`;
   do
      ktemp=`echo "4*$k" | bc`
      K=`printf "%3.2f" $ktemp`
      L=`printf "%3.2f" $l`
      sed -e "s/KAPPA/$K/g" -e "s/LAMBDA/$L/g" -e "s/TIME/$T/g" $PAR > params.tmp;
      sed -e "s/TEMP/$temp/g" -e "s/GAMMA/$gamma/g" -e "s/MATSU/$mat/g" params.tmp > params.cfg
      rm -f params.tmp
      ./${EXE} ${OPT} ${N} params.cfg

      IDX="-L$L-K$K"

      mv ${EXE}.out Data/${EXE}${IDX}.out
      mv params.cfg Data/params${N}${IDX}.cfg
      mv entropy.dat Data/entropy${IDX}.dat
      mv eta.dat Data/eta${IDX}.dat
      mv fidelity.dat Data/fidelity${IDX}.dat
      mv heat.dat Data/heat${IDX}.dat
      mv hermitian.dat Data/hermitian${IDX}.dat
      mv lambda.dat Data/lambda${IDX}.dat
      mv positivity.dat Data/positivity${IDX}.dat
      mv rho.dat Data/rho${IDX}.dat
      mv rho_E.dat Data/rho_E${IDX}.dat
      mv sprod.dat Data/sprod${IDX}.dat
      mv stats_E.dat Data/stats_E${IDX}.dat
      mv stats_XS.dat Data/stats_XS${IDX}.dat
      mv trace.dat Data/trace${IDX}.dat
      mv trace_distance.dat Data/trace_distance${IDX}.dat

      if [ $(echo "${N}==4" | bc) -eq 1 ]; then
         mv rhoA.dat Data/rhoA${IDX}.dat
         mv rhoB.dat Data/rhoB${IDX}.dat
         mv rhoA_E.dat Data/rhoA_E${IDX}.dat
         mv rhoB_E.dat Data/rhoB_E${IDX}.dat
      fi
   done
done
