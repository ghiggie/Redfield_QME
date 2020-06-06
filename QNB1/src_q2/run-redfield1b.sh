#!/bin/bash

# Name of the executable
EXE="Redfield1B"
# Name of the parameter file
PAR="params.template"
# Execution options
OPT="-G"
# Execution time
T=2000
# Temperate and gamma
temp=2.0
gamma=`echo "0.1*$temp" | bc -l`
# Number of Matsubara terms
mat=1
# Coupling strength between system and environment
LI=0.05
LF=0.05
DL=0.002
# Coupling strength between atoms
KI=5.00
KF=5.00
DK=0.02

mkdir Data

for l in `seq $LI $DL $LF`;
do
   for k in `seq $KI $DK $KF`;
   do
      ktemp=`echo "4*$k" | bc -l`
      K=`printf "%4.3f" $ktemp`
      L=`printf "%4.3f" $l`
      sed -e "s/KAPPA/$K/g" -e "s/LAMBDA/$L/g" -e "s/TIME/$T/g" $PAR > params.tmp;
      sed -e "s/TEMP/$temp/g" -e "s/GAMMA/$gamma/g" -e "s/MATSU/$mat/g" params.tmp > params.cfg
      rm -f params.tmp
      ./${EXE} ${OPT} params.cfg

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
      mv rhoA.dat Data/rhoA${IDX}.dat
      mv rhoB.dat Data/rhoB${IDX}.dat
      mv rhoA_E.dat Data/rhoA_E${IDX}.dat
      mv rhoB_E.dat Data/rhoB_E${IDX}.dat
      mv rhoA_fid.dat Data/rhoA_fid${IDX}.dat
      mv rhoB_fid.dat Data/rhoB_fid${IDX}.dat
      mv rhoA_traced.dat Data/rhoA_traced${IDX}.dat
      mv rhoB_traced.dat Data/rhoB_traced${IDX}.dat
   done
done
