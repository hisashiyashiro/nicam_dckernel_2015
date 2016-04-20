#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N dckernel_horizontal_adv_limiter
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_horizontal_adv_limiter/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_horizontal_adv_limiter.exe .
ln -svf ${HMDIR}/dckernel_horizontal_adv_limiter/data/snapshot.Horizontal_Adv_limiter.001 .

./dckernel_horizontal_adv_limiter.exe
