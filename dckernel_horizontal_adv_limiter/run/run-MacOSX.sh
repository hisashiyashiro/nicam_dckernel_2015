#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_horizontal_adv_limiter/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_horizontal_adv_limiter.exe .
ln -svf ${HMDIR}/dckernel_horizontal_adv_limiter/data/snapshot.Horizontal_Adv_limiter.001 .

./dckernel_horizontal_adv_limiter.exe
