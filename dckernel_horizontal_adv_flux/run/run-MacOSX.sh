#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_horizontal_adv_flux/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_horizontal_adv_flux.exe .

./dckernel_horizontal_adv_flux.exe
