#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/run
      sh run-MacOSX.sh || exit 1
   cd -
done
