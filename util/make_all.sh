#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/src
      make || exit 1
   cd -
done
