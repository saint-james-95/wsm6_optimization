#!/bin/bash
chunk=$1
nz=$2
echo "chunk = $chunk, nz = $nz";
#./makewsm6 arch=intelxeon threading=yes chunk=4 nz=32 fpmp=no >& make.xeon.out
run1="./makewsm6 arch=intelxeon threading=yes chunk=$chunk nz=$nz fpmp=no "
echo "Run: $run1"
$run1 >& make.xeon.out
run2="cd ../run_intelxeon_ompyes_CHUNK"$chunk"_NZ"$nz"_FPMPno"
echo "Run: $run2"
$run2
run3="./runintelxeon threads="$nz
echo "Run: $run3"
$run3 > out
run4="cd "$(head -1 out | cut -d':' -f 2)
echo "Dir: $run4"
$run4
run5="./runscript"
echo "Run: Rrun5"
$run5
run6="cat timing.summary"
echo "Run: $run6"
$run6
run7="cat stdout"
echo "Run: $run7"
$run7
