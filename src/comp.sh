#!/bin/bash
./makewsm6 clean
OUT="make.xeon.out"
ORIGINDIR=$(pwd)
if [ "$#" -eq 0 ]; then
	echo "Using default arguments: ./makewsm6 arch=intelxeon threading=yes chunk=4 nz=32 fpmp=no >& $OUT"
	./makewsm6 arch=intelxeon threading=yes chunk=4 nz=32 fpmp=no >& $OUT
else
	./makewsm6 $1 $2 $3 $4 $5 $6 $7 >& $OUT
fi

if [[ $(grep -i error $OUT) ]]; then
	echo "Error in compiling. See $OUT for details."
	exit -1
fi

SRCDIR=$(grep "SRCDIR = " $OUT | sed 's|.*= \(.*\)|\1|')
RUNDIR1=$(grep -A1 "Run scripts are in:" $OUT | tail -1)
RUNSCRPT=$(grep -A1 "Run script is" $OUT | tail -1 | sed 's/^\s*//')
cd $RUNDIR1
./$RUNSCRPT > $RUNSCRPT.out
RUNDIR2=$(grep "Made directory" $RUNSCRPT.out | sed 's/Made\sdirectory:\s//')
cd $RUNDIR2
./runscript
#cat stdout
echo "Running a diff on output with refstdout"
diff stdout $ORIGINDIR/refstdout
exit 0
