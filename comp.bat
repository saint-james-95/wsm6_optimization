chunk=$1
nz=$2
echo "chunk = $chunk, nz = $nz";
#./makewsm6 arch=intelxeon threading=yes chunk=4 nz=32 fpmp=no >& make.xeon.out
run1="./src/makewsm6 arch=intelxeon threading=yes chunk=$chunk nz=$nz fpmp=no >& make.xeon.out"
echo "Run: $run1"
$run1
