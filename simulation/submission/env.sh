module add languages/intel-compiler-16
export CFLAGS="-O3 -qopenmp -no-prec-div -ipo -xHOST -unroll-aggressive -m64 -auto-p32 -no-prec-sqrt -fp-model fast=2"
export CC="icc"
