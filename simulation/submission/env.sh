module add languages/intel-compiler-16
module add openmpi/intel/64/1.6.5
module add tools/intel_builds/tau-2.23.1-openmpi
export CFLAGS="-O3 -qopenmp -no-prec-div -ipo -xHOST -unroll-aggressive -m64 -auto-p32 -no-prec-sqrt -fp-model fast=2"
export CC="mpicc"
