init:
	. /opt/intel/oneapi/setvars.sh

compile:
	icc -lpapi -O3 -fp-model precise -o w_strassen_morton w_strassen_morton.c