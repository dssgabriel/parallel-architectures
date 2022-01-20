## ------------------------------------------------------------------------- ##
## This Makefile builds various versions of a 3D Nbody simulation.           ##
##                                                                           ##
## Each version of the Nbody simulation implements various optimizations     ##
## that can be compiled with the following compilers (provided that they are ##
## installed on your system):                                                ##
##   - GCC                                                                   ##
##   - Clang (LLVM and/or AOCC)                                              ##
##   - ICC                                                                   ##
##   - ICX                                                                   ##
##                                                                           ##
## Each of the programs can be compiled in both FP32 and FP64 precision by   ##
## adding the following variable definition in their make command (FP32 is   ##
## the default):                                                             ##
##   make PRECISION=-DFP64													 ##
## ------------------------------------------------------------------------- ##

## MAKEFILE MACROS
GCC=gcc
CLANG=clang
ICC=icc
ICX=icx

GCC_OFLAGS=-Wall -Wextra -g3 -march=native -mtune=native -mavx2 -Ofast -finline-functions -funroll-loops -fpeel-loops -floop-interchange -ftree-vectorize -ftree-loop-vectorize -flto
CLANG_OFLAGS=-Wall -Wextra -g3 -march=native -mtune=native -mavx2 -Ofast -finline-functions -funroll-loops
ICC_OFLAGS=-w2 -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto
ICX_OFLAGS=-Wall -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto

GCC_CLANG_LDFLAGS=-lm -fopenmp
ICC_ICX_LDFLAGS=-qmkl -qopenmp


## MAKEFILE RULES
all: base seq par

seq: base soa seq_soap seq_avx seq_avx512

par: par_soap par_avx par_avx512

gcc: target/gcc/nbody0 target/gcc/nbody_soa target/gcc/nbody_seq_soap target/gcc/nbody_par_soap target/gcc/nbody_seq_avx target/gcc/nbody_par_avx target/gcc/nbody_seq_avx512 target/gcc/nbody_par_avx512

clang: target/clang/nbody0  target/clang/nbody_soa target/clang/nbody_seq_soap target/clang/nbody_par_soap target/clang/nbody_seq_avx target/clang/nbody_par_avx target/clang/nbody_seq_avx512 target/clang/nbody_par_avx512

icc: target/icc/nbody0  target/icc/nbody_soa target/icc/nbody_seq_soap target/icc/nbody_par_soap target/icc/nbody_seq_avx target/icc/nbody_par_avx target/icc/nbody_seq_avx512 target/icc/nbody_par_avx512

icx: target/icx/nbody0  target/icx/nbody_soa target/icx/nbody_seq_soap target/icx/nbody_par_soap target/icx/nbody_seq_avx target/icx/nbody_par_avx target/icx/nbody_seq_avx512 target/icx/nbody_par_avx512

base: target/gcc/nbody0 target/clang/nbody0 target/icc/nbody0 target/icx/nbody0

soa: target/gcc/nbody_soa target/clang/nbody_soa target/icc/nbody_soa target/icx/nbody_soa

seq_soap: target/gcc/nbody_seq_soap target/clang/nbody_seq_soap target/icc/nbody_seq_soap target/icx/nbody_seq_soap

par_soap: target/gcc/nbody_par_soap target/clang/nbody_par_soap target/icc/nbody_par_soap target/icx/nbody_par_soap

seq_avx: target/gcc/nbody_seq_avx target/clang/nbody_seq_avx target/icc/nbody_seq_avx target/icx/nbody_seq_avx

par_avx: target/gcc/nbody_par_avx target/clang/nbody_par_avx target/icc/nbody_par_avx target/icx/nbody_par_avx

seq_avx512: target/gcc/nbody_seq_avx512 target/clang/nbody_seq_avx512 target/icc/nbody_seq_avx512 target/icx/nbody_seq_avx512

par_avx512: target/gcc/nbody_par_avx512 target/clang/nbody_par_avx512 target/icc/nbody_par_avx512 target/icx/nbody_par_avx512


## COMPILATION TARGETS
ref: src/nbody0.c src/utils.c
	@mkdir -p target/
	$(CLANG) $(PRECISION) -g3 -ffp-model=precise -march=native -mtune=native -O2 -mno-avx $^ -o target/$@ $(GCC_CLANG_LDFLAGS)
	target/ref $(PARAMS) -d -o data/ref.dat

target/gcc/nbody%: src/nbody%.c src/utils.c
	@mkdir -p target/gcc/
	$(GCC) $(PRECISION) $(GCC_OFLAGS) $^ -o $@ $(GCC_CLANG_LDFLAGS)

target/clang/nbody%: src/nbody%.c src/utils.c
	@mkdir -p target/clang/
	$(CLANG) $(PRECISION) $(CLANG_OFLAGS) $^ -o $@ $(GCC_CLANG_LDFLAGS)

target/icc/nbody%: src/nbody%.c src/utils.c
	@mkdir -p target/icc/
	$(ICC) $(PRECISION) $(ICC_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)
	
target/icx/nbody%: src/nbody%.c src/utils.c
	@mkdir -p target/icx/
	$(ICX) $(PRECISION) $(ICX_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)


## UTILS
clean:
	@rm -Rf target/ opt_reports/ maqao/**/* data/**/*
