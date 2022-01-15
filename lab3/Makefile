# Macros
GC=gcc
LC=clang
IC=icc
XC=icx

GFLAGS=-march=native -mtune=native -mavx2 -Ofast -finline-functions -funroll-loops -fpeel-loops -floop-interchange -ftree-vectorize -ftree-loop-vectorize -flto
LFLAGS=-march=native -mtune=native -mavx2 -Ofast -finline-functions -funroll-loops
IFLAGS=-w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -flto

GLLD=-lm -fopenmp
IXLD=-qmkl -qopenmp

# Rules
all: 0 1 2 3 _seq_avx _par_avx _seq_avx512 _par_avx512

seq: 0 1 2 _seq_avx _seq_avx512

par: 3 _par_avx _par_avx512

0: target/gcc/nbody0.g target/clang/nbody0.l target/icc/nbody0.i target/icx/nbody0.x

1: target/gcc/nbody1.g target/clang/nbody1.l target/icc/nbody1.i target/icx/nbody1.x

2: target/gcc/nbody2.g target/clang/nbody2.l target/icc/nbody2.i target/icx/nbody2.x

3: target/gcc/nbody3.g target/clang/nbody3.l target/icc/nbody3.i target/icx/nbody3.x

_seq_avx: target/gcc/nbody_seq_avx.g target/clang/nbody_seq_avx.l target/icc/nbody_seq_avx.i target/icx/nbody_seq_avx.x

_par_avx: target/gcc/nbody_par_avx.g target/clang/nbody_par_avx.l target/icc/nbody_par_avx.i target/icx/nbody_par_avx.x

_seq_avx512: target/gcc/nbody_seq_avx512.g target/clang/nbody_seq_avx512.l target/icc/nbody_seq_avx512.i target/icx/nbody_seq_avx512.x

_par_avx512: target/gcc/nbody_par_avx512.g target/clang/nbody_par_avx512.l target/icc/nbody_par_avx512.i target/icx/nbody_par_avx512.x
	
# Compilation
target/gcc/nbody%.g: src/nbody%.c
	@mkdir -p target/gcc/
	$(GC) $(GFLAGS) $< -o $@ $(GLLD)

target/clang/nbody%.l: src/nbody%.c
	@mkdir -p target/clang/
	$(LC) $(LFLAGS) $< -o $@ $(GLLD)

target/icc/nbody%.i: src/nbody%.c
	@mkdir -p target/icc/
	$(IC) $(IFLAGS) $< -o $@ $(IXLD)
	
target/icx/nbody%.x: src/nbody%.c
	@mkdir -p target/icx/
	$(XC) $(IFLAGS) $< -o $@ $(IXLD)

clean:
	@rm -Rf target/ opt_reports/ maqao/*
