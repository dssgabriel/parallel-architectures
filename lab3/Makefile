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
## Note: The variable precision has been disabled for AVX versions of the    ##
## code (FP32 only).                                                         ##
## ------------------------------------------------------------------------- ##

## MAKEFILE MACROS
GCC=gcc
CLANG=clang
ICC=icc
ICX=icx

GCC_OFLAGS=-Wall -Wextra -g3 -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -floop-interchange -ftree-vectorize -ftree-loop-vectorize -flto
CLANG_OFLAGS=-Wall -Wextra -g3 -march=native -mtune=native -Ofast -finline-functions -funroll-loops
ICC_OFLAGS=-w2 -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto -qopt-report=5 -qopt-report-file=$@.optrpt
ICX_OFLAGS=-Wall -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto

GCC_CLANG_LDFLAGS=-lm -fopenmp
ICC_ICX_LDFLAGS=-qmkl -qopenmp

TGCC=target/gcc
TCLG=target/clang
TICC=target/icc
TICX=target/icx


## MAKEFILE RULES
all: base seq par

seq: base soa ssoa+ savx savx512

par: psoa+ pavx pavx512

clang: $(TCLG)/nbody0 $(TCLG)/nbody1_soa $(TCLG)/nbody2_ssoa+ $(TCLG)/nbody3_savx $(TCLG)/nbody4_savx512 $(TCLG)/nbody5_psoa+ $(TCLG)/nbody6_pavx $(TCLG)/nbody7_pavx512

gcc: $(TGCC)/nbody0 $(TGCC)/nbody1_soa $(TGCC)/nbody2_ssoa+ $(TGCC)/nbody3_savx $(TGCC)/nbody4_savx512 $(TGCC)/nbody5_psoa+ $(TGCC)/nbody6_pavx $(TGCC)/nbody7_pavx512

icc: $(TICC)/nbody0 $(TICC)/nbody1_soa $(TICC)/nbody2_ssoa+ $(TICC)/nbody3_savx $(TICC)/nbody4_savx512 $(TICC)/nbody5_psoa+ $(TICC)/nbody6_pavx $(TICC)/nbody7_pavx512

icx: $(TICX)/nbody0 $(TICX)/nbody1_soa $(TICX)/nbody2_ssoa+ $(TICX)/nbody3_savx $(TICX)/nbody4_savx512 $(TICX)/nbody5_psoa+ $(TICX)/nbody6_pavx $(TICX)/nbody7_pavx512

base: $(TCLG)/nbody0 $(TGCC)/nbody0 $(TICC)/nbody0 $(TICX)/nbody0

soa: $(TCLG)/nbody1_soa $(TGCC)/nbody1_soa $(TICC)/nbody1_soa $(TICX)/nbody1_soa

ssoa+: $(TCLG)/nbody2_ssoa+ $(TGCC)/nbody2_ssoa+ $(TICC)/nbody2_ssoa+ $(TICX)/nbody2_ssoa+

savx: $(TCLG)/nbody3_savx $(TGCC)/nbody3_savx $(TICC)/nbody3_savx $(TICX)/nbody3_savx

savx512: $(TCLG)/nbody4_savx512 $(TGCC)/nbody4_savx512 $(TICC)/nbody4_savx512 $(TICX)/nbody4_savx512

psoa+: $(TCLG)/nbody5_psoa+ $(TGCC)/nbody5_psoa+ $(TICC)/nbody5_psoa+ $(TICX)/nbody5_psoa+

pavx: $(TCLG)/nbody6_pavx $(TGCC)/nbody6_pavx $(TICC)/nbody6_pavx $(TICX)/nbody6_pavx

pavx512: $(TCLG)/nbody7_pavx512 $(TGCC)/nbody7_pavx512 $(TICC)/nbody7_pavx512 $(TICX)/nbody7_pavx512


## COMPILATION TARGETS
ref: src/nbody0.c src/utils.c
	@mkdir -p target/
	$(CLANG) $(PRECISION) -g3 -ffp-model=precise -march=native -mtune=native -O2 -mno-avx $^ -o target/$@ $(GCC_CLANG_LDFLAGS)
	target/ref $(PARAMS) -d -o data/ref.dat

$(TCLG)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TCLG)
	$(CLANG) $(PRECISION) $(CLANG_OFLAGS) $^ -o $@ $(GCC_CLANG_LDFLAGS)

$(TGCC)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TGCC)
	$(GCC) $(PRECISION) $(GCC_OFLAGS) $^ -o $@ $(GCC_CLANG_LDFLAGS)

$(TICC)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TICC)
	$(ICC) $(PRECISION) $(ICC_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)
	
$(TICX)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TICX)
	$(ICX) $(PRECISION) $(ICX_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)


## UTILS
clean:
	@rm -Rf target/ opt_reports/ maqao/**/* data/**/*
