## ------------------------------------------------------------------------- ##
## This Makefile builds various versions of a 3D Nbody simulation.           ##
##                                                                           ##
## Each version of the Nbody simulation implements various optimizations     ##
## that can be compiled with the following compilers (provided that they are ##
## installed on your system):                                                ##
##   - GCC                                                                   ##
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
ICC=icc
ICX=icx

GCC_OFLAGS=-Wall -Wextra -g3 -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -floop-interchange -ftree-vectorize -ftree-loop-vectorize -flto
ICC_OFLAGS=-w2 -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto -qopt-report=5 -qopt-report-file=$@.optrpt
ICX_OFLAGS=-Wall -g3 -xHost -mavx -Ofast -funroll-loops -finline-functions -flto

GCC_LDFLAGS=-lm -fopenmp
ICC_ICX_LDFLAGS=-qmkl -qopenmp

TGCC=target/gcc
TICC=target/icc
TICX=target/icx


## MAKEFILE RULES
all: base seq par

seq: base soa ssoai savx savx512

par: psoai pavx pavx512

par_block: pbsoai pbavx pbavx512

gcc: $(TGCC)/nbody0 $(TGCC)/nbody1_soa $(TGCC)/nbody2_ssoai $(TGCC)/nbody3_savx $(TGCC)/nbody4_savx512 $(TGCC)/nbody5_psoai $(TGCC)/nbody6_pavx $(TGCC)/nbody7_pavx512 $(TGCC)/nbody8_pbsoai $(TGCC)/nbody9_pbavx $(TGCC)/nbodyA_pbavx512

icc: $(TICC)/nbody0 $(TICC)/nbody1_soa $(TICC)/nbody2_ssoai $(TICC)/nbody3_savx $(TICC)/nbody4_savx512 $(TICC)/nbody5_psoai $(TICC)/nbody6_pavx $(TICC)/nbody7_pavx512 $(TICC)/nbody8_pbsoai $(TICC)/nbody9_pbavx $(TICC)/nbodyA_pbavx512

icx: $(TICX)/nbody0 $(TICX)/nbody1_soa $(TICX)/nbody2_ssoai $(TICX)/nbody3_savx $(TICX)/nbody4_savx512 $(TICX)/nbody5_psoai $(TICX)/nbody6_pavx $(TICX)/nbody7_pavx512 $(TICX)/nbody8_pbsoai $(TICX)/nbody9_pbavx $(TICX)/nbodyA_pbavx512

base: $(TGCC)/nbody0 $(TICC)/nbody0 $(TICX)/nbody0

soa: $(TGCC)/nbody1_soa $(TICC)/nbody1_soa $(TICX)/nbody1_soa

ssoai: $(TGCC)/nbody2_ssoai $(TICC)/nbody2_ssoai $(TICX)/nbody2_ssoai

savx: $(TGCC)/nbody3_savx $(TICC)/nbody3_savx $(TICX)/nbody3_savx

savx512: $(TGCC)/nbody4_savx512 $(TICC)/nbody4_savx512 $(TICX)/nbody4_savx512

psoai: $(TGCC)/nbody5_psoai $(TICC)/nbody5_psoai $(TICX)/nbody5_psoai

pavx: $(TGCC)/nbody6_pavx $(TICC)/nbody6_pavx $(TICX)/nbody6_pavx

pavx512: $(TGCC)/nbody7_pavx512 $(TICC)/nbody7_pavx512 $(TICX)/nbody7_pavx512

pbsoai: $(TGCC)/nbody8_pbsoai $(TICC)/nbody8_pbsoai $(TICX)/nbody8_pbsoai

pbavx: $(TGCC)/nbody9_pbavx $(TICC)/nbody9_pbavx $(TICX)/nbody9_pbavx

pbavx512: $(TGCC)/nbodyA_pbavx512 $(TICC)/nbodyA_pbavx512 $(TICX)/nbodyA_pbavx512


## COMPILATION TARGETS
ref: src/nbody0.c src/utils.c
	@mkdir -p target/
	$(GCC) $(PRECISION) -g3 -march=native -mtune=native -O1 -mno-avx $^ -o target/$@ $(GCC_LDFLAGS)
	target/ref -d -i 1 -w 0 -o data/ref.dat

$(TGCC)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TGCC)
	$(GCC) $(PRECISION) $(GCC_OFLAGS) $^ -o $@ $(GCC_LDFLAGS)

$(TICC)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TICC)
	$(ICC) $(PRECISION) $(ICC_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)
	
$(TICX)/nbody%: src/nbody%.c src/utils.c
	@mkdir -p $(TICX)
	$(ICX) $(PRECISION) $(ICX_OFLAGS) $^ -o $@ $(ICC_ICX_LDFLAGS)


## UTILS
clean:
	@rm -Rf target/ opt_reports/ maqao/**/* data/**/* data/ref.dat bench/precision/**/* bench/perf/**/* 
