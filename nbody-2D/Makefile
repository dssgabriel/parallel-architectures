CC := gcc
CFLAGS := -Wall -Wextra -g
OFLAGS := -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -ftree-vectorize -ftree-loop-vectorize
LDFLAGS := -lm -lSDL2

all: build

build: target/nbody0 target/nbody1 target/nbody2 target/nbody3 target/nbody4 target/nbody5

bench: build
	@mkdir -p bench/
	taskset -c 7 target/nbody0 > bench/out_c.dat
	taskset -c 7 target/nbody1 > bench/out_sd.dat
	taskset -c 7 target/nbody2 > bench/out_pd.dat
	taskset -c 7 target/nbody3 > bench/out_soa.dat
	taskset -c 7 target/nbody5 > bench/out_intrin.dat
	@gnuplot scripts/plot.gp

target/nbody%: src/nbody%.c
	@mkdir -p target
	$(CC) $(CFLAGS) $(OFLAGS) $< -o $@ $(LDFLAGS)

clean:
	@rm -rf target/ bench/
