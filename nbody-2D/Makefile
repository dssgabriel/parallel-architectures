CC = gcc
CFLG = -Wall -Wextra -g -fopenmp
OFLG = -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -ftree-vectorize -ftree-loop-vectorize

all: build

build: target/nbody0 target/nbody1 target/nbody2 target/nbody3 target/nbody4

bench: build
	@mkdir -p bench/
	taskset -c 7 target/nbody0 > bench/out_c.dat
	taskset -c 7 target/nbody1 > bench/out_sd.dat
	taskset -c 7 target/nbody2 > bench/out_pd.dat
	taskset -c 7 target/nbody3 > bench/out_soa.dat
	taskset -c 7 target/nbody4_icc > bench/out_cheat.dat
	@gnuplot plot.gp

target/nbody%: src/nbody%.c
	@mkdir -p target
	$(CC) $(CFLG) $(OFLG) $< -o $@ -lm -lSDL2

target/nbody4_icc: src/nbody4.c
	@mkdir -p target
	icc -xHost -Ofast $< -o $@ -qmkl -qopenmp -lSDL2

clean:
	rm -rf target/ bench/
