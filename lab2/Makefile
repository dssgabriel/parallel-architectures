CC = gcc
CFLG = -Wall -Wextra -g -mno-avx
OFLG = -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -ftree-vectorize -ftree-loop-vectorize

all: build

build: target/nbody0 target/nbody1 target/nbody2 target/nbody3

bench: build
	@mkdir -p bench/
	taskset -c 7 target/nbody0 > bench/out_c.dat
	taskset -c 7 target/nbody1 > bench/out_sd.dat
	taskset -c 7 target/nbody2 > bench/out_pd.dat
	taskset -c 7 target/nbody3 > bench/out_soa.dat
	@gnuplot plot.gp

target/nbody%: src/nbody%.c
	@mkdir -p target
	$(CC) $(CFLG) $(OFLG) $< -o $@ -lm -lSDL2 

clean:
	rm -rf target/ bench/
