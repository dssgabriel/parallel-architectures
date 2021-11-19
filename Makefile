CC = gcc
CFLG = -Wall -Wextra -g -mno-avx
OFLG = -march=native -mtune=native -Ofast -finline-functions -funroll-loops -fpeel-loops -ftree-vectorize -ftree-loop-vectorize

all: target/nbody0

target/nbody%: src/nbody%.c
	@mkdir -p target
	$(CC) $(CFLG) $(OFLG) $< -o $@ -lm -lSDL2 

clean:
	rm -rf target/
