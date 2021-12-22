all: target/nbody0.g target/nbody0.i target/nbody1.g target/nbody1.i target/nbody2.g target/nbody2.i target/nbody2.l target/nbody3.g target/nbody3.i target/nbody_seq_avx.g target/nbody_seq_avx.i target/nbody_par_avx.g target/nbody_par_avx.i target/nbody_seq_avx512.g target/nbody_seq_avx512.i target/nbody_par_avx512.g target/nbody_par_avx512.i
	@mv target/*.optrpt opt_reports/

target/nbody0.g: src/nbody0.c
	@mkdir -p target opt_reports
	gcc -march=native -mtune=native -mavx2 -Ofast -fopt-info-all=opt_reports/nbody0.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody0.i: src/nbody0.c
	@mkdir -p target opt_reports
	icc -xHost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody1.g: src/nbody1.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx2 -Ofast -fopt-info-all=opt_reports/nbody1.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody1.i: src/nbody1.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody2.g: src/nbody2.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody2.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody2.i: src/nbody2.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
target/nbody2.l: src/nbody2.c
	@mkdir -p target opt_reports
	clang -Wall -g -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -fopenmp -funroll-loops -finline-functions -ftree-vectorize $< -o $@ -lm

target/nbody3.g: src/nbody3.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody3.i: src/nbody3.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
  
target/nbody_seq_avx.g: src/nbody_seq_avx.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody_seq_avx.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody_seq_avx.i: src/nbody_seq_avx.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
  
target/nbody_par_avx.g: src/nbody_par_avx.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody_par_avx.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody_par_avx.i: src/nbody_par_avx.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
  
target/nbody_seq_avx512.g: src/nbody_seq_avx512.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody_seq_avx512.i: src/nbody_seq_avx512.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody_par_avx512.g: src/nbody_par_avx512.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody_par_avx512.i: src/nbody_par_avx512.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	@rm -Rf target/ opt_reports/ maqao/*
