all: target/nbody0.g target/nbody0.i target/nbody1.g target/nbody1.i target/nbody2.g target/nbody2.i target/nbody2.l target/nbody3.g target/nbody3.i target/nbody4.g target/nbody4.i target/nbody5.g target/nbody5.i target/nbody6.g target/nbody6.i target/nbody7.g target/nbody7.i target/nbody8.g target/nbody8.i target/nbody9.g target/nbody9.i target/nbody10.g target/nbody10.i target/nbody11.g target/nbody11.i

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
	  
target/nbody4.g: src/nbody4.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody4.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody4.i: src/nbody4.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody5.g: src/nbody5.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody5.i: src/nbody5.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody6.g: src/nbody6.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody6.i: src/nbody6.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody7.g: src/nbody7.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody7.i: src/nbody7.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	
target/nbody8.g: src/nbody8.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody8.i: src/nbody8.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody9.g: src/nbody9.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody9.i: src/nbody9.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody10.g: src/nbody10.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody10.i: src/nbody10.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody11.g: src/nbody11.c
	@mkdir -p target opt_reports
	gcc -Wall -Wextra -g3 -march=native -mtune=native -mavx512f -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody5.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody11.i: src/nbody11.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	
clean:
	@rm -Rf target/ opt_reports/ maqao/*
