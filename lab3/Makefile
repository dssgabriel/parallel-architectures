all: target/nbody0.g target/nbody0.i target/nbody1.g target/nbody1.i target/nbody2.g target/nbody2.i target/nbody2.l target/nbody3.g target/nbody3.i 

target/nbody0.g: src/nbody0.c
	@mkdir -p target opt_reports
	gcc -march=native -mtune=native -mavx2 -Ofast -fopt-info-all=opt_reports/nbody0.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody0.i: src/nbody0.c
	@mkdir -p target opt_reports
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody1.g: src/nbody1.c
	@mkdir -p target opt_reports
	gcc -march=native -mtune=native -mavx2 -Ofast -fopt-info-all=opt_reports/nbody1.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody1.i: src/nbody1.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -axTIGERLAKE -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp

target/nbody2.g: src/nbody2.c
	@mkdir -p target opt_reports
	gcc -Wall -g -g3 -march=tigerlake -mtune=tigerlake -mavx512f -mavx512pf -mavx512er -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mavx512ifma -mavx512vbmi -mavx512vbmi2 -mavx512bf16 -mavx512bitalg -mavx512vpopcntdq -mavx5124fmaps  -mavx512vnni -mavx5124vnniw -mavx512vp2intersect -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody2.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody2.i: src/nbody2.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -axTIGERLAKE -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
target/nbody2.l: src/nbody2.c
	@mkdir -p target opt_reports
	clang -Wall -g -g3 -march=tigerlake -mtune=tigerlake -mavx2 -Ofast -fopenmp -funroll-loops -finline-functions -ftree-vectorize $< -o $@ -lm
	
target/nbody3.g: src/nbody3.c
	@mkdir -p target opt_reports
	gcc -Wall -g -g3 -march=tigerlake -mtune=tigerlake -mavx512f -mavx512pf -mavx512er -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mavx512ifma -mavx512vbmi -mavx512vbmi2 -mavx512bf16 -mavx512bitalg -mavx512vpopcntdq -mavx5124fmaps  -mavx512vnni -mavx5124vnniw -mavx512vp2intersect -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody3.i: src/nbody3.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -axTIGERLAKE -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody4.g: src/nbody4.c
	@mkdir -p target opt_reports
	gcc -Wall -g3 -march=tigerlake -mtune=tigerlake -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize $< -o $@ -lm -fopenmp
target/nbody4.i: src/nbody4.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -axTIGERLAKE -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
target/nbody5.g: src/nbody5.c
	@mkdir -p target opt_reports
	gcc -Wall -g -g3 -march=tigerlake -mtune=tigerlake -mavx512f -mavx512pf -mavx512er -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mavx512ifma -mavx512vbmi -mavx512vbmi2 -mavx512bf16 -mavx512bitalg -mavx512vpopcntdq -mavx5124fmaps  -mavx512vnni -mavx5124vnniw -mavx512vp2intersect -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -ftree-loop-vectorize -fopt-info-all=opt_reports/nbody3.gcc.optrpt $< -o $@ -lm -fopenmp
target/nbody5.i: src/nbody5.c
	@mkdir -p target opt_reports
	icc -w3 -g3 -xHost -axTIGERLAKE -Ofast -funroll-loops -finline-functions -qopt-report $< -o $@ -qmkl -qopenmp
	  
clean:
	@rm -Rf target/ opt_reports/ maqao/*

