#! /bin/bash
# ---------------------------------------------------------------------------- #
# This script needs to be modified when run on the KNL/KNM cluster.            #
# ---------------------------------------------------------------------------- #

mkdir -p output/gcc output/icc output/icx
printf "Compiling parallel versions of N-body... "
make par -j >/dev/null 2>&1
printf "\033[1;32mdone\033[0m\n\n"

start=$(date +%s.%N)

iter=60
warmup=20
# Small values
for i in $(seq 65536 65536 262144); do
    printf "\rBenchmarking with %d particles..." $i
    # SoA improved
    target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat >/dev/null 2>&1
    target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat >/dev/null 2>&1
    target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat >/dev/null 2>&1

    # AVX
    target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat >/dev/null 2>&1
    target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat >/dev/null 2>&1
    target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat >/dev/null 2>&1

    #AVX512
    target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat >/dev/null 2>&1
    target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat >/dev/null 2>&1
    target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat >/dev/null 2>&1
done

iter=18
warmup=6
# Medium values
for i in $(seq 524288 262144 1048576); do
    printf "\rBenchmarking with %d particles..." $i
    # SoA improved
    target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat >/dev/null 2>&1
    target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat >/dev/null 2>&1
    target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat >/dev/null 2>&1

    # AVX
    target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat >/dev/null 2>&1
    target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat >/dev/null 2>&1
    target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat >/dev/null 2>&1

    #AVX512
    target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat >/dev/null 2>&1
    target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat >/dev/null 2>&1
    target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat >/dev/null 2>&1
done

iter=6
warmup=2
# Medium values
for i in $(seq 1572854 524288 2097152); do
    printf "\rBenchmarking with %d particles..." $i
    # SoA improved
    target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat >/dev/null 2>&1
    target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat >/dev/null 2>&1
    target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat >/dev/null 2>&1

    # AVX
    target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat >/dev/null 2>&1
    target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat >/dev/null 2>&1
    target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat >/dev/null 2>&1

    #AVX512
    target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat >/dev/null 2>&1
    target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat >/dev/null 2>&1
    target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat >/dev/null 2>&1
done
end=$(date +%s.%N)

diff=$(echo "scale=3; $end - $start" | bc -l)
printf "\n\nFinished in %.2f seconds\n" $diff

exit 0
