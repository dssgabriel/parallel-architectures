#! /bin/bash
# ---------------------------------------------------------------------------- #
# This script needs to be modified when run on the KNL/KNM cluster.            #
# ---------------------------------------------------------------------------- #

mkdir -p output/gcc output/icc output/icx
printf "Compiling sequential versions of N-body... "
make seq -j >/dev/null 2>&1
printf "\033[1;32mdone\033[0m\n\n"

iter=10
warmup=3

start=$(date +%s.%N)
# Small values
for i in $(seq 16384 16384 65536); do
    printf "\rBenchmarking with %d particles..." $i
    # Base
    numactl -C 7 target/gcc/nbody0 -b -n $i -i $iter -w $warmup -o output/gcc/base_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody0 -b -n $i -i $iter -w $warmup -o output/icc/base_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody0 -b -n $i -i $iter -w $warmup -o output/icx/base_$i.dat >/dev/null 2>&1
    
    # SoA
    numactl -C 7 target/gcc/nbody1_soa -b -n $i -i $iter -w $warmup -o output/gcc/soa_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody1_soa -b -n $i -i $iter -w $warmup -o output/icc/soa_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody1_soa -b -n $i -i $iter -w $warmup -o output/icx/soa_$i.dat >/dev/null 2>&1

    # SoA improved
    numactl -C 7 target/gcc/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/gcc/ssoai_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/icc/ssoai_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/icx/ssoai_$i.dat >/dev/null 2>&1

    # AVX
    numactl -C 7 target/gcc/nbody3_savx -b -n $i -i $iter -w $warmup -o output/gcc/savx_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody3_savx -b -n $i -i $iter -w $warmup -o output/icc/savx_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody3_savx -b -n $i -i $iter -w $warmup -o output/icx/savx_$i.dat >/dev/null 2>&1

    #AVX512
    numactl -C 7 target/gcc/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/gcc/savx512_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/icc/savx512_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/icx/savx512_$i.dat >/dev/null 2>&1
done

iter=5
warmup=2
# Big values
for i in $(seq 98304 32768 131072); do
    # Base
    numactl -C 7 target/gcc/nbody0 -b -n $i -i $iter -w $warmup -o output/gcc/base_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody0 -b -n $i -i $iter -w $warmup -o output/icc/base_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody0 -b -n $i -i $iter -w $warmup -o output/icx/base_$i.dat >/dev/null 2>&1
    
    # SoA
    numactl -C 7 target/gcc/nbody1_soa -b -n $i -i $iter -w $warmup -o output/gcc/soa_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody1_soa -b -n $i -i $iter -w $warmup -o output/icc/soa_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody1_soa -b -n $i -i $iter -w $warmup -o output/icx/soa_$i.dat >/dev/null 2>&1

    # SoA improved
    numactl -C 7 target/gcc/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/gcc/ssoai_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/icc/ssoai_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody2_ssoai -b -n $i -i $iter -w $warmup -o output/icx/ssoai_$i.dat >/dev/null 2>&1

    # AVX
    numactl -C 7 target/gcc/nbody3_savx -b -n $i -i $iter -w $warmup -o output/gcc/savx_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody3_savx -b -n $i -i $iter -w $warmup -o output/icc/savx_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody3_savx -b -n $i -i $iter -w $warmup -o output/icx/savx_$i.dat >/dev/null 2>&1

    #AVX512
    numactl -C 7 target/gcc/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/gcc/savx512_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icc/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/icc/savx512_$i.dat >/dev/null 2>&1
    numactl -C 7 target/icx/nbody4_savx512 -b -n $i -i $iter -w $warmup -o output/icx/savx512_$i.dat >/dev/null 2>&1
done
end=$(date +%s.%N)

diff=$(echo "scale=3; $end - $start" | bc -l)
printf "\n\nFinished in %.2f seconds\n" $diff

exit 0
