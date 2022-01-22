#! /bin/bash

make seq -j
mkdir -p output/gcc output/icc output/icx

iter=5
warmup=2

for i in $(seq 16384 16384 65536); do
    # GCC sequential
    numactl -C 7 target/gcc/nbody0 -d -b -n $i -i $iter -w $warmup -o output/gcc/base_$i.dat
    numactl -C 7 target/gcc/nbody1_soa -d -b -n $i -i $iter -w $warmup -o output/gcc/soa_$i.dat
    numactl -C 7 target/gcc/nbody2_ssoai -d -b -n $i -i $iter -w $warmup -o output/gcc/ssoai_$i.dat
    numactl -C 7 target/gcc/nbody3_savx -d -b -n $i -i $iter -w $warmup -o output/gcc/savx_$i.dat
    numactl -C 7 target/gcc/nbody4_savx512 -d -b -n $i -i $iter -w $warmup -o output/gcc/savx512_$i.dat

    # ICC sequential
    numactl -C 7 target/icc/nbody0 -d -b -n $i -i $iter -w $warmup -o output/icc/base_$i.dat
    numactl -C 7 target/icc/nbody1_soa -d -b -n $i -i $iter -w $warmup -o output/icc/soa_$i.dat
    numactl -C 7 target/icc/nbody2_ssoai -d -b -n $i -i $iter -w $warmup -o output/icc/ssoai_$i.dat
    numactl -C 7 target/icc/nbody3_savx -d -b -n $i -i $iter -w $warmup -o output/icc/savx_$i.dat
    numactl -C 7 target/icc/nbody4_savx512 -d -b -n $i -i $iter -w $warmup -o output/icc/savx512_$i.dat

    # ICX sequential
    numactl -C 7 target/icx/nbody0 -d -b -n $i -i $iter -w $warmup -o output/icx/base_$i.dat
    numactl -C 7 target/icx/nbody1_soa -d -b -n $i -i $iter -w $warmup -o output/icx/soa_$i.dat
    numactl -C 7 target/icx/nbody2_ssoai -d -b -n $i -i $iter -w $warmup -o output/icx/ssoai_$i.dat
    numactl -C 7 target/icx/nbody3_savx -d -b -n $i -i $iter -w $warmup -o output/icx/savx_$i.dat
    numactl -C 7 target/icx/nbody4_savx512 -d -b -n $i -i $iter -w $warmup -o output/icx/savx512_$i.dat
done
exit 0
