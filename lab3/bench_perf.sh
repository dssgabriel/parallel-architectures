#! /bin/bash

make -j
mkdir -p output/gcc output/icc output/icx

step=16384
iter=60
warmup=15

for i in $(seq 16384 $step 2097152); do
    if [[ $i -lt 65536 ]]; then
        # GCC sequential runs
        numactl -C 71 target/gcc/nbody0 -b -n $i -o output/gcc/base_$i.dat
        numactl -C 71 target/gcc/nbody1_soa -b -n $i -o output/gcc/soa_$i.dat
        numactl -C 71 target/gcc/nbody2_ssoai -b -n $i -o output/gcc/ssoai_$i.dat
        numactl -C 71 target/gcc/nbody3_savx -b -n $i -o output/gcc/savx_$i.dat
        numactl -C 71 target/gcc/nbody4_savx512 -b -n $i -o output/gcc/savx512_$i.dat
        # GCC parallel runs
        target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat
        target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat
        target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat

        # ICC sequential runs
        numactl -C 71 target/icc/nbody0 -b -n $i -o output/icc/base_$i.dat
        numactl -C 71 target/icc/nbody1_soa -b -n $i -o output/icc/soa_$i.dat
        numactl -C 71 target/icc/nbody2_ssoai -b -n $i -o output/icc/ssoai_$i.dat
        numactl -C 71 target/icc/nbody3_savx -b -n $i -o output/icc/savx_$i.dat
        numactl -C 71 target/icc/nbody4_savx512 -b -n $i -o output/icc/savx512_$i.dat
        # ICC parallel runs
        target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat
        target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat
        target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat

        # ICX sequential runs
        numactl -C 7 target/icx/nbody0 -b -n $i -o output/icx/base_$i.dat
        numactl -C 7 target/icx/nbody1_soa -b -n $i -o output/icx/soa_$i.dat
        numactl -C 7 target/icx/nbody2_ssoai -b -n $i -o output/icx/ssoai_$i.dat
        numactl -C 7 target/icx/nbody3_savx -b -n $i -o output/icx/savx_$i.dat
        numactl -C 7 target/icx/nbody4_savx512 -b -n $i -o output/icx/savx512_$i.dat
        # ICX parallel runs
        target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat
        target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat
        target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat
    elif [[ $i -lt 262144 ]] && [[ $i -ge 65536 ]]; then
        step=32768
        iter=20
        warmup=5

        # GCC parallel runs
        target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat
        target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat
        target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat

        # ICC parallel runs
        target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat
        target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat
        target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat

        # ICX parallel runs
        target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat
        target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat
        target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat
    elif [[ $i -lt 1048576 ]] && [[ $i -ge 262144 ]]; then
        step=131072
        iter=10
        warmup=3

        # GCC parallel runs
        target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat
        target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat
        target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat

        # ICC parallel runs
        target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat
        target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat
        target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat

        # ICX parallel runs
        target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat
        target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat
        target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat
    else
        step=1048576
        iter=3
        warmup=1

        # GCC parallel runs
        target/gcc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/gcc/psoai_$i.dat
        target/gcc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/gcc/pavx_$i.dat
        target/gcc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/gcc/pavx512_$i.dat

        # ICC parallel runs
        target/icc/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icc/psoai_$i.dat
        target/icc/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icc/pavx_$i.dat
        target/icc/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icc/pavx512_$i.dat

        # ICX parallel runs
        target/icx/nbody5_psoai -b -n $i -i $iter -w $warmup -o output/icx/psoai_$i.dat
        target/icx/nbody6_pavx -b -n $i -i $iter -w $warmup -o output/icx/pavx_$i.dat
        target/icx/nbody7_pavx512 -b -n $i -i $iter -w $warmup -o output/icx/pavx512_$i.dat
    fi
done

exit 0
