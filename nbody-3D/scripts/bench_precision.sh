#! /bin/bash

make -j
make ref

# Bench GCC sequential
target/gcc/nbody0 -c -o data/gcc/base.dat
target/gcc/nbody1_soa -c -o data/gcc/soa.dat
target/gcc/nbody2_ssoai -c -o data/gcc/ssoai.dat
target/gcc/nbody3_savx -c -o data/gcc/savx.dat
target/gcc/nbody4_savx512 -c -o data/gcc/savx512.dat
# Bench GCC parallel
target/gcc/nbody5_psoai -c -o data/gcc/psoai.dat
target/gcc/nbody6_pavx -c -o data/gcc/pavx.dat
target/gcc/nbody7_pavx512 -c -o data/gcc/pavx512.dat

# Bench ICC sequential
target/icc/nbody0 -c -o data/icc/base.dat
target/icc/nbody1_soa -c -o data/icc/soa.dat
target/icc/nbody2_ssoai -c -o data/icc/ssoai.dat
target/icc/nbody3_savx -c -o data/icc/savx.dat
target/icc/nbody4_savx512 -c -o data/icc/savx512.dat
# Bench ICC parallel
target/icc/nbody5_psoai -c -o data/icc/psoai.dat
target/icc/nbody6_pavx -c -o data/icc/pavx.dat
target/icc/nbody7_pavx512 -c -o data/icc/pavx512.dat

# Bench ICX sequential
target/icx/nbody0 -c -o data/icx/base.dat
target/icx/nbody1_soa -c -o data/icx/soa.dat
target/icx/nbody2_ssoai -c -o data/icx/ssoai.dat
target/icx/nbody3_savx -c -o data/icx/savx.dat
target/icx/nbody4_savx512 -c -o data/icx/savx512.dat
# Bench ICX parallel
target/icx/nbody5_psoai -c -o data/icx/psoai.dat
target/icx/nbody6_pavx -c -o data/icx/pavx.dat
target/icx/nbody7_pavx512 -c -o data/icx/pavx512.dat

# Check GCC precision
python precision_checker.py data/gcc/base.dat > bench/precision/gcc/base
python precision_checker.py data/gcc/soa.dat > bench/precision/gcc/soa
python precision_checker.py data/gcc/ssoai.dat > bench/precision/gcc/ssoai
python precision_checker.py data/gcc/savx.dat > bench/precision/gcc/savx
python precision_checker.py data/gcc/savx512.dat > bench/precision/gcc/savx512
python precision_checker.py data/gcc/psoai.dat > bench/precision/gcc/psoai
python precision_checker.py data/gcc/pavx.dat > bench/precision/gcc/pavx
python precision_checker.py data/gcc/pavx512.dat > bench/precision/gcc/pavx512

# Check ICC precision
python precision_checker.py data/icc/base.dat > bench/precision/icc/base
python precision_checker.py data/icc/soa.dat > bench/precision/icc/soa
python precision_checker.py data/icc/ssoai.dat > bench/precision/icc/ssoai
python precision_checker.py data/icc/savx.dat > bench/precision/icc/savx
python precision_checker.py data/icc/savx512.dat > bench/precision/icc/savx512
python precision_checker.py data/icc/psoai.dat > bench/precision/icc/psoai
python precision_checker.py data/icc/pavx.dat > bench/precision/icc/pavx
python precision_checker.py data/icc/pavx512.dat > bench/precision/icc/pavx512

# Check ICX precision
python precision_checker.py data/icx/base.dat > bench/precision/icx/base
python precision_checker.py data/icx/soa.dat > bench/precision/icx/soa
python precision_checker.py data/icx/ssoai.dat > bench/precision/icx/ssoai
python precision_checker.py data/icx/savx.dat > bench/precision/icx/savx
python precision_checker.py data/icx/savx512.dat > bench/precision/icx/savx512
python precision_checker.py data/icx/psoai.dat > bench/precision/icx/psoai
python precision_checker.py data/icx/pavx.dat > bench/precision/icx/pavx
python precision_checker.py data/icx/pavx512.dat > bench/precision/icx/pavx512
