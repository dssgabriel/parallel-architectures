#1 /bin/bash

make seq
make ref PARAMS="-i1 -w0"
target/gcc/nbody0 -i1 -w0 -0data/gcc/base.dat
target/gcc/nbody_soa -i1 -w0 -0data/gcc/soa.dat
target/gcc/nbody_seq_soap -i1 -w0 -0data/gcc/ssoa+.dat
target/gcc/nbody_seq_avx -i1 -w0 -0data/gcc/savx.dat
target/gcc/nbody_seq_avx512 -i1 -w0 -0data/gcc/savx512.dat

target/icc/nbody0 -i1 -w0 -0data/icc/base.dat
target/icc/nbody_soa -i1 -w0 -0data/icc/soa.dat
target/icc/nbody_seq_soap -i1 -w0 -0data/icc/ssoa+.dat
target/icc/nbody_seq_avx -i1 -w0 -0data/icc/savx.dat
target/icc/nbody_seq_avx512 -i1 -w0 -0data/icc/savx512.dat

target/icx/nbody0 -i1 -w0 -0data/icx/base.dat
target/icx/nbody_soa -i1 -w0 -0data/icx/soa.dat
target/icx/nbody_seq_soap -i1 -w0 -0data/icx/ssoa+.dat
target/icx/nbody_seq_avx -i1 -w0 -0data/icx/savx.dat
target/icx/nbody_seq_avx512 -i1 -w0 -0data/icx/savx512.dat

python precision_checker.py data/gcc/base.dat
python precision_checker.py data/gcc/soa.dat
python precision_checker.py data/gcc/ssoa+.dat
python precision_checker.py data/gcc/savx.dat
python precision_checker.py data/gcc/savx512.dat

python precision_checker.py data/icc/base.dat
python precision_checker.py data/icc/soa.dat
python precision_checker.py data/icc/ssoa+.dat
python precision_checker.py data/icc/savx.dat
python precision_checker.py data/icc/savx512.dat

python precision_checker.py data/icx/base.dat
python precision_checker.py data/icx/soa.dat
python precision_checker.py data/icx/ssoa+.dat
python precision_checker.py data/icx/savx.dat
python precision_checker.py data/icx/savx512.dat
