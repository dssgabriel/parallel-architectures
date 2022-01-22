#! /bin/bash

make target/icx/nbodyA_pbavx512

target/icx/nbodyA_pbavx512 -d -n 524288 -i3 -w0 -b -o data/icx/pbavx512.dat
target/icx/nbodyA_pbavx512 -d -n 1048576 -i3 -w0 -b -o data/icx/pbavx512.dat
target/icx/nbodyA_pbavx512 -d -n 2097152 -i3 -w0 -b -o data/icx/pbavx512.dat
target/icx/nbodyA_pbavx512 -d -n 4194304 -i3 -w0 -b -o data/icx/pbavx512.dat
target/icx/nbodyA_pbavx512 -d -n 6291456 -i3 -w0 -b -o data/icx/pbavx512.dat
