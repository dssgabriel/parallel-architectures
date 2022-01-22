#! /usr/bin/python
import numpy as np
import statistics
import sys
from pathlib import Path

bt = "\033[1m"
br = "\033[1;31m"
by = "\033[1;33m"
nt = "\033[0m"

def helper():
    print("precision_checker — Python script that checks the precision of a 3D N-body simulation\n")
    print(f"{bt}USAGE:{nt}\n    ./precision_check.py <FILE>\n")
    print(f"{bt}FILE:{nt}\n    Path to the output data file containing the particles' positions.")
    print("    NOTE: It is advised to run simulations with the same number of particles and iterations"
          "as they were in the reference (located at `data/ref.dat`)")


def check_files(ref, out):
    if not ref.is_file():
        print(f"{br}error:{nt} reference data file `{ref}` does not exist")
        raise SystemExit

    if not out.is_file():
        print(f"{br}error:{nt} data file `{out}` does not exist")
        raise SystemExit


def read_files(ref, out):
    ref_cfg = ref.readline().split()
    ref.readline()

    out_cfg = out.readline().split()
    out.readline()
    
    ref_contents = ref.readlines()
    out_contents = out.readlines()

    if ref_cfg[0] != out_cfg[0]:
        print(f"{br}error:{nt} reference and comparative simulations do not have the same number of particles (reference is {ref_cfg[0]}, comparative is {out_cfg[0]})")
        raise SystemExit

    if ref_cfg[1] != out_cfg[1]:
        print(f"{by}warning:{nt} reference and comparative simulations do not have the same number of iterations (reference is {ref_cfg[1]}, comparative is {out_cfg[1]})")

    if ref_cfg[2] != out_cfg[2]:
        print(f"{by}warning:{nt} reference and comparative simulations do not have the same time step (reference is {ref_cfg[2]}, comparative is {out_cfg[2]})")

    if ref_cfg[3] != out_cfg[3]:
        print(f"{by}warning:{nt} reference and comparative simulations do not have the same precision (reference is {ref_cfg[3]}, comparative is {out_cfg[3]})")

    for idx, coords in enumerate(ref_contents):
        ref_contents[idx] = [float(p) for p in coords.split()]

    for idx, coords in enumerate(out_contents):
        out_contents[idx] = [float(p) for p in coords.split()]

    return ref_cfg, ref_contents, out_cfg, out_contents


def compute_stats(ref, out):
    err1 = 0.0
    err2 = 0.0
    max = 0.0
    min = np.inf
    for i in range(len(ref)):
        x1 = abs(ref[i][0] - out[i][0])
        y1 = abs(ref[i][1] - out[i][1])
        z1 = abs(ref[i][2] - out[i][2])
        x2 = (ref[i][0] - out[i][0]) * (ref[i][0] - out[i][0])
        y2 = (ref[i][1] - out[i][1]) * (ref[i][1] - out[i][1])
        z2 = (ref[i][2] - out[i][2]) * (ref[i][2] - out[i][2])
        sum1 = x1 + y1 + z1
        sum2 = x2 + y2 + z2
        if sum1 > max:
            max = sum1

        if sum1 < min:
            min = sum1

        err1 += sum1
        err2 += np.sqrt(sum2)
        
    err1 /= len(ref)
    err2 /= len(ref)
    return min, max, err1, err2


def results(ref_cfg, out_cfg, min, max, err1, err2):
    print("\t\t\tReference\tComparative\n"
          f"Particles\t{ref_cfg[0]}\t\t{out_cfg[0]}\n"
          f"Iterations\t{ref_cfg[1]}\t\t\t{out_cfg[1]}\n"
          f"Time step\t{ref_cfg[2]}\t{out_cfg[2]}\n"
          f"Precision\t{ref_cfg[3]}\t\t{out_cfg[3]}\n")
    print(f"Min error: {min}")
    print(f"Max error: {max}")
    print(f"Relative error (norm 1): {err1}")
    print(f"Relative error (norm 2): {err2}") 
    

def main():
    if len(sys.argv) != 2:
        print("error: wrong number of arguments\nSee help below\n")
        helper()
        raise SystemExit

    ref = Path("data/ref.dat")
    out = Path(sys.argv[1])
    check_files(ref, out)
    ref_file = open(ref, "r") 
    out_file = open(out, "r")

    ref_cfg, ref_contents, out_cfg, out_contents = read_files(ref_file, out_file)
    min, max, err1, err2 = compute_stats(ref_contents, out_contents)
    results(ref_cfg, out_cfg, min, max, err1, err2)
   
if __name__ == "__main__":
    main()
