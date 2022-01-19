#/bin/python

import functools
import statistics
import sys
from pathlib import Path

def helper():
    print("precision_checker — Python script that checks the precision of a 3D Nbody simulation\n")
    print("USAGE:\n    ./precision_check.py <FILE>\n")
    print("FILE:\n    Path to the output data file containing the particles' positions.")
    print("    NOTE: Simulation had to be ran with the same number of bodies and iterations as the",
          "reference (located at `data/ref.dat`)")


def check_files(ref, out):
    if not ref.is_file():
        print(f"error: reference data file `{ref}` does not exist")
        raise SystemExit

    if not out.is_file():
        print(f"error: data file `{out}` does not exist")
        raise SystemExit


def read_files(ref, out):
    ref_contents = ref.readlines()
    out_contents = out.readlines()

    if len(ref_contents) != len(out_contents):
        print(f"error: simulation do not have the same number of particles (ref is {len(ref_contents)}, out is {len(out_contents)})")
        raise SystemExit

    for i in range(0, len(ref_contents)):
        ref_contents[i] = ref_contents[i].split()
        out_contents[i] = out_contents[i].split()
    
    return ref_contents, out_contents


def process_files(ref, out):
    deltas = []
    for i in range(0, len(ref)):
        ref[i] = [float(x) for x in ref[i]]
        out[i] = [float(x) for x in out[i]]
        functools.reduce(lambda x, y: x+y, ref[i])
        functools.reduce(lambda x, y: x+y, out[i])
        ref[i] = ref[i][0]
        out[i] = out[i][0]
        deltas.append(ref[i] - out[i])

    return deltas


def compute_stats(deltas):
    mean = statistics.fmean(deltas)
    return mean


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

    ref_contents, out_contents = read_files(ref_file, out_file)
    deltas = process_files(ref_contents, out_contents)
    mean = compute_stats(deltas)

    print(f"Average precision loss from reference: {mean}")
   
if __name__ == "__main__":
    main()
