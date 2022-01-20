import numpy as np
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
    for i in range(0, len(ref)):
        ref[i] = [float(x) for x in ref[i]]
        out[i] = [float(x) for x in out[i]]


def compute_stats(ref, out):
    norm1 = 0.0
    norm2 = 0.0

    max = 0.0
    min = np.inf

    for i in range(0, len(ref)):
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

        norm1 += sum1
        norm2 += np.sqrt(sum2)

    print(f"Maximum variation is: {max}")
    print(f"Minimum variation is: {min}")
    print(f"Error in 1-norm: {norm1 / len(ref)}")
    print(f"Error in 2-norm: {norm2 / len(ref)}")
    
    # maxd = max(deltas)
    # mind = min(deltas)
    # mean = statistics.mean(deltas)
    # print(f"Maximum precision loss: {maxd}")
    # print(f"Minimum precision loss: {mind}")
    # print(f"Average precision loss: {abs(mean)}")


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
    process_files(ref_contents, out_contents)
    compute_stats(ref_contents, out_contents)
   
if __name__ == "__main__":
    main()
