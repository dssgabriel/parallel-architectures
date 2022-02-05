import argparse
import numpy as np
import statistics
import sys
from pathlib import Path

# Define escape characters for bold, bold red, bold yellow and normal text.
BTXT = "\033[1m"
BRED = "\033[1;31m"
BYLW = "\033[1;33m"
NTXT = "\033[0m"

def check_configs(ref, cmp):
    if ref[0] != cmp[0]:
        print(f"{BRED}error:{NTXT} simulations do not have the same number",
            f"of particles (reference is {ref[0]}, comparative is {cmp[0]})")
        return False
    if ref[1] != cmp[1]:
        print(f"{BYLW}warning:{NTXT} simulations do not have the same number",
            f"of iterations (reference is {ref[1]}, comparative is {cmp[1]})")
    if ref[2] != cmp[2]:
        print(f"{BYLW}warning:{NTXT} simulations do not use the same time step",
            f"(reference is {ref[2]}, comparative is {cmp[2]})")
    if ref[3] != cmp[3]:
        print(f"{BYLW}warning:{NTXT} simulations do not use the same precision",
            f"(reference is {ref[3]}, comparative is {cmp[3]})")
    return True


def compute_stats(ref, cmp, nb_part):
    min = np.inf
    max = 0.0
    err1 = 0.0
    err2 = 0.0
    for idx, (r, c) in enumerate(zip(ref, cmp)):
        x1 = abs(r[0] - c[0])
        y1 = abs(r[1] - c[1])
        z1 = abs(r[2] - c[2])
        x2 = (r[0] - c[0]) ** 2
        y2 = (r[1] - c[1]) ** 2
        z2 = (r[2] - c[2]) ** 2
        sum1 = x1 + y1 + z1
        sum2 = x2 + y2 + z2
        if sum1 > max:
            max = sum1
        if sum1 < min:
            min = sum1
        err1 += sum1
        err2 += np.sqrt(sum2)
        
    err1 /= nb_part
    err2 /= nb_part
    return min, max, err1, err2


def print_results(ref, cmp, min, max, err1, err2):
    print(f"Comparison of {ref} and {cmp}:")
    print(f"Min error: {min}")
    print(f"Max error: {max}")
    print(f"Relative error (norm 1): {err1}")
    print(f"Relative error (norm 2): {err2}") 
    

def main():
    parser = argparse.ArgumentParser(description="A script that checks the\
    precision between two N-body simulations")
    parser.add_argument("ref", metavar="REFERENCE", nargs=1,
        help="path to reference simulation output")
    parser.add_argument("cmp", metavar="COMPARATIVE", nargs=1,
        help="path to comparative simulation output")

    args = parser.parse_args()
    ref = Path(args.ref[0])
    cmp = Path(args.cmp[0])
    if not ref.is_file() or not cmp.is_file():
        raise ValueError("No such file or directory")

    with open(ref, "r") as f:
        r_cfg = f.readline().split()
        r_cts = []
        for line in f.readlines():
            r_cts.append([float(x) for x in line.split()])
    with open(cmp, "r") as f:
        c_cfg = f.readline().split()
        c_cts = []
        for line in f.readlines():
            c_cts.append([float(x) for x in line.split()])
    if not check_configs(r_cfg, c_cfg):
        raise ValueError("Incompatible configurations")

    min, max, err1, err2 = compute_stats(r_cts, c_cts, int(r_cfg[0]))
    print_results(ref, cmp, min, max, err1, err2)
   
if __name__ == "__main__":
    main()
