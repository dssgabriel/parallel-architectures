#! /usr/bin/python
from matplotlib import pyplot as plt
import numpy as np
import sys

bt = "\033[1m"
br = "\033[1;31m"
by = "\033[1;33m"
nt = "\033[0m"

def compute_mean_time(file_times):
    mean = 0.0
    dev = 0.0
    for time in file_times:
        mean += time[1]
        dev += time[1] * time[1]

    mean /= file_times[-1][0]
    dev = np.sqrt(dev / file_times[-1][0] - mean * mean)
    return mean, dev


def arrange_times(configs, times):
    gcc = []
    icc = []
    icx = []

    for i, t in enumerate(times):
        if configs[i][-1] == "gcc":
            gcc.append(t)
        elif configs[i][-1] == "icc":
            icc.append(t)
        elif configs[i][-1] == "icx":
            icx.append(t)

    return gcc, icc, icx


def arrange_flops(configs, flops):
    gcc = []
    icc = []
    icx = []

    for i, f in enumerate(flops):
        if configs[i][-1] == "gcc":
            gcc.append(f)
        elif configs[i][-1] == "icc":
            icc.append(f)
        elif configs[i][-1] == "icx":
            icx.append(f)

    return gcc, icc, icx


def main():
    if len(sys.argv) == 1:
        print(f"{br}error:{nt} wrong number of arguments\nNeeds at list one file\n")
        raise SystemExit

    files = [open(fp) for fp in sys.argv[1:-1]]
    configs = [fp.readline().split() for fp in files]
    for i, cfg in enumerate(configs):
        cfg.append(sys.argv[i + 1][7:10])

    flops = [fp.readline().split() for fp in files]
    for i, f in enumerate(flops):
        flops[i] = [float(x) for x in f]

    times = [fp.readlines() for fp in files]
    for i, f in enumerate(times):
        for j, t in enumerate(f):
            t = t.split()
            it = int(t[0]) + 1
            time = float(t[1])
            times[i][j] = [it, time]

    for i, file_times in enumerate(times):
        mean, dev = compute_mean_time(file_times)
        times[i] = [mean, dev]

    gcc, icc, icx = arrange_times(configs, times)
    gfs, ifs, xfs = arrange_flops(configs, flops)

    bwidth = 0.28
    fig = plt.subplots(figsize=(16, 9))
    plt.grid(linestyle='-', linewidth=0.5, alpha = 0.2)

    gccbar = np.arange(len(gcc))
    iccbar = [x + bwidth for x in gccbar]
    icxbar = [x + bwidth for x in iccbar]

    plt.bar(gccbar, [i[0] for i in gcc], color="indianred", width=bwidth, edgecolor="black", label="GCC", yerr=[i[1] for i in gcc])
    for i in range(len(gccbar)):
        plt.text(gccbar[i] - 0.12, gcc[i][0] + 0.33, f"{round(gfs[i][0], 1)} GFLOP/s", fontsize=8)

    plt.bar(iccbar, [i[0] for i in icc], color="mediumseagreen" , width=bwidth, edgecolor="black", label="ICC", yerr=[i[1] for i in icc])
    for i in range(len(iccbar)):
        plt.text(iccbar[i] - 0.12, icc[i][0] + 0.33, f"{round(ifs[i][0], 1)} GFLOP/s", fontsize=8)

    plt.bar(icxbar, [i[0] for i in icx], color="steelblue" , width=bwidth, edgecolor="black", label="ICX", yerr=[i[1] for i in icx])
    for i in range(len(icxbar)):
        plt.text(icxbar[i] - 0.12, icx[i][0] + 0.33, f"{round(xfs[i][0], 1)} GFLOP/s", fontsize=8)

    # Adding Xticks
    plt.xlabel('Number of particles', fontweight ='bold', fontsize = 15)
    plt.ylabel('Latency in seconds', fontweight ='bold', fontsize = 15)
    plt.xticks([r + bwidth for r in range(len(gcc))], [c[0] for c in configs])
    plt.legend()
    plt.savefig(sys.argv[-1])
    plt.show()
  
if __name__ == "__main__":
    main()
