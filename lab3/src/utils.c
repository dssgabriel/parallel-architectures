#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utils.h"

config_t config_new()
{
    config_t cfg = {
        .nb_bodies = 16384,
        .nb_iter = 10,
        .nb_warmups = 3,
        .dt = 0.01f,
        .output = "stdout",
        .debug = false,
    };

    return cfg;
}

config_t config_from(int argc, char **argv)
{
    config_t cfg = config_new();

    int curr;
    while ((curr = getopt(argc, argv, "n:i:w:t:o:dh")) != -1) {
        switch (curr) {
            case 'n':
                cfg.nb_bodies = atoll(optarg);
                break;
            case 'i':
                cfg.nb_iter = atol(optarg);
                break;
            case 'w':
                cfg.nb_warmups = atol(optarg);
                break;
            case 't':
                cfg.dt = atof(optarg);
                break;
            case 'o':
                cfg.output = optarg;
                break;
            case 'd':
                cfg.debug = true;
                break;
            case 'h':
                help(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                printf("See help below for more information\n\n");
                help(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    return cfg;
}

void config_print(const config_t cfg)
{
    char *precision;
    if (sizeof(real) == 8)
        precision = "f64";
    else
        precision = "f32";

    printf("CONFIG:\n\t%s:\t%llu\n\t%s:\t%u\n\t%s:\t%u\n\t%s:\t%lf\n\t%s:\t\t%s\n\t%s:\t\t%s\n",
           "Number of bodies", cfg.nb_bodies,
           "Number of iterations", cfg.nb_iter,
           "Number of warmups", cfg.nb_warmups,
           "Size of time step", cfg.dt,
           "Output file", cfg.output,
           "Precision used", precision);
}

void help(const char *exec)
{
    printf("%s — 3D Nbody simulation\n\n", exec);
    printf("USAGE:\n\t./%s [OPTIONS] [FLAGS]\n\n", exec);
    printf("FLAGS:\n\t-h\t\tPrints help information\n\t-d\t\tEnables debug mode\n\n");
    printf("OPTIONS:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n",
           "-n <nb_bodies>\tNumber of bodies in the simulation",
           "-i <nb_iter>\tNumber of iterations performed in total",
           "-w <nb_warmups>\tNumber of warmup iterations",
           "-t <dt>\t\tSize of time step (discretization)",
           "-o <output>\tOutput file");
}
