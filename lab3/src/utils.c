#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "utils.h"

config_t config_new()
{
    config_t cfg = {
        .nb_bodies = 16384,
        .nb_iter = 10,
        .nb_warmups = 3,
        .dt = 0.01f,
        .output = "none",
        .debug = false,
        .bench = false,
        .check = false,
    };

    return cfg;
}

config_t config_from(int argc, char **argv)
{
    config_t cfg = config_new();

    int curr;
    while ((curr = getopt(argc, argv, "n:i:w:t:o:dbch")) != -1) {
        switch (curr) {
            case 'n':
                cfg.nb_bodies = atoll(optarg);
                break;
            case 'i':
                if (cfg.check == true) {
                    fprintf(stderr, "\033[1;33warning:\033[0m flag `-c` overrides `-i` option\n");
                    cfg.nb_iter = 1;
                    cfg.nb_warmups = 0;
                } else {
                    cfg.nb_iter = atol(optarg);
                }
                break;
            case 'w':
                if (cfg.check == true) {
                    fprintf(stderr, "\033[1;33warning:\033[0m flag `-c` overrides `-i` option\n");
                    cfg.nb_iter = 1;
                    cfg.nb_warmups = 0;
                } else {
                    cfg.nb_warmups = atol(optarg);
                }
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
            case 'b':
                if (cfg.check == true) {
                    fprintf(stderr, "\033[1;31merror:\033[0m \033[1mbenchmark\033[0m mode is not compatible with \033[1mprecision check\033[0m mode\n");
                    exit(EXIT_FAILURE);
                }
                cfg.bench = true;
                break;
            case 'c':
                if (cfg.bench == true) {
                    fprintf(stderr, "\033[1;31merror:\033[0m \033[1mprecision check\033[0m mode is not compatible with \033[1mbenchmark\033[0m mode\n");
                    exit(EXIT_FAILURE);
                }
                cfg.check = true;
                cfg.nb_iter = 1;
                cfg.nb_warmups = 0;
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

    if ((!strcmp(cfg.output, "none")) && (cfg.check == true || cfg.bench == true)) {
        fprintf(stderr, "\033[1;31merror:\033[0m flag `-c` requires the `-o` option\n");
        exit(EXIT_FAILURE);
    }

    return cfg;
}

void config_print(const config_t cfg)
{
    char *precision;
    if (sizeof(real) == 8)
        precision = "fp64";
    else
        precision = "fp32";

    fprintf(stderr,
            "\033[1mCONFIGURATION:\033[0m\n\t%s:\t%llu\n\t%s:\t%u\n\t%s:\t%u\n\t%s:\t%lf\n\t%s:\t\t%s\n\t%s:\t\t%s\n",
            "Number of bodies", cfg.nb_bodies,
            "Number of iterations", cfg.nb_iter,
            "Number of warmups", cfg.nb_warmups,
            "Size of time step", cfg.dt,
            "Output file", cfg.output,
            "Precision used", precision);
}

void help(const char *exec)
{
    fprintf(stderr, "%s — 3D N-body simulation\n\n", exec);
    fprintf(stderr, "\033[1mUSAGE:\033[0m\n\t./%s [FLAGS] [OPTIONS]\n\n", exec);
    fprintf(stderr,
            "\033[1mFLAGS:\033[0m\n\t%s\n\t%s\n\t%s\n\t%s\n\n",
            "-h\t\tPrints help information.",
            "-d\t\tEnables debug mode.",
            "-b\t\tEnables benchmark mode. Not compatible with `-c` flag. Requires `-o` option.",
            "-c\t\tEnables precision check mode. Not compatible with `-b` flag. Requires `-o` option.");
    fprintf(stderr,
            "\033[1mOPTIONS:\033[0m\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n",
            "-n <nb_bodies>\tNumber of particles in the simulation.",
            "-i <nb_iter>\tNumber of iterations performed in total.",
            "-w <nb_warmups>\tNumber of warmup iterations.",
            "-t <dt>\t\tSize of time step (discretization).",
            "-o <output>\tPath to the output file.");
}
