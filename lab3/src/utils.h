#pragma once

#include <stdbool.h>

#include "types.h"

#ifdef FP64
    typedef f64 real;
#else
    typedef f32 real;
#endif

typedef struct {
    u64 nb_bodies;
    u32 nb_iter;
    u32 nb_warmups;
    real dt;
    char *output;
    bool debug;
} config_t;

config_t config_new();
config_t config_from(int argc, char **argv);
void config_print(const config_t cfg);

void help(const char *executable);
