#pragma once

#include <stdbool.h>

#include "types.h"

#ifdef FP64
    typedef f64 real;
    typedef __m128d f128;
    typedef __m256d f256;
    typedef __m512d f512;
#else
    typedef f32 real;
    typedef __m128 f128;
    typedef __m256 f256;
    typedef __m512 f512;
#endif

typedef struct {
    u64 nb_bodies;
    u32 nb_iter;
    u32 nb_warmups;
    real dt;
    char *output;
    bool debug;
    bool bench;
    bool check;
} config_t;

config_t config_new();
config_t config_from(int argc, char **argv);
void config_print(const config_t cfg);

void help(const char *executable);
