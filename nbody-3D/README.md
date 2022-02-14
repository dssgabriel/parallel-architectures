# Optimizing a 3D N-body simulation

## Introduction
The aim of this lab is to optimize the base version of the algorithm (`src/nbody0.c`) and benchmark the improvements obtained by each iteration of optimizations.

The _N-body_ algorithm is quite simple, it simple consists in computing every particle-particle interactions for each step of the simulation using Newton's laws of classical mechanics. Given a _N_ objects system, this results in a _O(n²)_ arithmetic complexity.

The base implementation uses an _Array of Structures_ (AoS) to represent a particle:
```c
struct particle_t {
    float pos_x;
    float pos_y;
    float pos_z;
    float vel_x;
    float vel_y;
    float vel_z;
};
```

## Contribution
### Converting into _Structure of Arrays_ (SoA)
The first optimization that was done was to convert the basic AoS reprensentation into SoA, which is much more efficient in terms of memory usage.
```c
struct particles_t {
    float *pos_x;
    float *pos_y;
    float *pos_z;
    float *vel_x;
    float *vel_y;
    float *vel_z;
};
```
See `src/nbody1_soa.c` for the SoA version.

### Optimizing the hot function: `move_particles`
When looking at the inner loop of the function, we can notice that we are re-computing the address of `i`-bound values at each iteration, which is seriously inefficient:
```c
for (size_t i = 0; i < N; i++) {
    // -- snip --
    
    for (size_t j = 0; j < N; j++) {
        float dx = p->pos_x[j] - p->pos_x[i]; // {
        float dy = p->pos_y[j] - p->pos_y[i]; // { inefficient
        float dz = p->pos_z[j] - p->pos_z[i]; // {
        
        // -- snip --
    }
}
```
Removing these and computing them in the outer loop is a first step towards optimization:
```c
for (size_t i = 0; i < N; i++) {
    // -- snip --
    
    float pos_x_i = p->pos_x[i];
    float pos_y_i = p->pos_y[i];
    float pos_z_i = p->pos_z[i];
    
    for (size_t j = 0; j < N; j++) {
        float dx = p->pos_x[j] - pos_x_i;
        float dy = p->pos_y[j] - pos_y_i;
        float dz = p->pos_z[j] - pos_z_i;
        
        // -- snip --
    }
}
```

The inner loop can be optimized further by removing the call to the `pow` function from C's `math.h` header.
```c
float d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening;
float d_3_over_2 = pow(d_2, 3.0f / 2.0f); // very slow

float fx += dx / d_3_over_2;
float fy += dy / d_3_over_2;
float fz += dz / d_3_over_2;
```
Replacing it by a call to `sqrtf` and an inversion allows for two major optimizations:
1. it makes the computation of `fx`, `fy` and `fz` a multiply-add operation, which can take advantage of the AVX `vfmadd` instructions;
2. it allows the compiler to compute the reciprocal square root using the hardware-optimized `vrsqrt` instructions, which *significantly* improve the performance.
```c
float d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening;
d_2 = sqrtf(d_2);
d_2 = 1.0f / d_2;
float d_3_over_2 = d_2 * d_2 * d_2;

float fx += dx * d_3_over_2;
float fy += dy * d_3_over_2;
float fz += dz * d_3_over_2;
```
See `src/nbody2_ssoai.c` for the sequential SoA improved version.

### Using AVX intrinsics to guarantee loop vectorization
By translating the entire `move_particles` function into calls to Intel provided _intrinsics_ function (`immintrin.h` header), we can guarantee that the compiler generates fully vectorized code that takes advantage of the entire AVX/AVX512 registers length (depending on the intrinsics used).

Both AVX and AVX512 have been implemented (sequentially), see `src/nbody3_savx.c` and `src/nbody4_savx512.c` respectively.

## Bonus
### Parallelization
The improved SoA and both intrinsics versions have been updated to be parallelized using `OpenMP`. The `src/nbody5_psoai.c`, `src/nbody6_pavx.c` and `src/nbody7_pavx512.c` are simply parallelized by adding a _pragma_ directive over the outer loop:
```c
#pragma omp parallel for
for (size_t i = 0; i < N; i++) {
    // -- snip --
    
    for (size_t j = 0; j < N; j++) {
        // -- snip --
    }
}
```

We also implemented parallel-block versions specifically designed to get maximum performance on Intel Xeon Phis processors (Knights Landing/Knights Mill) by exploiting their core architecture and memory layout. See `src/nbody8_pbsoai.c`, `src/nbody9_pbavx.c` and `src/nbodyA_pbavx512.c`.

Finally, an hybrid implementation using both `MPI` and `OpenMP` was designed to try and leverage of the performance on multiple cluster nodes. See `src/nbodyB_mpi.c`.
