#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "types.h"

typedef struct particle_s {
    f32 *x, *y, *z;
    f32 *vx, *vy, *vz;
} particle_t;

particle_t *init(u64 n)
{
    particle_t *p = aligned_alloc(64, sizeof(particle_t));

    p->x = aligned_alloc(64, n * sizeof(f32));
    p->y = aligned_alloc(64, n * sizeof(f32));
    p->z = aligned_alloc(64, n * sizeof(f32));

    p->vx = aligned_alloc(64, n * sizeof(f32));
    p->vy = aligned_alloc(64, n * sizeof(f32));
    p->vz = aligned_alloc(64, n * sizeof(f32));

    for (u64 i = 0; i < n; i++) {
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        f32 sign = (r1 > r2) ? 1 : -1;

        p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->y[i] = (f32)rand() / (f32)RAND_MAX;
        p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

        p->vx[i] = (f32)rand() / (f32)RAND_MAX;
        p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }

    return p;
}

void move_particles(particle_t *p, const f32 dt, u64 n)
{
    const f32 softening = 1e-20;

    for (u64 i = 0; i < n; i++) {
        f32 __attribute__((aligned(64))) fx[4] = { 0.0 };
        f32 __attribute__((aligned(64))) fy[4] = { 0.0 };
        f32 __attribute__((aligned(64))) fz[4] = { 0.0 };

        // 23 floating-point operations
        for (u64 j = 0; j < n; j += 4) {
            // Newton's law
            f32 __attribute__((aligned(64))) dx[4];
            f32 __attribute__((aligned(64))) dy[4];
            f32 __attribute__((aligned(64))) dz[4];
            f32 __attribute__((aligned(64))) d_2[4];
            f32 __attribute__((aligned(64))) d_3_over_2[4];

            dx[0] = p->x[j] - p->x[i]; // 1
            dx[1] = p->x[j + 1] - p->x[i + 1]; // 1
            dx[2] = p->x[j + 2] - p->x[i + 2]; // 1
            dx[3] = p->x[j + 3] - p->x[i + 3]; // 1

            dy[0] = p->y[j] - p->y[i]; // 1
            dy[1] = p->y[j + 1] - p->y[i + 1]; // 1
            dy[2] = p->y[j + 2] - p->y[i + 2]; // 1
            dy[3] = p->y[j + 3] - p->y[i + 3]; // 1

            dz[0] = p->z[j] - p->z[i]; // 1
            dz[1] = p->z[j + 1] - p->z[i + 1]; // 1
            dz[2] = p->z[j + 2] - p->z[i + 2]; // 1
            dz[3] = p->z[j + 3] - p->z[i + 3]; // 1

            d_2[0] = (dx[0] * dx[0]) + (dy[0] * dy[0]) + (dz[0] * dz[0]) + softening; // 9
            d_2[1] = (dx[1] * dx[1]) + (dy[1] * dy[1]) + (dz[1] * dz[1]) + softening; // 9
            d_2[2] = (dx[2] * dx[2]) + (dy[2] * dy[2]) + (dz[2] * dz[2]) + softening; // 9
            d_2[3] = (dx[3] * dx[3]) + (dy[3] * dy[3]) + (dz[3] * dz[3]) + softening; // 9

            d_3_over_2[0] = sqrtf(d_2[0]) * sqrt(d_2[0]) * sqrt(d_2[0]); // 11
            d_3_over_2[1] = sqrtf(d_2[1]) * sqrt(d_2[1]) * sqrt(d_2[1]); // 11
            d_3_over_2[2] = sqrtf(d_2[2]) * sqrt(d_2[2]) * sqrt(d_2[2]); // 11
            d_3_over_2[3] = sqrtf(d_2[3]) * sqrt(d_2[3]) * sqrt(d_2[3]); // 11

            // Net force
            fx[0] += dx[0] * (1 / d_3_over_2[0]); // 13
            fx[1] += dx[1] * (1 / d_3_over_2[1]); // 13
            fx[2] += dx[2] * (1 / d_3_over_2[2]); // 13
            fx[3] += dx[3] * (1 / d_3_over_2[3]); // 13

            fy[0] += dy[0] * (1 / d_3_over_2[0]); // 15
            fy[1] += dy[1] * (1 / d_3_over_2[1]); // 15
            fy[2] += dy[2] * (1 / d_3_over_2[2]); // 15
            fy[3] += dy[3] * (1 / d_3_over_2[3]); // 15

            fz[0] += dz[0] * (1 / d_3_over_2[0]); // 17
            fz[1] += dz[1] * (1 / d_3_over_2[1]); // 17
            fz[2] += dz[2] * (1 / d_3_over_2[2]); // 17
            fz[3] += dz[3] * (1 / d_3_over_2[3]); // 17
    	}

        p->vx[i] += dt * (fx[0] + fx[1] + fx[2] + fx[3]); // 19
        p->vy[i] += dt * (fy[0] + fy[1] + fy[2] + fy[3]); // 21
        p->vz[i] += dt * (fz[0] + fz[1] + fz[2] + fz[3]); // 23
    }

    // 3 floating-point operations
    for (u64 i = 0; i < n; i++) {
        p->x[i] += dt * p->vx[i];
        p->y[i] += dt * p->vy[i];
        p->z[i] += dt * p->vz[i];
    }
}

int main(int argc, char **argv)
{
    const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
    const u64 steps = 10;
    const f32 dt = 0.01;

    f64 rate = 0.0, drate = 0.0;

    // Steps to skip for warm up
    const u64 warmup = 3;

    particle_t *p = init(n);
    const u64 s = sizeof(particle_t) * n;

    printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
    printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

    for (u64 i = 0; i < steps; i++) {
        // Measure
        const f64 start = omp_get_wtime();
        move_particles(p, dt, n);
        const f64 end = omp_get_wtime();

        // Number of interactions/iterations
        const f32 h1 = (f32)(n) * (f32)(n - 1);
        // GFLOPS
        const f32 h2 = (23.0 * h1 + 3.0 * (f32)n) * 1e-9;

        if (i >= warmup) {
            rate += h2 / (end - start);
            drate += (h2 * h2) / ((end - start) * (end - start));
	    }

        printf("%5llu %10.3e %10.3e %8.1f %s\n",
               i,
               (end - start),
               h1 / (end - start),
               h2 / (end - start),
               (i < warmup) ? "*" : "");
      
        fflush(stdout);
    }

    rate /= (f64)(steps - warmup);
    drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

    printf("-----------------------------------------------------\n");
    printf("\033[1m%s\033[30;42m      %10.1lf +- %.1lf GFLOP/s \033[0m\n",
           "Average performance:", rate, drate);
    printf("-----------------------------------------------------\n");

    free(p);
    return 0;
}
