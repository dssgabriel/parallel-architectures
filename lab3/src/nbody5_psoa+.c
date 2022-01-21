#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "utils.h"

static const u64 BLOCK = 4096;

typedef struct particles_s {
    real *px, *py, *pz;
    real *vx, *vy, *vz;
} particles_t;

particles_t *particles_new(const u64 nb_bodies)
{
    particles_t *p = aligned_alloc(64, nb_bodies * sizeof(particles_t));
    if (!p)
        goto particles_failed_alloc;

    p->px = aligned_alloc(64, nb_bodies * sizeof(real));
    p->py = aligned_alloc(64, nb_bodies * sizeof(real));
    p->pz = aligned_alloc(64, nb_bodies * sizeof(real));
    p->vx = aligned_alloc(64, nb_bodies * sizeof(real));
    p->vy = aligned_alloc(64, nb_bodies * sizeof(real));
    p->vz = aligned_alloc(64, nb_bodies * sizeof(real));
    if (!p->px || !p->py || !p->pz || !p->px || !p->py || !p->pz)
        goto particles_failed_alloc;

    return p;

particles_failed_alloc:
    return NULL;
}

void particles_init(particles_t *p, const u64 nb_bodies)
{
    srand(0);

    for (u64 i = 0; i < nb_bodies; i++) {
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        real sign = (r1 > r2) ? 1 : -1;

        p->px[i] = sign * (real)rand() / (real)RAND_MAX;
        p->py[i] = (real)rand() / (real)RAND_MAX;
        p->pz[i] = sign * (real)rand() / (real)RAND_MAX;

        p->vx[i] = (real)rand() / (real)RAND_MAX;
        p->vy[i] = sign * (real)rand() / (real)RAND_MAX;
        p->vz[i] = (real)rand() / (real)RAND_MAX;
    }
}

void particles_update(particles_t *p, const u64 nb_bodies, const real dt)
{
    const real softening = 1e-20;

    #pragma omp parallel
    for (usize blk_start = 0; blk_start < nb_bodies; blk_start += BLOCK) {
        usize blk_end = blk_start + BLOCK;

        // 6 floating-point operations
        #pragma omp for schedule(guided)
        for (usize i = 0; i < nb_bodies; i++) {
            real fx = 0.0f;
            real fy = 0.0f;
            real fz = 0.0f;

            const real pxi = p->px[i];
            const real pyi = p->py[i];
            const real pzi = p->pz[i];

            // 19 floating-point operations
            #pragma omp simd
            for (usize j = blk_start; j < blk_end; j++) {
                // Newton's law
                const real dx = p->px[j] - pxi; // 1
                const real dy = p->py[j] - pyi; // 2
                const real dz = p->pz[j] - pzi; // 3

                const real d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; // 9
                const real d_2_sqrt = sqrtf(d_2); // 10
                const real d_2_rsqrt = 1.0f / d_2_sqrt; // 11
                const real d_3_over_2 = d_2_rsqrt * d_2_rsqrt * d_2_rsqrt; // 13

                // Net force
                fx += dx * d_3_over_2; // 15
                fy += dy * d_3_over_2; // 17
                fz += dz * d_3_over_2; // 19
        	}

            p->vx[i] += dt * fx; // 2
            p->vy[i] += dt * fy; // 4
            p->vz[i] += dt * fz; // 6
        }

        // 6 floating-point operations
        #pragma omp for schedule(guided)
        for (usize i = 0; i < nb_bodies; i++) {
            p->px[i] += dt * p->vx[i]; // 8
            p->py[i] += dt * p->vy[i]; // 10
            p->pz[i] += dt * p->vz[i]; // 12
        }
    }
}

void particles_print(particles_t *p, const config_t cfg)
{
    FILE *fp;
    if (strcmp(cfg.output, "stdout")) {
        fp = fopen(cfg.output, "wb");
        if (!fp) {
            printf("warning: failed to open file %s\n", cfg.output);
            fp = stdout;
        }
    } else {
        fp = stdout;
    }

    if (sizeof(real) == 4) {
        for (usize i = 0; i < cfg.nb_bodies; i++)
            fprintf(fp, "%.8e\t%.8e\t%.8e\n", p->px[i], p->py[i], p->pz[i]);
    } else {
        for (usize i = 0; i < cfg.nb_bodies; i++)
            fprintf(fp, "%.16e\t%.16e\t%.16e\n", p->px[i], p->py[i], p->pz[i]);
    }

    if (strcmp(cfg.output, "stdout") && fp != stdout)
        fclose(fp);
}

void particles_drop(particles_t *p)
{
    free(p->px);
    free(p->py);
    free(p->pz);
    free(p->vx);
    free(p->vy);
    free(p->vz);
    free(p);
}

int main(int argc, char **argv)
{
    config_t cfg = (argc > 1) ? config_from(argc, argv) : config_new();
    if (cfg.debug)
        config_print(cfg);

    f64 rate = 0.0, drate = 0.0;
    particles_t *p = particles_new(cfg.nb_bodies);
    if (!p)
        return printf("error: failed to allocate particles\n"), 1;
    particles_init(p, cfg.nb_bodies);

    const u64 mem_size = cfg.nb_bodies * 6 * sizeof(real);
    fprintf(stderr,
            "\n\033[1mTotal memory size:\033[0m %llu B, %.2lf KiB, %.2lf MiB\n\n",
            mem_size, (f64)mem_size / 1024.0f, (f64)mem_size / 1048576.0f);
    fprintf(stderr,
            "\033[1m%5s %10s %10s %8s\033[0m\n",
            "Iter", "Time (s)", "Interact/s", "GFLOP/s");
    fflush(stderr);

    for (usize i = 0; i < cfg.nb_iter; i++) {
        // Measure
        const f64 start = omp_get_wtime();
        particles_update(p, cfg.nb_bodies, cfg.dt);
        const f64 end = omp_get_wtime();

        // Number of interactions/iterations
        const f64 h1 = (f64)(cfg.nb_bodies) * (f64)(cfg.nb_bodies - 1);
        // GFLOPS
        const f64 h2 = (19.0f * h1 + 12.0f * (f64)(cfg.nb_bodies)) * 1e-9;

        if (i >= cfg.nb_warmups) {
            rate += h2 / (end - start);
            drate += (h2 * h2) / ((end - start) * (end - start));
	    }

        fprintf(stderr,
                "%5zu %10.3e %10.3e %8.1f %s\n",
                i + 1,
                (end - start),
                h1 / (end - start),
                h2 / (end - start),
                (i < cfg.nb_warmups) ? "*" : "");
        fflush(stderr);
    }
    rate /= (f64)(cfg.nb_iter - cfg.nb_warmups);
    drate = sqrtl(drate / (f64)(cfg.nb_iter - cfg.nb_warmups) - (rate * rate));

    fprintf(stderr, "-----------------------------------------------------\n");
    fprintf(stderr, "\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
           "Average performance:", "", rate, drate);
    fprintf(stderr, "-----------------------------------------------------\n");

    particles_print(p, cfg);

    particles_drop(p);
    return 0;
}
