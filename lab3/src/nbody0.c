#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "utils.h"

typedef struct particle_s {
    real px, py, pz;
    real vx, vy, vz;
} particle_t;

particle_t *particles_new(const u64 nb_bodies)
{
    return malloc(nb_bodies * sizeof(particle_t));
}

void particles_init(particle_t *p, const u64 nb_bodies)
{
    srand(0);

    for (usize i = 0; i < nb_bodies; i++) {
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        real sign = (r1 > r2) ? 1.0f : -1.0f;

        p[i].px = sign * (real)(rand()) / (real)(RAND_MAX);
        p[i].py = (real)(rand()) / (real)(RAND_MAX);
        p[i].pz = sign * (real)(rand()) / (real)(RAND_MAX);

        p[i].vx = (real)(rand()) / (real)(RAND_MAX);
        p[i].vy = sign * (real)(rand()) / (real)(RAND_MAX);
        p[i].vz = (real)(rand()) / (real)(RAND_MAX);
    }
}

void particles_update(particle_t *p, const u64 nb_bodies, const real dt)
{
    const real softening = 1e-20;

    // 6 floating-point operations
    for (usize i = 0; i < nb_bodies; i++) {
        real fx = 0.0f;
        real fy = 0.0f;
        real fz = 0.0f;

        // 17 floating-point operations
        for (usize j = 0; j < nb_bodies; j++) {
            // Newton's law
            const real dx = p[j].px - p[i].px; // 1
            const real dy = p[j].py - p[i].py; // 2
            const real dz = p[j].pz - p[i].pz; // 3
            const real d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; // 9
            const real d_3_over_2 = pow(d_2, 3.0f / 2.0f); // 11

            // Net force
            fx += dx / d_3_over_2; // 13
            fy += dy / d_3_over_2; // 15
            fz += dz / d_3_over_2; // 17
    	}

        p[i].vx += dt * fx; // 2
        p[i].vy += dt * fy; // 4
        p[i].vz += dt * fz; // 6
    }

    // 6 floating-point operations
    for (usize i = 0; i < nb_bodies; i++) {
        p[i].px += dt * p[i].vx; // 8
        p[i].py += dt * p[i].vy; // 10
        p[i].pz += dt * p[i].vz; // 12
    }
}

void particles_print(particle_t *p, const config_t cfg)
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
            fprintf(fp, "%.8e\t%.8e\t%.8e\n", p[i].px, p[i].py, p[i].pz);
    } else {
        for (usize i = 0; i < cfg.nb_bodies; i++)
            fprintf(fp, "%.16e\t%.16e\t%.16e\n", p[i].px, p[i].py, p[i].pz);
    }

    if (strcmp(cfg.output, "stdout") && fp != stdout)
        fclose(fp);
}

int main(int argc, char **argv)
{
    config_t cfg = (argc > 1) ? config_from(argc, argv) : config_new();
    if (cfg.debug)
        config_print(cfg);

    f64 rate = 0.0f, drate = 0.0f;
    particle_t *p = particles_new(cfg.nb_bodies);
    if (!p)
        return printf("error: failed to allocate particles\n"), 1;
    particles_init(p, cfg.nb_bodies);

    const u64 mem_size = cfg.nb_bodies * sizeof(particle_t);
    fprintf(stderr,
            "\n\033[1mTotal memory size:\033[0m %llu B, %.2lf kB, %.2lf MB\n\n",
            mem_size, (f64)mem_size / 1000.0f, (f64)mem_size / 1000000.0f);
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
        const f64 h2 = (17.0f * h1 + 12.0f * (f64)(cfg.nb_bodies)) * 1e-9;

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

    if (strcmp(cfg.output, "none"))
        particles_print(p, cfg);

    free(p);
    return 0;
}
