#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "utils.h"

typedef struct particles_s {
    f32 *px, *py, *pz;
    f32 *vx, *vy, *vz;
} particles_t;

particles_t *particles_new(const u64 nb_bodies)
{
    particles_t *p = aligned_alloc(64, nb_bodies * sizeof(particles_t));
    if (!p)
        goto particles_failed_alloc;

    p->px = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    p->py = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    p->pz = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    p->vx = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    p->vy = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    p->vz = aligned_alloc(64, (nb_bodies + (nb_bodies % 128)) * sizeof(f32));
    if (!p->px || !p->py || !p->pz || !p->px || !p->py || !p->pz)
        goto particles_failed_alloc;

    return p;

particles_failed_alloc:
    return NULL;
}

void particles_init(particles_t *p, const u64 nb_bodies)
{
    srand(0);

    for (usize i = 0; i < nb_bodies; i++) {
        u64 r1 = (u64)rand();
        u64 r2 = (u64)rand();
        f32 sign = (r1 > r2) ? 1 : -1;

        p->px[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->py[i] = (f32)rand() / (f32)RAND_MAX;
        p->pz[i] = sign * (f32)rand() / (f32)RAND_MAX;

        p->vx[i] = (f32)rand() / (f32)RAND_MAX;
        p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
        p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

void particles_update(particles_t *p, const u64 nb_bodies, const f32 dt)
{
    const f256 softening = _mm256_set1_ps(1e-20);
    const f256 vdt = _mm256_set1_ps(dt);

    // 6 floating-point operations
    for (usize i = 0; i < nb_bodies; i++) {
        f256 fx = _mm256_setzero_ps();
        f256 fy = _mm256_setzero_ps();
        f256 fz = _mm256_setzero_ps();

        // Load all i-bound values
        f256 pxi = _mm256_set1_ps(p->px[i]);
        f256 pyi = _mm256_set1_ps(p->py[i]);
        f256 pzi = _mm256_set1_ps(p->pz[i]);
        f256 vxi = _mm256_set1_ps(p->vx[i]);
        f256 vyi = _mm256_set1_ps(p->vy[i]);
        f256 vzi = _mm256_set1_ps(p->vz[i]);

        // 19 floating-point operations
        for (usize j = 0; j < nb_bodies; j += 8) {
            // From here: p_i = d_, p_j = d_^2
            const f256 dx = _mm256_sub_ps(_mm256_loadu_ps(p->px + j), pxi); // 1
            const f256 dy = _mm256_sub_ps(_mm256_loadu_ps(p->py + j), pyi); // 2
            const f256 dz = _mm256_sub_ps(_mm256_loadu_ps(p->pz + j), pzi); // 3

            const f256 d_2 =
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(dx, dx),
                            _mm256_mul_ps(dy, dy)
                        ),
                        _mm256_mul_ps(dz, dz)
                    ),
                    softening); // 9

            const f256 invd = _mm256_rsqrt_ps(d_2); // 11
            const f256 invd_3 = _mm256_mul_ps(_mm256_mul_ps(invd, invd), invd); // 13

            fx = _mm256_fmadd_ps(dx, invd_3, fx); // 15
            fy = _mm256_fmadd_ps(dy, invd_3, fy); // 17
            fz = _mm256_fmadd_ps(dz, invd_3, fz); // 19
        }

        p->vx[i] = _mm256_cvtss_f32(_mm256_fmadd_ps(vdt, fx, vxi)); // 2
        p->vy[i] = _mm256_cvtss_f32(_mm256_fmadd_ps(vdt, fy, vyi)); // 4
        p->vz[i] = _mm256_cvtss_f32(_mm256_fmadd_ps(vdt, fz, vzi)); // 6
    }

    // 6 floating-point operations
    for (usize i = 0; i < nb_bodies; i += 8) {
        // Reload v_i values
    	f256 pxi = _mm256_loadu_ps(p->px + i);
    	f256 pyi = _mm256_loadu_ps(p->py + i);
    	f256 pzi = _mm256_loadu_ps(p->pz + i);
    	f256 vxi = _mm256_loadu_ps(p->vx + i);
    	f256 vyi = _mm256_loadu_ps(p->vy + i);
    	f256 vzi = _mm256_loadu_ps(p->vz + i);

        pxi = _mm256_fmadd_ps(vdt, vxi, pxi); // 8
        pyi = _mm256_fmadd_ps(vdt, vyi, pyi); // 10
        pzi = _mm256_fmadd_ps(vdt, vzi, pzi); // 12

        _mm256_storeu_ps(p->px + i, pxi);
        _mm256_storeu_ps(p->py + i, pyi);
        _mm256_storeu_ps(p->pz + i, pzi);
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

    for (usize i = 0; i < cfg.nb_bodies; i++)
        fprintf(fp, "%.8e\t%.8e\t%.8e\n", p->px[i], p->py[i], p->pz[i]);

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

    const u64 mem_size = cfg.nb_bodies * 6 * sizeof(f32);
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
