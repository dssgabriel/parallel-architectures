#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "utils.h"

static const usize BLK_SIZE = 4096;

static inline
f32 hsum(const f512 x)
{
    return _mm512_reduce_add_ps(x);
}

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
    const f512 softening = _mm512_set1_ps(1e-20);
    const f512 vdt = _mm512_set1_ps(dt);

    #pragma omp parallel
    for (usize blk_start = 0; blk_start < nb_bodies; blk_start += BLK_SIZE) {
        const usize blk_end = blk_start + BLK_SIZE;

        // 30 floating-point operations
        #pragma omp for schedule(guided)
        for (usize i = 0; i < nb_bodies; i++) {
            f512 fx = _mm512_setzero_ps();
            f512 fy = _mm512_setzero_ps();
            f512 fz = _mm512_setzero_ps();

            // Load i-bound values
            f512 pxi = _mm512_set1_ps(p->px[i]);
            f512 pyi = _mm512_set1_ps(p->py[i]);
            f512 pzi = _mm512_set1_ps(p->pz[i]);

            // 19 floating-point operations
            #pragma omp simd reduction(+:fx,fy,fz) 
            for (usize j = blk_start; j < blk_end; j += 16) {
                const f512 dx = _mm512_sub_ps(_mm512_loadu_ps(p->px + j), pxi); // 1
                const f512 dy = _mm512_sub_ps(_mm512_loadu_ps(p->py + j), pyi); // 2
                const f512 dz = _mm512_sub_ps(_mm512_loadu_ps(p->pz + j), pzi); // 3

                const f512 d_2 =
                    _mm512_add_ps(
                        _mm512_add_ps(
                            _mm512_add_ps(
                                _mm512_mul_ps(dx, dx),
                                _mm512_mul_ps(dy, dy)
                            ),
                            _mm512_mul_ps(dz, dz)
                        ),
                        softening); // 9

            #ifdef __AVX512ER__
                const f512 invd = _mm512_rsqrt28_ps(d_2); // 11
            #else
                const f512 invd = _mm512_rsqrt14_ps(d_2); // 11
            #endif
                const f512 invd_3 = _mm512_mul_ps(_mm512_mul_ps(invd, invd), invd); // 13

                fx = _mm512_fmadd_ps(dx, invd_3, fx); // 15
                fy = _mm512_fmadd_ps(dy, invd_3, fy); // 17
                fz = _mm512_fmadd_ps(dz, invd_3, fz); // 19
            }

            f32 hfx = hsum(fx); // 8
            f32 hfy = hsum(fy); // 16
            f32 hfz = hsum(fz); // 24
            p->vx[i] += dt * hfx; // 26
            p->vy[i] += dt * hfy; // 28
            p->vz[i] += dt * hfz; // 30
        }

        // 6 floating-point operations
        #pragma omp for schedule(guided)
        for (usize i = 0; i < nb_bodies; i += 16) {
            // Reload v_i values
        	f512 pxi = _mm512_loadu_ps(p->px + i);
        	f512 pyi = _mm512_loadu_ps(p->py + i);
        	f512 pzi = _mm512_loadu_ps(p->pz + i);
        	f512 vxi = _mm512_loadu_ps(p->vx + i);
        	f512 vyi = _mm512_loadu_ps(p->vy + i);
        	f512 vzi = _mm512_loadu_ps(p->vz + i);

            pxi = _mm512_fmadd_ps(vdt, vxi, pxi); // 8
            pyi = _mm512_fmadd_ps(vdt, vyi, pyi); // 10
            pzi = _mm512_fmadd_ps(vdt, vzi, pzi); // 12

            _mm512_storeu_ps(p->px + i, pxi);
            _mm512_storeu_ps(p->py + i, pyi);
            _mm512_storeu_ps(p->pz + i, pzi);
        }
    }
}

void particles_print(particles_t *p, const config_t cfg)
{
    if (!strcmp(cfg.output, "none"))
        return;

    FILE *fp;
    if (strcmp(cfg.output, "stdout")) {
        fp = fopen(cfg.output, "wb");
        if (!fp)
            fp = stdout;
    } else {
        fp = stdout;
    }

    char *precision = sizeof(real) == 4 ? "fp32" : "fp64";
    fprintf(fp, "%llu\t%u\t%lf\t%s\n\n", cfg.nb_bodies, cfg.nb_iter, cfg.dt, precision);

    for (usize i = 0; i < cfg.nb_bodies; i++)
        fprintf(fp, "%.8e\t%.8e\t%.8e\n", p->px[i], p->py[i], p->pz[i]);

    if (strcmp(cfg.output, "stdout") && fp != stdout)
        fclose(fp);

    if (strcmp(cfg.output, "stdout") && fp == stdout)
        fprintf(stderr, "\033[1;33mwarning:\033[0m failed to open file %s\n", cfg.output);
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

    f64 rate = 0.0f, drate = 0.0f;
    particles_t *p = particles_new(cfg.nb_bodies);
    if (!p)
        return fprintf(stderr, "\033[1;31merror:\033[0m failed to allocate particles\n"), 1;
    particles_init(p, cfg.nb_bodies);

    const u64 mem_size = cfg.nb_bodies * 6 * sizeof(f32);
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
