#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "utils.h"

static const usize offset = 7;

typedef struct particles_s {
    real *px, *py, *pz;
    real *vx, *vy, *vz;
} particles_t;

particles_t *particles_new(const u64 nb_bodies)
{
    particles_t *p = aligned_alloc(64, nb_bodies * sizeof(particles_t));
    if (!p)
        goto particles_failed_alloc;

    p->px = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    p->py = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    p->pz = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    p->vx = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    p->vy = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    p->vz = aligned_alloc(64, (nb_bodies + (offset * 2 + 1)) * sizeof(real));
    if (!p->px || !p->py || !p->pz || !p->px || !p->py || !p->pz)
        goto particles_failed_alloc;

    memset(p->px, 0, (nb_bodies + 15) * sizeof(real)); 
    memset(p->py, 0, (nb_bodies + 15) * sizeof(real)); 
    memset(p->pz, 0, (nb_bodies + 15) * sizeof(real)); 
    memset(p->vx, 0, (nb_bodies + 15) * sizeof(real)); 
    memset(p->vy, 0, (nb_bodies + 15) * sizeof(real)); 
    memset(p->vz, 0, (nb_bodies + 15) * sizeof(real)); 

    return p;

particles_failed_alloc:
    return NULL;
}

void particles_init(particles_t *p, const u64 nb_bodies)
{
    srand(0);

    for (u64 i = offset; i < nb_bodies + offset; i++) {
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

#ifdef FP64
void particles_update(particles_t *p, const u64 nb_bodies, const real dt)
{
    const f256 softening = _mm256_set1_pd(1e-20);
    const f256 vdt = _mm256_set1_pd(dt);
    const f256 one = _mm256_set1_pd(1.0f);

    // 25 floating-point operations
    for (usize i = 0; i < nb_bodies + (offset * 2); i++) {
        f256 fx = _mm256_setzero_pd();
        f256 fy = _mm256_setzero_pd();
        f256 fz = _mm256_setzero_pd();

        // Load all i-bound values
        f256 pxi = _mm256_loadu_pd(p->px + i);
        f256 pyi = _mm256_loadu_pd(p->py + i);
        f256 pzi = _mm256_loadu_pd(p->pz + i);
        f256 vxi = _mm256_loadu_pd(p->vx + i);
        f256 vyi = _mm256_loadu_pd(p->vy + i);
        f256 vzi = _mm256_loadu_pd(p->vz + i);

        for (usize j = offset; j < nb_bodies + offset; j += 4) {
            f256 pxj = _mm256_loadu_pd(p->px + j);
            f256 pyj = _mm256_loadu_pd(p->py + j);
            f256 pzj = _mm256_loadu_pd(p->pz + j);

            // From here: p_i = d_, p_j = d_^2
            pxi = _mm256_sub_pd(pxj, pxi); // 1
            pyi = _mm256_sub_pd(pyj, pyi); // 2
            pzi = _mm256_sub_pd(pzj, pzi); // 3

            pxj = _mm256_mul_pd(pxi, pxi); // 4
            pyj = _mm256_mul_pd(pyi, pyi); // 5
            pzj = _mm256_mul_pd(pzi, pzi); // 6

            f256 d_2 = _mm256_add_pd(pxj, pyj); // 7
            d_2 = _mm256_add_pd(d_2, pzj); // 8
            d_2 = _mm256_add_pd(d_2, softening); // 9
            d_2 = _mm256_sqrt_pd(d_2); // 10
            d_2 = _mm256_div_pd(one, d_2); // 11
            const f256 tmp = _mm256_mul_pd(d_2, d_2); // 12
            d_2 = _mm256_mul_pd(tmp, d_2); // 13

            fx = _mm256_fmadd_pd(pxi, d_2, fx); // 15
            fy = _mm256_fmadd_pd(pyi, d_2, fy); // 17
            fz = _mm256_fmadd_pd(pzi, d_2, fz); // 19
        }

        vxi = _mm256_fmadd_pd(vdt, fx, vxi); // 21
        vyi = _mm256_fmadd_pd(vdt, fy, vyi); // 23
        vzi = _mm256_fmadd_pd(vdt, fz, vzi); // 25

        _mm256_storeu_pd(p->vx + i, vxi);
        _mm256_storeu_pd(p->vy + i, vyi);
        _mm256_storeu_pd(p->vz + i, vzi);
    }

    // 6 floating-point operations
    for (usize i = offset; i < nb_bodies + offset; i += 4) {
        // Reload v_i values
    	f256 pxi = _mm256_loadu_pd(p->px + i);
    	f256 pyi = _mm256_loadu_pd(p->py + i);
    	f256 pzi = _mm256_loadu_pd(p->pz + i);
    	f256 vxi = _mm256_loadu_pd(p->vx + i);
    	f256 vyi = _mm256_loadu_pd(p->vy + i);
    	f256 vzi = _mm256_loadu_pd(p->vz + i);

        pxi = _mm256_fmadd_pd(vdt, vxi, pxi); // 2
        pyi = _mm256_fmadd_pd(vdt, vyi, pyi); // 4
        pzi = _mm256_fmadd_pd(vdt, vzi, pzi); // 6

        _mm256_storeu_pd(p->px + i, pxi);
        _mm256_storeu_pd(p->py + i, pyi);
        _mm256_storeu_pd(p->pz + i, pzi);
    }
}
#else
void particles_update(particles_t *p, const u64 nb_bodies, const real dt)
{
    const f256 softening = _mm256_set1_ps(1e-20);
    const f256 vdt = _mm256_set1_ps(dt);

    // 25 floating-point operations
    for (usize i = 0; i < nb_bodies + (offset * 2); i++) {
        f256 fx = _mm256_setzero_ps();
        f256 fy = _mm256_setzero_ps();
        f256 fz = _mm256_setzero_ps();

        // Load all i-bound values
        f256 pxi = _mm256_loadu_ps(p->px + i);
        f256 pyi = _mm256_loadu_ps(p->py + i);
        f256 pzi = _mm256_loadu_ps(p->pz + i);
        f256 vxi = _mm256_loadu_ps(p->vx + i);
        f256 vyi = _mm256_loadu_ps(p->vy + i);
        f256 vzi = _mm256_loadu_ps(p->vz + i);

        for (usize j = 0; j < nb_bodies + (offset * 2); j += 8) {
            f256 pxj = _mm256_loadu_ps(p->px + j);
            f256 pyj = _mm256_loadu_ps(p->py + j);
            f256 pzj = _mm256_loadu_ps(p->pz + j);

            // From here: p_i = d_, p_j = d_^2
            pxi = _mm256_sub_ps(pxj, pxi); // 1
            pyi = _mm256_sub_ps(pyj, pyi); // 2
            pzi = _mm256_sub_ps(pzj, pzi); // 3

            pxj = _mm256_mul_ps(pxi, pxi); // 4
            pyj = _mm256_mul_ps(pyi, pyi); // 5
            pzj = _mm256_mul_ps(pzi, pzi); // 6

            f256 d_2 = _mm256_add_ps(pxj, pyj); // 7
            d_2 = _mm256_add_ps(d_2, pzj); // 8
            d_2 = _mm256_add_ps(d_2, softening); // 9
            d_2 = _mm256_rsqrt_ps(d_2); // 11
            const f256 tmp = _mm256_mul_ps(d_2, d_2); // 12
            d_2 = _mm256_mul_ps(tmp, d_2); // 13

            fx = _mm256_fmadd_ps(pxi, d_2, fx); // 15
            fy = _mm256_fmadd_ps(pyi, d_2, fy); // 17
            fz = _mm256_fmadd_ps(pzi, d_2, fz); // 19
        }

        vxi = _mm256_fmadd_ps(vdt, fx, vxi); // 21
        vyi = _mm256_fmadd_ps(vdt, fy, vyi); // 23
        vzi = _mm256_fmadd_ps(vdt, fz, vzi); // 25

        _mm256_storeu_ps(p->vx + i, vxi);
        _mm256_storeu_ps(p->vy + i, vyi);
        _mm256_storeu_ps(p->vz + i, vzi);
    }

    // 6 floating-point operations
    for (usize i = offset; i < nb_bodies + offset; i += 8) {
        // Reload v_i values
    	f256 pxi = _mm256_loadu_ps(p->px + i);
    	f256 pyi = _mm256_loadu_ps(p->py + i);
    	f256 pzi = _mm256_loadu_ps(p->pz + i);
    	f256 vxi = _mm256_loadu_ps(p->vx + i);
    	f256 vyi = _mm256_loadu_ps(p->vy + i);
    	f256 vzi = _mm256_loadu_ps(p->vz + i);

        pxi = _mm256_fmadd_ps(vdt, vxi, pxi); // 2
        pyi = _mm256_fmadd_ps(vdt, vyi, pyi); // 4
        pzi = _mm256_fmadd_ps(vdt, vzi, pzi); // 6

        _mm256_storeu_ps(p->px + i, pxi);
        _mm256_storeu_ps(p->py + i, pyi);
        _mm256_storeu_ps(p->pz + i, pzi);
    }
}
#endif

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
        const f64 h2 = (25.0f * h1 + 6.0f * (f64)(cfg.nb_bodies)) * 1e-9;

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
