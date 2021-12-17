#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#include "types.h"

void init(f32 *x, f32 *y, f32 *z,
          f32 *vx, f32 *vy, f32 *vz,
          f32 *softening, f32 softening_value,
          f32 *dt, f32 dt_value, u64 n)
{
	for (u64 i = 0; i < n; i++) {
		u64 r1 = (u64)rand();
		u64 r2 = (u64)rand();
		f32 sign = (r1 > r2) ? 1 : -1;

		x[i] = sign * (f32)rand() / (f32)RAND_MAX;
		y[i] = (f32)rand() / (f32)RAND_MAX;
		z[i] = sign * (f32)rand() / (f32)RAND_MAX;

		vx[i] = (f32)rand() / (f32)RAND_MAX;
		vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
		vz[i] = (f32)rand() / (f32)RAND_MAX;
	}

	for (u64 i = 0; i < 16 ; i++) {
		softening[i] = softening_value;
		dt[i] = dt_value;
	}
}

void move_particles(f32 *x, f32 *y, f32 *z,
                    f32 *vx, f32 *vy, f32 *vz,
                    f32 *softening, f32 *dt, u64 n)
{
    __m512 pxi, pyi, pzi;
    __m512 pxj, pyj, pzj;
    __m512 vxi, vyi, vzi; 
    __m512 fx, fy, fz;
    __m512 dx, dy, dz;
    __m512 d_2, d_2xy, d_2zs, d_2_rsqrt, d_3_over_2;
    __m512 vsoft, vdt;

    pxi = _mm512_setzero_ps();
    pyi = _mm512_setzero_ps();
    pzi = _mm512_setzero_ps();

    pxj = _mm512_setzero_ps();
    pyj= _mm512_setzero_ps();
    pzj = _mm512_setzero_ps();

    vxi = _mm512_setzero_ps();
    vyi = _mm512_setzero_ps();
    vzi = _mm512_setzero_ps();

    dx = _mm512_setzero_ps();
    dy = _mm512_setzero_ps();
    dz = _mm512_setzero_ps();
    d_2 = _mm512_setzero_ps();
    d_2xy = _mm512_setzero_ps();
    d_2zs = _mm512_setzero_ps();
    d_2_rsqrt = _mm512_setzero_ps();
    d_3_over_2 = _mm512_setzero_ps();

    vsoft = _mm512_load_ps(softening);
    vdt = _mm512_load_ps(dt);

    for (u64 i = 0; i < n; i++) {
		fx = _mm512_setzero_ps();
		fy = _mm512_setzero_ps();
		fz = _mm512_setzero_ps();

		// Load all i-bound values
		pxi = _mm512_loadu_ps(x + i);
		pyi = _mm512_loadu_ps(y + i);
		pzi = _mm512_loadu_ps(z + i);
		vxi = _mm512_loadu_ps(vx + i);
		vyi = _mm512_loadu_ps(vy + i);
		vzi = _mm512_loadu_ps(vz + i);

		for (u64 j = 0; j < n; j += 16) {
			pxj = _mm512_load_ps(x + j);
			pyj= _mm512_load_ps(y + j);
			pzj = _mm512_load_ps(z + j);

			dx = _mm512_sub_ps(pxj, pxi);
			dy = _mm512_sub_ps(pyj, pyi);
			dz = _mm512_sub_ps(pzj, pzi);

			dx = _mm512_mul_ps(dx, dx);
			dy = _mm512_mul_ps(dy, dy);
			dz = _mm512_mul_ps(dz, dz);

            d_2xy = _mm512_add_ps(dx, dy);
            d_2zs = _mm512_add_ps(dz, vsoft);
            d_2 = _mm512_add_ps(d_2xy, d_2zs);
            d_2_rsqrt = _mm512_rsqrt14_ps(d_2);
            d_3_over_2 = _mm512_mul_ps(d_2_rsqrt, d_2_rsqrt);
            d_3_over_2 = _mm512_mul_ps(d_3_over_2, d_2_rsqrt);

            fx = _mm512_fmadd_ps(dx, d_3_over_2, fx);
            fy = _mm512_fmadd_ps(dy, d_3_over_2, fy);
            fz = _mm512_fmadd_ps(dz, d_3_over_2, fz);
		}

        vxi = _mm512_fmadd_ps(vdt, fx, vxi);
        vyi = _mm512_fmadd_ps(vdt, fy, vyi);
        vzi = _mm512_fmadd_ps(vdt, fz, vzi);

		_mm512_storeu_ps(vx + i, vxi);
		_mm512_storeu_ps(vy + i, vyi);
		_mm512_storeu_ps(vz + i, vzi);
	}

	// 3 f32ing-point operations
	for (u64 i = 0; i < n; i += 16) {
        // Reload v_i values
    	pxi = _mm512_loadu_ps(x + i);
    	pyi = _mm512_loadu_ps(y + i);
    	pzi = _mm512_loadu_ps(z + i);
    	vxi = _mm512_loadu_ps(vx + i);
    	vyi = _mm512_loadu_ps(vy + i);
    	vzi = _mm512_loadu_ps(vz + i);

        pxi = _mm512_fmadd_ps(vdt, vxi, pxi);
        pyi = _mm512_fmadd_ps(vdt, vyi, pyi);
        pzi = _mm512_fmadd_ps(vdt, vzi, pzi);

		_mm512_storeu_ps(x + i, pxi);
		_mm512_storeu_ps(y + i, pyi);
		_mm512_storeu_ps(z + i, pzi);
	}
}

int main(int argc, char **argv)
{
	const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
	const u64 steps = 10;
	const f32 dt_value = 0.01;
	const f32 softening_value = 1e-20;
	f64 rate = 0.0, drate = 0.0;

	// Steps to skip for warm up
	const u64 warmup = 3;

	f32 *x = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *y = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *z = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *vx = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *vy = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *vz = aligned_alloc(64, sizeof(f32) * (n + 64));
	f32 *softening = aligned_alloc(64, sizeof(f32) * 16);
	f32 *dt = aligned_alloc(64, sizeof(f32) * 16);

	init(x, y,z, vx, vy, vz, softening, softening_value, dt, dt_value, n);

	const u64 sp = 6 * sizeof(f64) * n;
	printf("\n\033[1mTotal memory size:\033[0m %llu B, %.2lf KiB, %.2lf MiB\n\n",
	       sp, (f64)sp / 1024.0f, (f64)sp / 1048576.0f);
	printf("\033[1m%5s %10s %10s %8s\033[0m\n",
	       "Step", "Time, s", "Interact/s", "GFLOP/s");
	fflush(stdout);

	for (u64 i = 0; i < steps; i++) {
        // Measure
        const f64 start = omp_get_wtime();
        move_particles(x, y, z, vx, vy, vz, softening, dt, n);
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
	drate = sqrtl(drate / (f64)(steps - warmup) - (rate * rate));

	printf("-----------------------------------------------------\n");
	printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	       "Average performance:", "", rate, drate);
	printf("-----------------------------------------------------\n");

	free(x);
	free(y);
	free(z);
	free(vx);
	free(vy);
	free(vz);
	free(softening);
	free(dt);

	return 0;
}
