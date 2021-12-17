#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#include "types.h"

void init(float *x, float *y, float *z,
          float *vx, float *vy, float *vz,
          float *softening, float softening_value,
          float *dt, float dt_value, u64 n)
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

void move_particles(float *x, float *y, float *z,
                    float *vx, float *vy, float *vz,
                    float *softening, float *dt, u64 n)
{
    __m512 rxi, ryi, rzi, rxj, ryj, rzj;
    __m512 rvxi, rvyi, rvzi; 
    __m512 rfx, rfy, rfz;
    __m512 rdx, rdy, rdz;
    __m512 rd_2, rd_2xy, rd_2zs, rd_2_rsqrt, rsoft, rdt;

    rxi = _mm512_setzero_ps();
    ryi = _mm512_setzero_ps();
    rzi = _mm512_setzero_ps();

    rxj = _mm512_setzero_ps();
    ryj = _mm512_setzero_ps();
    rzj = _mm512_setzero_ps();

    rvxi = _mm512_setzero_ps();
    rvyi = _mm512_setzero_ps();
    rvzi = _mm512_setzero_ps();

    rdx = _mm512_setzero_ps();
    rdy = _mm512_setzero_ps();
    rdz = _mm512_setzero_ps();
    rd_2 = _mm512_setzero_ps();
    rd_2xy = _mm512_setzero_ps();
    rd_2zs = _mm512_setzero_ps();
    rd_2_rsqrt = _mm512_setzero_ps();

    rsoft = _mm512_load_ps(softening);
    rdt = _mm512_load_ps(dt);

    for (u64 i = 0; i < n; i++) {
		rfx = _mm512_setzero_ps();
		rfy = _mm512_setzero_ps();
		rfz = _mm512_setzero_ps();

		// Load i-bound values
		rxi = _mm512_loadu_ps(&x[i]);
		ryi = _mm512_loadu_ps(&y[i]);
		rzi = _mm512_loadu_ps(&z[i]);
		rvxi = _mm512_loadu_ps(&vx[i]);
		rvyi = _mm512_loadu_ps(&vy[i]);
		rvzi = _mm512_loadu_ps(&vz[i]);

		for (u64 j = 0; j < n; j += 16) {
			rxj = _mm512_load_ps(&x[j]);
			ryj = _mm512_load_ps(&y[j]);
			rzj = _mm512_load_ps(&z[j]);

			rdx = _mm512_sub_ps(rxj, rxi);
			rdy = _mm512_sub_ps(ryj, ryi);
			rdz = _mm512_sub_ps(rzj, rzi);

			rdx = _mm512_sub_ps(rdx, rdx);
			rdy = _mm512_sub_ps(rdy, rdy);
			rdz = _mm512_add_ps(rdz, rdz);

            rd_2xy = _mm512_add_ps(rdx, rdy);
            rd_2zs = _mm512_add_ps(rdz, rsoft);
            rd_2 = _mm512_add_ps(rd_2xy, rd_2zs);
            rd_2_rsqrt = _mm512_rsqrt14_ps(rd_2);

            rfx = _mm512_fmadd_ps(rdx, rd_2_rsqrt, rfx);
            rfy = _mm512_fmadd_ps(rdy, rd_2_rsqrt, rfy);
            rfz = _mm512_fmadd_ps(rdz, rd_2_rsqrt, rfz);
		}

        rvxi = _mm512_fmadd_ps(rdt, rfx, rvxi);
        rvyi = _mm512_fmadd_ps(rdt, rfy, rvyi);
        rvzi = _mm512_fmadd_ps(rdt, rfz, rvzi);

		_mm512_storeu_ps(&vx[i], rvxi);
		_mm512_storeu_ps(&vy[i], rvyi);
		_mm512_storeu_ps(&vz[i], rvzi);
	}

	// 3 floating-point operations
	for (u64 i = 0; i < n; i++) {
        // Reload v_i values
    	rxi = _mm512_loadu_ps(&x[i]);
    	ryi = _mm512_loadu_ps(&y[i]);
    	rzi = _mm512_loadu_ps(&z[i]);
    	rvxi = _mm512_loadu_ps(&vx[i]);
    	rvyi = _mm512_loadu_ps(&vy[i]);
    	rvzi = _mm512_loadu_ps(&vz[i]);

        rxi = _mm512_fmadd_ps(rdt, rvxi, rxi);
        ryi = _mm512_fmadd_ps(rdt, rvyi, ryi);
        rzi = _mm512_fmadd_ps(rdt, rvzi, rzi);

		_mm512_storeu_ps(&x[i], rxi);
		_mm512_storeu_ps(&y[i], ryi);
		_mm512_storeu_ps(&z[i], rzi);
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

	float *x = aligned_alloc(64, sizeof(float) * (n + 64));
	float *y = aligned_alloc(64, sizeof(float) * (n + 64));
	float *z = aligned_alloc(64, sizeof(float) * (n + 64));
	float *vx = aligned_alloc(64, sizeof(float) * (n + 64));
	float *vy = aligned_alloc(64, sizeof(float) * (n + 64));
	float *vz = aligned_alloc(64, sizeof(float) * (n + 64));
	float *softening = aligned_alloc(64, sizeof(float) * 16);
	float *dt = aligned_alloc(64, sizeof(float) * 16);

	init(x, y,z, vx, vy, vz, softening, softening_value, dt, dt_value, n);

	const u64 sx = 6 * sizeof(double) * n;
	printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", sx, sx >> 10, sx >> 20);
	printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s");
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
