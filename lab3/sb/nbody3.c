//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
//typedef struct particle_s {
//
//  f32 x, y, z;
//  f32 vx, vy, vz;
//  
//} particle_t;

//
void init(float *x, float *y, float *z, float *vx, float *vy, float *vz, float *softening, float softening_value, float *dt, float dt_value,  u64 n)
{
	for (u64 i = 0; i < n; i++)
	{
		//
		u64 r1 = (u64)rand();
		u64 r2 = (u64)rand();
		f32 sign = (r1 > r2) ? 1 : -1;

		//
		x[i] = sign * (f32)rand() / (f32)RAND_MAX;
		y[i] = (f32)rand() / (f32)RAND_MAX;
		z[i] = sign * (f32)rand() / (f32)RAND_MAX;

		//
		vx[i] = (f32)rand() / (f32)RAND_MAX;
		vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
		vz[i] = (f32)rand() / (f32)RAND_MAX;
	}
	for (u64 i  = 0  ; i < 64 ; i++){
		softening[i] = softening_value;
		dt[i] = dt_value ; 
	}
}

//
void move_particles(float *x, float *y, float *z, float *vx, float *vy, float *vz, float *softening, float *dt, u64 n)
{
	//
	__m256 rxi, ryi, rzi, rxj, ryj, rzj , rfx, rfy, rfz;
	__m256 rvxi, rvyi, rvzi ; 
	__m256 rdx  , rdy , rdz , rdxyz , tmp, tmp2, rsoft, rt;
	rxi = _mm256_setzero_ps() ;
	ryi = _mm256_setzero_ps() ;	
	rzi = _mm256_setzero_ps() ;	
	rvxi = _mm256_setzero_ps() ;	
	rvyi = _mm256_setzero_ps() ;	
	rvzi = _mm256_setzero_ps() ;	
	
	rxj = _mm256_setzero_ps() ;
	ryj = _mm256_setzero_ps() ;	
	rzj = _mm256_setzero_ps() ;	

	
	rdx = _mm256_setzero_ps() ;
	rdy = _mm256_setzero_ps() ;	
	rzj = _mm256_setzero_ps() ;	
	rdz = _mm256_setzero_ps() ;	
	rdxyz = _mm256_setzero_ps() ;	
	tmp = _mm256_setzero_ps() ;		
	tmp2 = _mm256_setzero_ps() ;
	rsoft = _mm256_load_ps(&softening[0]) ;
	rt = _mm256_load_ps(&dt[0]) ;	
	for (u64 i = 0; i < n; i++)
	{
	//
		rfx = _mm256_setzero_ps() ;	
		rfy = _mm256_setzero_ps() ;	
		rfz = _mm256_setzero_ps() ;

		// load i values
		rxi = _mm256_loadu_ps(&x[i]);
		ryi = _mm256_loadu_ps(&y[i]);
		rzi = _mm256_loadu_ps(&z[i]);			
		
		for (u64 j = 0; j < n; j+=8)
		{
			//printf("j : %lld\n", j);
			
			rxj = _mm256_load_ps(&x[j]);
			ryj = _mm256_load_ps(&y[j]);
			rzj = _mm256_load_ps(&z[j]);
			rdx = _mm256_sub_ps(rxj, rxi) ; 
			rdy = _mm256_sub_ps(ryj, ryi) ; 
			rdz = _mm256_sub_ps(rzj, rzi) ; 
			tmp = _mm256_sub_ps(rdx, rdx) ;
			tmp2 = _mm256_sub_ps(rdz, rdz) ;
			rdxyz = _mm256_add_ps(tmp, tmp2) ;
			tmp =   _mm256_sub_ps(rdz, rdz) ;
			rdxyz = _mm256_add_ps(tmp, rdxyz) ;	
			rdxyz = _mm256_add_ps(rsoft, rdxyz) ;	
			tmp = _mm256_sqrt_ps(rdxyz) ;
			rdxyz = _mm256_mul_ps(tmp, tmp) ;	
			rdxyz = _mm256_mul_ps(rdxyz, tmp) ;	
			
			rdx = _mm256_div_ps(rdx, rdxyz) ; 			
			rdy = _mm256_div_ps(rdy, rdxyz) ; 
			rdz = _mm256_div_ps(rdz, rdxyz) ; 			
			
			rfx = _mm256_add_ps(rfx, rdx) ; 	
			rfy = _mm256_add_ps(rfy, rdy) ; 	
			rfz = _mm256_add_ps(rfz, rdz) ; 				
		}
		
		rfx = _mm256_mul_ps(rfx, rt) ; 	
		rfy = _mm256_mul_ps(rfy, rt) ; 	
		rfz = _mm256_mul_ps(rfz, rt) ; 		
		
		rvxi = _mm256_loadu_ps(&vx[i]);
		rvyi = _mm256_loadu_ps(&vy[i]);
		rvzi = _mm256_loadu_ps(&vz[i]);
		
		rvxi = _mm256_add_ps(rvxi, rfx) ; 	
		rvyi = _mm256_add_ps(rvyi, rfy) ; 
		rvzi = _mm256_add_ps(rvzi, rfz) ; 
		
		_mm256_storeu_ps(&vx[i], rvxi) ;				
		_mm256_storeu_ps(&vy[i], rvyi) ;	
		_mm256_storeu_ps(&vz[i], rvzi) ;	

	}

	//3 floating-point operations
	for (u64 i = 0; i < n; i++)
	{
	x[i] += dt[0] * vx[i];
	y[i] += dt[0] * vy[i];
	z[i] += dt[0] * vz[i];
	}
}

//
int main(int argc, char **argv)
{
	//
	const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
	const u64 steps= 10;
	const f32 dt_value = 0.01;
	const f32 softening_value = 1e-20;
	//
	f64 rate = 0.0, drate = 0.0;

	//Steps to skip for warm up
	const u64 warmup = 3;

	//
	//particle_t *p = malloc(sizeof(particle_t) * n);
	float *x = aligned_alloc(64, sizeof(float) * n) ; 
	float *y = aligned_alloc(64, sizeof(float) * n) ; 	
	float *z = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vx = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vy = aligned_alloc(64, sizeof(float) * n) ; 	
	float *vz = aligned_alloc(64, sizeof(float) * n) ; 	
	float *softening = aligned_alloc(64, sizeof(float) * 64) ; 
	float *dt = aligned_alloc(64, sizeof(float) * 64) ; 	
	//
	init(x, y,z, vx, vy, vz, softening, softening_value, dt, dt_value,  n);

	const u64 sx = 6 *  sizeof(double) * n;

	printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", sx, sx >> 10, sx >> 20);

	//
	printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

	//
	for (u64 i = 0; i < steps; i++)
	{
	//Measure
	const f64 start = omp_get_wtime();

	move_particles(x, y, z, vx, vy, vz, softening, dt, n);

	const f64 end = omp_get_wtime();

	//Number of interactions/iterations
	const f32 h1 = (f32)(n) * (f32)(n - 1);

	//GFLOPS
	const f32 h2 = (23.0 * h1 + 3.0 * (f32)n) * 1e-9;

	if (i >= warmup)
	{
		rate += h2 / (end - start);
		drate += (h2 * h2) / ((end - start) * (end - start));
	}

	//
	printf("%5llu %10.3e %10.3e %8.1f %s\n",
	i,
	(end - start),
	h1 / (end - start),
	h2 / (end - start),
	(i < warmup) ? "*" : "");

	fflush(stdout);
	}

	//
	rate /= (f64)(steps - warmup);
	drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

	printf("-----------------------------------------------------\n");
	printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	"Average performance:", "", rate, drate);
	printf("-----------------------------------------------------\n");

	//
	free(x);
	free(y);
	free(z);
	free(vx);
	free(vy);
	free(vz);

	//
	return 0;
}
