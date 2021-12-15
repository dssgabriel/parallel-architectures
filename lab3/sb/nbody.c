//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

  f32 x, y, z;
  f32 vx, vy, vz;
  
} particle_t;

//
void init(particle_t *p, u64 n)
{
  for (u64 i = 0; i < n; i++)
    {
      //
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p[i].x = sign * (f32)rand() / (f32)RAND_MAX;
      p[i].y = (f32)rand() / (f32)RAND_MAX;
      p[i].z = sign * (f32)rand() / (f32)RAND_MAX;

      //
      p[i].vx = (f32)rand() / (f32)RAND_MAX;
      p[i].vy = sign * (f32)rand() / (f32)RAND_MAX;
      p[i].vz = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t *p, const f32 dt, u64 n)
{
  //
  const f32 softening = 1e-20;

  //
  for (u64 i = 0; i < n; i++)
    {
      //
      f32 fx = 0.0;
      f32 fy = 0.0;
      f32 fz = 0.0;

      //23 floating-point operations
      for (u64 j = 0; j < n; j++)
	{
	  //Newton's law
	  const f32 dx = p[j].x - p[i].x; //1
	  const f32 dy = p[j].y - p[i].y; //2
	  const f32 dz = p[j].z - p[i].z; //3
	  const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9
	  const f32 d_3_over_2 = pow(d_2, 3.0 / 2.0); //11

	  //Net force
	  fx += dx / d_3_over_2; //13
	  fy += dy / d_3_over_2; //15
	  fz += dz / d_3_over_2; //17
	}

      //
      p[i].vx += dt * fx; //19
      p[i].vy += dt * fy; //21
      p[i].vz += dt * fz; //23
    }

  //3 floating-point operations
  for (u64 i = 0; i < n; i++)
    {
      p[i].x += dt * p[i].vx;
      p[i].y += dt * p[i].vy;
      p[i].z += dt * p[i].vz;
    }
}

//
int main(int argc, char **argv)
{
  //
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 steps= 10;
  const f32 dt = 0.01;

  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t *p = malloc(sizeof(particle_t) * n);

  //
  init(p, n);

  const u64 s = sizeof(particle_t) * n;
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  //
  for (u64 i = 0; i < steps; i++)
    {
      //Measure
      const f64 start = omp_get_wtime();

      move_particles(p, dt, n);

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
  free(p);

  //
  return 0;
}
