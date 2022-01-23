#include <cassert>
#include <cmath>
#include <cstdio>
#include <mpi.h>
#include <omp.h>
#include <time.h>

#include "types.h"
#include "utils.h"

struct particles_t {
    f32 *px, *py, *pz;
    f32 *vx, *vy, *vz;
};

void particles_update(particles_t &p,
                      const u64 nb_particles,
                      const f32 dt,
		              const int mpi_rank,
                      const int mpi_world_size,
                      const int particles)
{
    const size_t TILE_SIZE = 16;
    const int start_particle = mpi_rank * particles;
    const int end_particle = (mpi_rank + 1) * particles;

    assert(particles % TILE_SIZE == 0);

    #pragma omp parallel for
    for (size_t t = start_particle; t < end_particle; t += TILE_SIZE) {
        f32 fx[TILE_SIZE], fy[TILE_SIZE], fz[TILE_SIZE]; 
        fx[:] = fy[:] = fz[:] = 0.0f;

        const f32 softening = 1e-20;

        #pragma unroll(TILE_SIZE)
        for (size_t j = 0; j < nb_particles; j++) {
            #pragma vector aligned
            for (size_t i = t; i < t + TILE_SIZE; i++) {
                const f32 dx = p.px[j] - p.px[i];
                const f32 dy = p.py[j] - p.py[i];
                const f32 dz = p.pz[j] - p.pz[i];
                const f32 dr_squared = ((dx * dx) + (dy * dy)) + ((dz * dz) + softening);

                const f32 dr_sqrted = sqrtf(dr_squared);
                const f32 inv_dr = 1.0f / dr_sqrted;
                const f32 inv_dr_p3 = inv_dr * inv_dr * inv_dr;

                fx[i - t] += dx * inv_dr_p3;
                fy[i - t] += dy * inv_dr_p3;
                fz[i - t] += dz * inv_dr_p3;
            }
        }

        p.vx[t : TILE_SIZE] += dt * fx [0 : TILE_SIZE]; 
        p.vy[t : TILE_SIZE] += dt * fy [0 : TILE_SIZE]; 
        p.vz[t : TILE_SIZE] += dt * fz [0 : TILE_SIZE];
    }

    #pragma omp parallel for
    #pragma vector aligned
    for (size_t i = 0 ; i < nb_particles; i++) { 
        p.px[i] += p.vx[i] * dt;
        p.py[i] += p.vy[i] * dt;
        p.pz[i] += p.vz[i] * dt;
    }

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.px, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.py, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.pz, particles, MPI_FLOAT, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init((int*)&argc, (char***)&argv);
    int mpi_world_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const u64 nb_particles = (argc > 1 ? atoi(argv[1]) : 16384);

    assert(nb_particles % mpi_world_size == 0);

    const int particles = nb_particles / mpi_world_size;
    const int nb_iter = 10;
    const int nb_warm = 3; 
    const f32 dt = 0.01f;
    f64 rate = 0.0f, drate = 0.0f, duration = 0.0f, dduration = 0.0f;

    particles_t p;
    p.px = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);
    p.py = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);
    p.pz = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);
    p.vx = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);
    p.vy = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);
    p.vz = (f32*)_mm_malloc(nb_particles * sizeof(f32), 64);

    srand(0);
    for (size_t i = 0; i < nb_particles; i++) {
        p.px[i] = rand() / RAND_MAX;
        p.py[i] = rand() / RAND_MAX;
        p.pz[i] = rand() / RAND_MAX;
        p.vx[i] = rand() / RAND_MAX;
        p.vy[i] = rand() / RAND_MAX;
        p.vz[i] = rand() / RAND_MAX;
    }
  
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.px, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.py, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.pz, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.vx, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.vy, particles, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, p.vz, particles, MPI_FLOAT, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        const u32 num_threads = omp_get_max_threads();
        printf("\nPropagating %llu particles using %d threads, distributed across %d nodes on %s...\n\n",
               nb_particles, num_threads, mpi_world_size, "CPU");
        printf("\033[1m%5s %10s %10s %8s\033[0m\n",
               "Iter", "Time (s)", "Interact/s", "GFLOP/s");
        fflush(stdout);
    }

    for (size_t i = 1; i <= nb_iter; i++) {
        const f64 start = omp_get_wtime();
        particles_update(p, nb_particles, dt, mpi_rank, mpi_world_size, particles);
        const f64 end = omp_get_wtime();

        const f64 h1= (f64)(nb_particles) * (f64)(nb_particles - 1);
        const f64 h2 = 20.0f * 1e-9 * (f64)(nb_particles) * (f64)(nb_particles - 1);

        if (nb_iter > nb_warm) {
            rate += h1 / (end - start);
            drate += h2 * h2 / ((end - start) * (end - start)); 
            duration += (end - start);
            dduration += (end - start) * (end - start);
        }

        if (mpi_rank == 0) {
            printf("%5d %10.3e %10.3e %8.1f %s\n", 
                   nb_iter, (end - start),
                   h1 / (end - start),
                   h2 / (end - start),
                   (nb_iter <= nb_warm ? "*" : ""));
        }

        fflush(stdout);
  }

    rate /= (f64)(nb_iter - nb_warm); 
    drate = sqrtl(drate / (f64)(nb_iter - nb_warm) - (rate * rate));
    duration /= (f64)(nb_iter - nb_warm);
    dduration = sqrtl(dduration / (f64)(nb_iter - nb_warm) - (duration * duration));

    if (mpi_rank == 0) {
        printf("-----------------------------------------------------\n");
        printf("\033[1m%s %4s \033[42m%10.1f +- %.1f GFLOP/s\033[0m\n",
               "Average performance:", "", rate, drate);
        printf("\033[1m%s %4s \033[42m%10.3e +- %10.3e s\033[0m\n",
               "Average duration:", "", duration, dduration);
        printf("\033[1m%s %4s \033[42m%.3f +- %.3f ms\033[0m\n",
               "Average duration:", "", duration * 1000.0, dduration * 1000.0);
        printf("-----------------------------------------------------\n");
    }

    _mm_free(p.px);
    _mm_free(p.py);
    _mm_free(p.pz);
    _mm_free(p.vx);
    _mm_free(p.vy);
    _mm_free(p.vz);

    MPI_Finalize();
}
