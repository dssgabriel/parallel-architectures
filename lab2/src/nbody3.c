#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <SDL2/SDL.h>
#include "types.h"

typedef struct {
    f64 *x;
    f64 *y;
} vectors;

static const usize width = 800;
static const usize height = 800;
static const usize nbodies = 500;
static const usize time_step = 1000;
static const f64 mass = 5.0;
static const f64 grav_cst = 1.0;
static const f64 cst = grav_cst * mass;
static usize hit = 0;

vectors positions, velocities, accelerations;

static inline
usize rdtsc()
{
    usize a, d;

    __asm__ volatile(
        "rdtsc"
        : "=a" (a), "=d" (d));

    return (d << 32) | a;
}

i32 randxy(i32 x, i32 y)
{
    return (rand() % (y - x + 1)) + x; 
}

double randreal()
{
    i32 s = (randxy(0, 1)) ? 1 : -1;
    i32 a = randxy(1, RAND_MAX), b = randxy(1, RAND_MAX);
    return s * ((f64)a / (f64)b); 
}

void init_system()
{
    positions.x = aligned_alloc(64, nbodies * sizeof(vectors));
    positions.y = aligned_alloc(64, nbodies * sizeof(vectors));
    velocities.x = aligned_alloc(64, nbodies * sizeof(vectors));
    velocities.y = aligned_alloc(64, nbodies * sizeof(vectors));
    accelerations.x = aligned_alloc(64, nbodies * sizeof(vectors));
    accelerations.y = aligned_alloc(64, nbodies * sizeof(vectors));

    for (usize i = 0; i < nbodies; i++) {
        positions.x[i] = randxy(10, width);
        positions.y[i] = randxy(10, height);

        velocities.x[i] = randreal();
        velocities.y[i] = randreal();
    }
}

void deinit_system()
{
    free(positions.x);
    free(positions.y);
    free(velocities.x);
    free(velocities.y);
    free(accelerations.x);
    free(accelerations.y);
}

void compute_accelerations()
{ 
    for (usize i = 0; i < nbodies; i++) {
        accelerations.x[i] = 0;
        accelerations.y[i] = 0;

        f64 __attribute__((aligned(64))) posx[8];
        f64 __attribute__((aligned(64))) posy[8];
        f64 __attribute__((aligned(64))) modv[8];
        f64 __attribute__((aligned(64))) dist[8];

        for (usize j = i + 1; j < nbodies; j += 8) {
            posx[0] = positions.x[j] - positions.x[i];
            posx[1] = positions.x[j + 1] - positions.x[i];
            posx[2] = positions.x[j + 2] - positions.x[i];
            posx[3] = positions.x[j + 3] - positions.x[i];
            posx[4] = positions.x[j + 4] - positions.x[i];
            posx[5] = positions.x[j + 5] - positions.x[i];
            posx[6] = positions.x[j + 6] - positions.x[i];
            posx[7] = positions.x[j + 7] - positions.x[i];
            
            posy[0] = positions.y[j] - positions.y[i];
            posy[1] = positions.y[j + 1] - positions.y[i];
            posy[2] = positions.y[j + 2] - positions.y[i];
            posy[3] = positions.y[j + 3] - positions.y[i];
            posy[4] = positions.y[j + 4] - positions.y[i];
            posy[5] = positions.y[j + 5] - positions.y[i];
            posy[6] = positions.y[j + 6] - positions.y[i];
            posy[7] = positions.y[j + 7] - positions.y[i];

            modv[0] = sqrt(posx[0] * posx[0] + posy[0] * posy[0]);
            modv[1] = sqrt(posx[1] * posx[1] + posy[1] * posy[1]);
            modv[2] = sqrt(posx[2] * posx[2] + posy[2] * posy[2]);
            modv[3] = sqrt(posx[3] * posx[3] + posy[3] * posy[3]);
            modv[4] = sqrt(posx[4] * posx[4] + posy[4] * posy[4]);
            modv[5] = sqrt(posx[5] * posx[5] + posy[5] * posy[5]);
            modv[6] = sqrt(posx[6] * posx[6] + posy[6] * posy[6]);
            modv[7] = sqrt(posx[7] * posx[7] + posy[7] * posy[7]);
            
            dist[0] = cst / (modv[0] * modv[0] * modv[0] + 1e7);
            dist[1] = cst / (modv[1] * modv[1] * modv[1] + 1e7);
            dist[2] = cst / (modv[2] * modv[2] * modv[2] + 1e7);
            dist[3] = cst / (modv[3] * modv[3] * modv[3] + 1e7);
            dist[4] = cst / (modv[4] * modv[4] * modv[4] + 1e7);
            dist[5] = cst / (modv[5] * modv[5] * modv[5] + 1e7);
            dist[6] = cst / (modv[6] * modv[6] * modv[6] + 1e7);
            dist[7] = cst / (modv[7] * modv[7] * modv[7] + 1e7);
        }

        for (usize j = i + 1; j < nbodies; j += 8) {
            accelerations.x[i] += (dist[0] * posx[0]) + (dist[1] * posx[1]) +
                                  (dist[2] * posx[2]) + (dist[3] * posx[3]) + 
                                  (dist[4] * posx[4]) + (dist[5] * posx[5]) +
                                  (dist[6] * posx[6]) + (dist[7] * posx[7]);
            accelerations.y[i] += (dist[0] * posy[0]) + (dist[1] * posy[1]) +
                                  (dist[2] * posy[2]) + (dist[3] * posy[3]) +
                                  (dist[4] * posy[4]) + (dist[5] * posy[5]) +
                                  (dist[6] * posy[6]) + (dist[7] * posy[7]);

            accelerations.x[j] = -accelerations.x[i];
            accelerations.y[j] = -accelerations.y[i];
        }
    }
}

void compute_velocities()
{  
    for (usize i = 0; i < nbodies; i++) {
        velocities.x[i] += accelerations.x[i];
        velocities.y[i] += accelerations.y[i];
    }
}

void compute_positions()
{
    for (usize i = 0; i < nbodies; i++) {
        positions.x[i] += (velocities.x[i] + 0.5 * accelerations.x[i]);
        positions.y[i] += (velocities.y[i] + 0.5 * accelerations.y[i]);
    }
}

void simulate()
{
    compute_accelerations();
    compute_positions();
    compute_velocities();
}

int headless()
{
    u8 quit = 0;
    usize before, after;

    srand(getpid());
    init_system();

    // Main loop
    for (usize i = 0; !quit && i < time_step; i++) {	  
        before = rdtsc();
        simulate();
        after = rdtsc();

        printf("%zu %zu\n", i, (after - before));
    }

    deinit_system();

    return 0;
}

int visual()
{
    u8 quit = 0;
    usize before, after;

    SDL_Event event;
    SDL_Window *window;
    SDL_Renderer *renderer;

    srand(getpid());
    init_system();
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(width, height,SDL_WINDOW_OPENGL,
                                &window, &renderer);

    // Main loop
    for (usize i = 0; !quit && i < time_step; i++) {	  
        before = rdtsc();
        simulate();
        after = rdtsc();

        printf("%zu %zu\n", i, (after - before));

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        for (usize j = 0; j < nbodies; j++) {
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderDrawPoint(renderer, positions.x[j], positions.y[j]);
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(10);
        while (SDL_PollEvent(&event)) {
        	if (event.type == SDL_QUIT) {
        	    quit = 1;
        	} else if (event.type == SDL_KEYDOWN) {
        	    if (event.key.keysym.sym == SDLK_q)
        	        quit = 1;
        	}
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    deinit_system();
    printf("%zu\n", hit);

    return 0;
}

int main()
{
#define SDL 0
#if SDL == 0 
    return headless();
#else
    return visual();
#endif
}
