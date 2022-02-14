#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <SDL2/SDL.h>
#include "types.h"

typedef struct {
    f64 *x;
    f64 *y;
} Vec2D;

static const usize width = 800;
static const usize height = 800;
static const usize nbodies = 500;
static const usize time_step = 1000;
static const f64 mass = 5.0;
static const f64 grav_cst = 1.0;

Vec2D positions, velocities, accelerations;

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

f64 randreal()
{
    i32 s = (randxy(0, 1)) ? 1 : -1;
    i32 a = randxy(1, RAND_MAX), b = randxy(1, RAND_MAX);
    return s * ((f64)a / (f64)b); 
}

void init_system()
{
    positions.x = malloc(nbodies * sizeof(f32));
    positions.y = malloc(nbodies * sizeof(f32));

    velocities.x = malloc(nbodies * sizeof(f32));
    velocities.y = malloc(nbodies * sizeof(f32));

    accelerations.x = malloc(nbodies * sizeof(f32));
    accelerations.y = malloc(nbodies * sizeof(f32));

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
        accelerations.x[i] = 0.0f;
        accelerations.y[i] = 0.0f;

        for (usize j = i + 1; j < nbodies; j++) {
        	if (i != j) {
                const f32 dx = positions.x[i] - positions.x[j];
                const f32 dy = positions.y[i] - positions.y[j];
                const f32 d_2 = (dx * dx) + (dy * dy);
                const f32 d_2_sqrt = sqrtf(d_2);
                const f32 d_2_rsqrt = 1.0f / d_2_sqrt;
                const f32 d_3_over_2 = (d_2_rsqrt * d_2_rsqrt) * d_2_rsqrt + 1e7;
                const f32 f = d_3_over_2 * mass;

            	accelerations.x[i] += f * dx;
            	accelerations.y[i] += f * dy;

            	accelerations.x[j] = -accelerations.x[i];
            	accelerations.y[j] = -accelerations.y[i];
        	}
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
        positions.x[i] += velocities.x[i] + (0.5 * accelerations.x[i]);
        positions.y[i] += velocities.y[i] + (0.5 * accelerations.y[i]);
    }
}

void simulate()
{
    compute_accelerations();
    compute_velocities();
    compute_positions();
}

int main()
{
    u8 quit = 0;
    usize before, after;

    // SDL_Event event;
    // SDL_Window *window;
    // SDL_Renderer *renderer;

    srand(getpid());
    init_system();
    // SDL_Init(SDL_INIT_VIDEO);
    // SDL_CreateWindowAndRenderer(width, height,SDL_WINDOW_OPENGL,
    //                             &window, &renderer);

    // Main loop
    for (usize i = 0; !quit && i < time_step; i++) {	  
        before = rdtsc();
        simulate();
        after = rdtsc();

        fprintf(stderr, "%zu %zu\n", i, (after - before));

    //     SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    //     SDL_RenderClear(renderer);

    //     for (usize j = 0; j < nbodies; j++) {
    //         SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    //         SDL_RenderDrawPoint(renderer, positions.x[j], positions.y[j]);
    //     }

    //     SDL_RenderPresent(renderer);
    //     SDL_Delay(10);
    //     while (SDL_PollEvent(&event)) {
    //     	if (event.type == SDL_QUIT) {
    //     	    quit = 1;
    //     	} else if (event.type == SDL_KEYDOWN) {
    //     	    if (event.key.keysym.sym == SDLK_q)
    //     	        quit = 1;
    //     	}
    //     }
    }

    // SDL_DestroyRenderer(renderer);
    // SDL_DestroyWindow(window);
    // SDL_Quit();
    // deinit_system();

    return 0;
}
