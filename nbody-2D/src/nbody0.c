#include "types.h"

#include <math.h>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    double x;
    double y;
} vector;

static const size_t width = 800;
static const size_t height = 800;
static const size_t nbodies = 500;
static const size_t time_step = 1000;
static const double mass = 5.0;
static const double grav_cst = 1.0;

vector *positions, *velocities, *accelerations;

static inline size_t rdtsc() {
    size_t a, d;

    __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));

    return (d << 32) | a;
}

int32_t randxy(int32_t x, int32_t y) {
    return (rand() % (y - x + 1)) + x;
}

double randreal() {
    int32_t s = (randxy(0, 1)) ? 1 : -1;
    int32_t a = randxy(1, RAND_MAX), b = randxy(1, RAND_MAX);
    return s * ((double)a / (double)b);
}

vector add_vectors(vector a, vector b) {
    vector c = { a.x + b.x, a.y + b.y };
    return c;
}

vector scale_vector(double b, vector a) {
    vector c = { b * a.x, b * a.y };
    return c;
}

vector sub_vectors(vector a, vector b) {
    vector c = { a.x - b.x, a.y - b.y };
    return c;
}

double mod(vector a) {
    return sqrt(a.x * a.x + a.y * a.y);
}

void init_system() {
    positions = malloc(nbodies * sizeof(vector));
    velocities = malloc(nbodies * sizeof(vector));
    accelerations = malloc(nbodies * sizeof(vector));

    for (size_t i = 0; i < nbodies; i++) {
        positions[i].x = randxy(10, width);
        positions[i].y = randxy(10, height);

        velocities[i].x = randreal();
        velocities[i].y = randreal();
    }
}

void deinit_system() {
    free(positions);
    free(velocities);
    free(accelerations);
}

void resolve_collisions() {
    for (size_t i = 0; i < nbodies - 1; i++) {
        for (size_t j = i + 1; j < nbodies; j++) {
            if (positions[i].x == positions[j].x && positions[i].y == positions[j].y) {
                vector tmp = velocities[i];
                velocities[i] = velocities[j];
                velocities[j] = tmp;
            }
        }
    }
}

void compute_accelerations() {
    for (size_t i = 0; i < nbodies; i++) {
        accelerations[i].x = 0;
        accelerations[i].y = 0;

        for (size_t j = 0; j < nbodies; j++) {
            if (i != j)
                accelerations[i] = add_vectors(
                    accelerations[i],
                    scale_vector(grav_cst * mass /
                                     (pow(mod(sub_vectors(positions[i], positions[j])), 3) + 1e7),
                                 sub_vectors(positions[j], positions[i])));
        }
    }
}

void compute_velocities() {
    for (size_t i = 0; i < nbodies; i++)
        velocities[i] = add_vectors(velocities[i], accelerations[i]);
}

void compute_positions() {
    for (size_t i = 0; i < nbodies; i++)
        positions[i] = add_vectors(positions[i],
                                   add_vectors(velocities[i], scale_vector(0.5, accelerations[i])));
}

void simulate() {
    compute_accelerations();
    compute_positions();
    compute_velocities();
    resolve_collisions();
}

int main() {
    uint8_t quit = 0;
    size_t before, after;

    SDL_Event event;
    SDL_Window *window;
    SDL_Renderer *renderer;

    srand(getpid());
    init_system();
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(width, height, SDL_WINDOW_OPENGL, &window, &renderer);

    // Main loop
    for (size_t i = 0; !quit && i < time_step; i++) {
        before = rdtsc();
        simulate();
        after = rdtsc();

        printf("%zu %zu\n", i, (after - before));

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        for (size_t j = 0; j < nbodies; j++) {
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderDrawPoint(renderer, positions[j].x, positions[j].y);
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(10);
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
            else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_q)
                    quit = 1;
            }
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    deinit_system();

    return 0;
}
