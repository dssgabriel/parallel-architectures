#include "types.h"

#include <immintrin.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define BOLD "\033[1m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define RESET "\033[0m"

#define DEFAULT_DIMENSIONS 800
#define DEFAULT_BODIES 500
#define DEFAULT_ITER 1000
#define DEFAULT_MASS 1.0
#define DEFAULT_GRAVITY 1.0

#define CORRECTION 1e7
#define FP32_IN_ZMM (512 / (8 * sizeof(float)))
#define MAX_RGB 255
#define ALIGNMENT 64

void print_help(char const *bin) {
    printf(BOLD "%s - 2D N-body simulation\n" RESET, bin);
    printf(BOLD "\nUsage: " RESET "%s [FLAGS] [OPTIONS]\n", bin);
    printf(BOLD "\nFlags:\n" RESET);
    printf(GREEN "    -h    \t" RESET "Prints this help.\n");
    printf(GREEN "    -s    \t" RESET "Enables SDL output.\n");
    printf(BOLD "\nOptions:\n" RESET);
    printf(GREEN "    -d <N>\t" RESET "Dimensions in pixels (default = %ux%u).\n",
           DEFAULT_DIMENSIONS, DEFAULT_DIMENSIONS);
    printf(GREEN "    -n <N>\t" RESET "Number of particules (default = %u).\n", DEFAULT_BODIES);
    printf(GREEN "    -i <N>\t" RESET "Number of iterations (default = %u).\n", DEFAULT_ITER);
    printf(GREEN "    -m <N>\t" RESET "Mass of the particules (default = %f).\n", DEFAULT_MASS);
    printf(GREEN "    -g <N>\t" RESET "Gravity constant (default = %f).\n", DEFAULT_GRAVITY);
}

typedef struct {
    float *x;
    float *y;
    size_t len;
} Vec2D;

typedef struct {
    size_t dimensions;
    size_t bodies;
    size_t iterations;
    float mass;
    float gravity;
    bool has_SDL;
} Config;

Config config_default() {
    Config cfg = {
        .dimensions = DEFAULT_DIMENSIONS,
        .bodies = DEFAULT_BODIES,
        .iterations = DEFAULT_ITER,
        .mass = DEFAULT_MASS,
        .gravity = DEFAULT_GRAVITY,
        .has_SDL = false,
    };
    return cfg;
}

Config config_from_args(int32_t const argc, char *argv[argc + 1]) {
    Config cfg = config_default();

    int32_t current;
    while ((current = getopt(argc, argv, "d:n:i:m:g:sh")) != -1) {
        switch (current) {
            case 'd':
                cfg.dimensions = (size_t)(atoi(optarg));
                break;
            case 'n':
                cfg.bodies = (size_t)(atoi(optarg));
                break;
            case 'i':
                cfg.iterations = (size_t)(atoi(optarg));
                break;
            case 'm':
                cfg.mass = (double)(atof(optarg));
                break;
            case 'g':
                cfg.gravity = (double)(atof(optarg));
                break;
            case 's':
                cfg.has_SDL = true;
                break;
            case 'h':
                print_help(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                fprintf(stderr, BOLD RED "error: " RESET "unknown argument %s.\n", optarg);
                print_help(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    return cfg;
}

static inline uint64_t rdtsc() {
    uint64_t a;
    uint64_t d;
    __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));
    return (d << 32) | a;
}

int32_t rand_i32_in_range(int32_t min, int32_t max) {
    return (random() % (max - min + 1)) + min;
}

double rand_f64() {
    int32_t sign = rand_i32_in_range(0, 1) ? 1 : -1;
    int32_t a = rand_i32_in_range(1, RAND_MAX);
    int32_t b = rand_i32_in_range(1, RAND_MAX);
    return sign * ((double)(a) / (double)(b));
}

void init_system(Vec2D *pos, Vec2D *vel, Vec2D *acc, size_t bodies, size_t dimensions) {
    pos->x = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    pos->y = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    pos->len = bodies;

    vel->x = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    vel->y = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    vel->len = bodies;

    acc->x = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    acc->y = aligned_alloc(ALIGNMENT, bodies * sizeof(float));
    acc->len = bodies;

    for (size_t i = 0; i < bodies; ++i) {
        pos->x[i] = (float)(rand_i32_in_range(1, dimensions));
        pos->y[i] = (float)(rand_i32_in_range(1, dimensions));
        vel->x[i] = (float)(rand_f64());
        vel->y[i] = (float)(rand_f64());
    }
}

void deinit_system(Vec2D *pos, Vec2D *vel, Vec2D *acc) {
    free(pos->x);
    free(pos->y);
    free(vel->x);
    free(vel->y);
    free(acc->x);
    free(acc->y);
}

void resolve_collisions(Vec2D const *pos, Vec2D *vel) {
    for (size_t i = 0; i < pos->len - 1; ++i) {
        for (size_t j = i + 1; j < pos->len; ++j) {
            if (pos->x[i] == pos->x[j] && pos->y[i] == pos->y[j]) {
                float tmp_x = vel->x[i];
                float tmp_y = vel->y[i];
                vel->x[i] = vel->x[j];
                vel->y[i] = vel->y[j];
                vel->x[j] = tmp_x;
                vel->y[j] = tmp_y;
            }
        }
    }
}

void compute_acc(Vec2D *acc, Vec2D const *pos, float mass, float gravity) {
    float force = mass * gravity;
    for (size_t i = 0; i < acc->len; ++i) {
        acc->x[i] = 0.0;
        acc->y[i] = 0.0;
        for (size_t j = 0; j < acc->len; ++j) {
            if (i != j) {
                float dx = pos->x[j] - pos->x[i];
                float dy = pos->y[j] - pos->y[i];
                float dsq = sqrtf(dx * dx + dy * dy);
                dsq = dsq * dsq * dsq + CORRECTION;
                acc->x[i] += force / dsq * (dx);
                acc->y[i] += force / dsq * (dy);
            }
        }
    }
}

void compute_vel(Vec2D *vel, Vec2D const *acc) {
    size_t len = vel->len - (vel->len % FP32_IN_ZMM);
    for (size_t i = 0; i < len; i += FP32_IN_ZMM) {
        __m512 vx = _mm512_loadu_ps(&vel->x[i]);
        __m512 vy = _mm512_loadu_ps(&vel->y[i]);
        __m512 ax = _mm512_loadu_ps(&acc->x[i]);
        __m512 ay = _mm512_loadu_ps(&acc->y[i]);
        vx = _mm512_add_ps(vx, ax);
        vy = _mm512_add_ps(vy, ay);
        _mm512_storeu_ps(&vel->x[i], vx);
        _mm512_storeu_ps(&vel->y[i], vy);
    }

    for (size_t i = len; i < vel->len; ++i) {
        vel->x[i] += acc->x[i];
        vel->y[i] += acc->y[i];
    }
}

void compute_pos(Vec2D *pos, Vec2D const *vel, Vec2D const *acc) {
    __m512 half = _mm512_set1_ps(0.5);
    size_t len = pos->len - (pos->len % FP32_IN_ZMM);
    for (size_t i = 0; i < len; i += FP32_IN_ZMM) {
        __m512 px = _mm512_loadu_ps(&pos->x[i]);
        __m512 py = _mm512_loadu_ps(&pos->y[i]);
        __m512 vx = _mm512_loadu_ps(&vel->x[i]);
        __m512 vy = _mm512_loadu_ps(&vel->y[i]);
        __m512 ax = _mm512_loadu_ps(&acc->x[i]);
        __m512 ay = _mm512_loadu_ps(&acc->y[i]);
        px = _mm512_add_ps(px, _mm512_fmadd_ps(half, ax, vx));
        py = _mm512_add_ps(py, _mm512_fmadd_ps(half, ay, vy));
        _mm512_storeu_ps(&pos->x[i], px);
        _mm512_storeu_ps(&pos->y[i], py);
    }

    for (size_t i = len; i < pos->len; ++i) {
        pos->x[i] += vel->x[i] + 0.5 * acc->x[i];
        pos->y[i] += vel->y[i] + 0.5 * acc->y[i];
    }
}

void simulate_with_graphics(Config const *cfg, Vec2D *pos, Vec2D *vel, Vec2D *acc) {
    uint8_t quit = 0;
    SDL_Event event;
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(cfg->dimensions, cfg->dimensions, SDL_WINDOW_OPENGL, &window,
                                &renderer);

    for (size_t i = 0; !quit && i < cfg->iterations; ++i) {
        uint64_t before = rdtsc();
        compute_acc(acc, pos, cfg->mass, cfg->gravity);
        compute_vel(vel, acc);
        compute_pos(pos, vel, acc);
        resolve_collisions(pos, vel);
        uint64_t after = rdtsc();
        printf("%lu %lu\n", i + 1, (after - before));

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, MAX_RGB);
        SDL_RenderClear(renderer);
        for (size_t j = 0; j < cfg->bodies; ++j) {
            SDL_SetRenderDrawColor(renderer, MAX_RGB, MAX_RGB, MAX_RGB, MAX_RGB);
            SDL_RenderDrawPoint(renderer, pos->x[j], pos->y[j]);
        }
        SDL_RenderPresent(renderer);
        // SDL_Delay(10);
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
            else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_q) {
                    quit = 1;
                }
            }
        }
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void simulate(Config const *cfg, Vec2D *pos, Vec2D *vel, Vec2D *acc) {
    for (size_t i = 0; i < cfg->iterations; ++i) {
        uint64_t before = rdtsc();
        compute_acc(acc, pos, cfg->mass, cfg->gravity);
        compute_vel(vel, acc);
        compute_pos(pos, vel, acc);
        resolve_collisions(pos, vel);
        uint64_t after = rdtsc();
        printf("%lu %lu\n", i + 1, (after - before));
    }
}

int32_t main(int32_t argc, char *argv[argc + 1]) {
    printf("%lu\n", FP32_IN_ZMM);
    Config cfg;
    if (argc < 2) {
        cfg = config_default();
    }
    else {
        cfg = config_from_args(argc, argv);
    }

    Vec2D pos;
    Vec2D vel;
    Vec2D acc;
    srand(getpid() + time(NULL));
    init_system(&pos, &vel, &acc, cfg.bodies, cfg.dimensions);
    if (cfg.has_SDL) {
        simulate_with_graphics(&cfg, &pos, &vel, &acc);
    }
    else {
        simulate(&cfg, &pos, &vel, &acc);
    }
    deinit_system(&pos, &vel, &acc);

    return 0;
}
