/*
   Multivariable optimization without derivatives/gradients.
   Implementation of Nelder-Mead simplex algorithm.
   http://www.scholarpedia.org/article/Nelder-Mead_algorithm

Do this:

#define NELDER_MEAD_IMPLEMENTATION

before you include this file in *one* C or C++ file to create the implementation.
You can disable the debug_log feature to remove some extra checks and includes.

#define NM_NO_DEBUG_LOG

This version is written by Justin Meiners (2021), derived from Matteo Maggioni's (2017) work.
The full MIT License is listed at the end of the file.
 */

#ifndef NM_OPTIMIZER_H
#define NN_OPTIMIZER_H

#ifdef __cplusplus
extern "C"
{
#endif

// Types
//-----------------------------------------------------------------------------

#ifndef NM_REAL
#define NM_REAL double
#endif

// Cost function.

// Parameters:
// - number of variables
// - point (array of n values)
// - args is user  data 
typedef NM_REAL (*nm_multivar_real_func_t)(int, const NM_REAL *, void *);

// Parmeters:
// - tol_x: Terminate if any dimension of the simplex is smaller. 
// - tol_fx: Terminate if we see improvement less than this amount.
// - max_iterations: Terminate if we exceed this number of iterations.
// - restarts: How many time to try improving after a termination.
// - debug_log: whether to show debug text.

typedef struct {
    NM_REAL tol_x;
    NM_REAL tol_fx;
    int max_iterations;
    int restarts;
    int debug_log;
} nm_params_t;

typedef struct {
    int tol_satisfied;
    int iterations;
    NM_REAL min_fx;
} nm_result_t;

// CONVENIENCE API 
//-----------------------------------------------------------------------------

// These functions try to find the minimum using a variety of tricks and techniques on top of a simplex.
// They are supposed to be as much of a "black box" as possible.

// Parameters:
// - dimension: number of variables
// - initial: point to start search around (array of n values)
// - func: is a pointer to a cost function to optimize 
// - args: optional user arguments passed to the function.
// - params: tolerance and iteration control. 
// - out: the point which minimizes function (array of n values)
//
nm_result_t nm_multivar_optimize(
        int dimension,
        const NM_REAL *initial,
        const NM_REAL *initial_search_size,
        nm_multivar_real_func_t func,
        void *args,
        const nm_params_t *params,
        NM_REAL *out
        );


//-----------------------------------------------------------------------------
// DETAILED API
//-----------------------------------------------------------------------------

// If you want to write your own seach/restart strategy using simplex iteration
// this detailed API can help.


// An n-dimensional point x with it's function value fx. 
typedef struct {
    NM_REAL *x;
    NM_REAL fx;
} nm_simplex_pt_t;

// An n-dimensional simplex with n+1 simplex points.
typedef struct {
    int dimension;
    nm_simplex_pt_t *p;
    NM_REAL* p_buffer;
    NM_REAL* temp_buffer;
} nm_simplex_t;

void nm_simplex_init(nm_simplex_t* simplex, int dimension);
void nm_simplex_shutdown(nm_simplex_t* simplex);

// Guess the simplex size in each dimension using only the initial value.
void nm_guess_simplex_size(int n, const NM_REAL* initial, NM_REAL* out_size);

// Place the simplex in a standard position around the initial point.
// This is roughly offseting the intiial point by each vector in the standard basis.
/*
   |\
   | \
   |  \
   0---
 */

void nm_simplex_position_around(nm_simplex_t* simplex, const NM_REAL* initial, const NM_REAL* size);

int nm_simplex_iterate(
        nm_simplex_t* simplex,
        nm_multivar_real_func_t func,
        void *args,
        const nm_params_t *params
        );


#ifdef __cplusplus
}
#endif

#endif // NELDER_MEAD_H

#ifdef NM_OPTIMIZER_IMPLEMENTATION

// You can override these with #define if you really want to.
#ifndef NM_RHO
#define NM_RHO 1.0
#endif

#ifndef NM_CHI
#define NM_CHI 2.0
#endif 

#ifndef NM_GAMMA
#define NM_GAMMA 0.5
#endif

#ifndef NM_SIGMA
#define NM_SIGMA 0.5
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TEMP_POINT_COUNT_ 4

void nm_params_init_default(nm_params_t* params, int dimension) {
    params->tol_x = 0;
    params->tol_fx = 0.00001;
    params->max_iterations = 1000 + 1000 * (int)sqrt((double)dimension);

    params->restarts = 3;
    params->debug_log = 0;
}

void nm_simplex_init(nm_simplex_t* simplex, int dimension) {
    int n = dimension;
    simplex->p = malloc((n + 1) * sizeof(nm_simplex_pt_t));
    simplex->dimension = n;

    // Simplex has n + 1 points where n is the dimension pf the data.
    simplex->p_buffer = malloc((n + 1) * n * sizeof(NM_REAL));
    for (int i = 0; i < n + 1; i++) {
        simplex->p[i].x = simplex->p_buffer + i * n; 
    }

    simplex->temp_buffer = malloc(TEMP_POINT_COUNT_ * n * sizeof(NM_REAL));
}

void nm_simplex_shutdown(nm_simplex_t* simplex) {
    free(simplex->p_buffer);
    free(simplex->p);
    free(simplex->temp_buffer);
}

void nm_guess_simplex_size(int n, const NM_REAL* initial, NM_REAL* out_size) {
    for (int j = 0; j < n; ++j) {
        out_size[j] = (0.1 * fabs(initial[j])) + 0.005;
    }
}

void nm_simplex_position_around(nm_simplex_t* simplex, const NM_REAL* initial, const NM_REAL* size) {
    int n = simplex->dimension;
    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < n; ++j) {
            simplex->p[i].x[j] = initial[j] + ((i - 1 == j) ? size[j] : 0.0);
        }
    }
}

// Simplex sorting
static int cmp_pt_value_(const void *arg1, const void *arg2) {
    const nm_simplex_pt_t* p1 = arg1;
    const nm_simplex_pt_t* p2 = arg2;

    if (p1->fx < p2->fx) {
        return -1;
    } else if (p2->fx < p1->fx) {
        return 1;
    } else {
        return 0;
    }
}

static void simplex_sort(nm_simplex_t *simplex) {
    qsort(simplex->p, simplex->dimension + 1, sizeof(nm_simplex_pt_t), cmp_pt_value_);
}


// Get centroid (average position) of simplex
static void simplex_centroid(const nm_simplex_t *simplex, nm_simplex_pt_t *centroid) {
    int n = simplex->dimension;
    NM_REAL scale = 1.0 / (NM_REAL)n;
    // (last point is not included)
    for (int j = 0; j < n; j++) {
        centroid->x[j] = 0;
        for (int i = 0; i < n; i++) {
            centroid->x[j] += simplex->p[i].x[j];
        }
        centroid->x[j] *= scale; 
    }
}

static void update_images(nm_simplex_pt_t* start, nm_simplex_pt_t* end, int dimension, nm_multivar_real_func_t func, void* args) {
    while (start != end) {
        start->fx = func(dimension, start->x, args);
        ++start;
    }
}

static void simplex_shrink(nm_simplex_t* simplex) {
    int n = simplex->dimension;
    for (int i = 1; i < n + 1; i++) {
        for (int j = 0; j < n; j++) {
            simplex->p[i].x[j] = simplex->p[0].x[j] + NM_SIGMA * (simplex->p[i].x[j] - simplex->p[0].x[j]);
        }
    }
}

// Assess if simplex satisfies the minimization requirements

static int should_stop_(const nm_simplex_t *simplex, int iterations, const nm_params_t *params) {

    if (iterations >= params->max_iterations && params->max_iterations > 0) return 1;

    int n = simplex->dimension;

    if (params->tol_fx > 0.0) {
        NM_REAL best = fabs(simplex->p[0].fx);
        NM_REAL worst = fabs(simplex->p[n].fx);
        NM_REAL fractional_difference = 2.0 * (worst - best) / (worst + best);

        //printf("%f\n", fractional_difference);

        if (fractional_difference < params->tol_fx) return 1;
    }

    NM_REAL max_width = 0.0;
    if (params->tol_x > 0.0) {
        for (int i = 1; i < n + 1; i++) {
            for (int j = 0; j < n; j++) {
                const NM_REAL width = fabs(simplex->p[0].x[j] - simplex->p[i].x[j]);
                if (width > max_width) max_width = width;
            }
        }
        return max_width < params->tol_x; 
    } else {
        return 0;
    }
}

// Update current point
static void update_point(const nm_simplex_t *simplex, const nm_simplex_pt_t *centroid,
        NM_REAL lambda, nm_simplex_pt_t *point) {
    int n = simplex->dimension;
    for (int j = 0; j < n; j++) {
        point->x[j] = (1.0 + lambda) * centroid->x[j] - lambda * simplex->p[n].x[j];
    }
}

static void copy_point_(int n, const nm_simplex_pt_t *src, nm_simplex_pt_t *dst) {
    memcpy(dst->x, src->x, sizeof(NM_REAL) * n);
    dst->fx = src->fx;
}

static void swap_points_(nm_simplex_pt_t *p1, nm_simplex_pt_t *p2) {
    nm_simplex_pt_t temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}


int nm_simplex_iterate(
        nm_simplex_t* simplex,
        nm_multivar_real_func_t func,
        void *args,
        const nm_params_t *params
        ) {

#ifdef NM_NO_DEBUG_LOG
    int verbose = 0;
#else
    int verbose = params->debug_log;
#endif

    int n = simplex->dimension;
    // internal points
    nm_simplex_pt_t point_r;
    nm_simplex_pt_t point_e;
    nm_simplex_pt_t point_c;
    nm_simplex_pt_t centroid;

    point_r.x = simplex->temp_buffer + 0 * n;
    point_e.x = simplex->temp_buffer + 1 * n; 
    point_c.x = simplex->temp_buffer + 2 * n;
    centroid.x = simplex->temp_buffer + 3 * n; 

    // sort points in the simplex so that simplex.p[0] is the point having
    // minimum fx and simplex.p[n] is the one having the maximum fx
    update_images(simplex->p, simplex->p + (n + 1), n, func, args);
    simplex_sort(simplex);

    // compute the simplex centroid
    simplex_centroid(simplex, &centroid);

    int iterations = 0;

    // continue minimization until stop conditions are met
    while (!should_stop_(simplex, iterations, params)) {
        int shrink = 0;

        if (verbose) printf("Iteration %04d     ", iterations);

        update_point(simplex, &centroid, NM_RHO, &point_r);
        point_r.fx = func(n, point_r.x, args);

        if (point_r.fx < simplex->p[0].fx) {
            update_point(simplex, &centroid, NM_RHO * NM_CHI, &point_e);

            point_e.fx = func(n, point_e.x, args);
            if (point_e.fx < point_r.fx) {
                // expand
                if (verbose) printf("expand          ");
                copy_point_(n, &point_e, simplex->p + n);
            } else {
                // reflect
                if (verbose) printf("reflect         ");
                copy_point_(n, &point_r, simplex->p + n);
            }
        } else {
            if (point_r.fx < simplex->p[n - 1].fx) {
                // reflect
                if (verbose) printf("reflect         ");
                copy_point_(n, &point_r, simplex->p + n);
            } else {
                if (point_r.fx < simplex->p[n].fx) {
                    update_point(simplex, &centroid, NM_RHO * NM_GAMMA, &point_c);
                    point_c.fx = func(n, point_c.x, args);
                    if (point_c.fx <= point_r.fx) {
                        // contract outside
                        if (verbose) printf("contract out    ");
                        copy_point_(n, &point_c, simplex->p + n);
                    } else {
                        // shrink
                        if (verbose) printf("shrink         ");
                        shrink = 1;
                    }
                } else {
                    update_point(simplex, &centroid, -NM_GAMMA, &point_c);
                    point_c.fx = func(n, point_c.x, args);
                    if (point_c.fx <= simplex->p[n].fx) {
                        // contract inside
                        if (verbose) printf("contract in     ");
                        copy_point_(n, &point_c, simplex->p + n);
                    } else {
                        // shrink
                        if (verbose) printf("shrink         ");
                        shrink = 1;
                    }
                }
            }
        }
        if (shrink) {
            simplex_shrink(simplex);
            update_images(simplex->p + 1, simplex->p + (n + 1), n, func, args);
            simplex_sort(simplex);
        } else {
            for (int i = n - 1; i >= 0 && simplex->p[i + 1].fx < simplex->p[i].fx; i--) {
                swap_points_(simplex->p + (i + 1), simplex->p + i);
            }
        }
        simplex_centroid(simplex, &centroid);
        ++iterations;

        if (verbose) {
            // print current minimum
            printf("[ ");
            for (int i = 0; i < n; i++) {
                printf("%.2f ", simplex->p[0].x[i]);
            }
            printf("]    %.2f \n", simplex->p[0].fx);
        }
    }

    return iterations;
}

nm_result_t nm_multivar_optimize(
        int dimension,
        const NM_REAL *initial,
        const NM_REAL *initial_search_size,
        nm_multivar_real_func_t func,
        void *args,
        const nm_params_t *params,
        NM_REAL *out
        ) {

    nm_simplex_t simplex;
    nm_simplex_init(&simplex, dimension);

    nm_simplex_position_around(&simplex, initial, initial_search_size);

    nm_result_t r;
    r.iterations = nm_simplex_iterate(&simplex, func, args, params);
    r.tol_satisfied = r.iterations < params->max_iterations;

    for (int restarts = 0; restarts <= params->restarts; ++restarts) {
        nm_simplex_position_around(&simplex, simplex.p[0].x, initial_search_size);
        int iterations = nm_simplex_iterate(&simplex, func, args, params);
        r.iterations += iterations;
        if (!r.tol_satisfied) r.tol_satisfied = r.iterations < params->max_iterations;
    }

    r.min_fx = simplex.p[0].fx;
    // save solution in output argument
    if (out) memcpy(out, simplex.p[0].x, dimension * sizeof(NM_REAL));
    nm_simplex_shutdown(&simplex);
    return r;
}

#undef TEMP_POINT_COUNT

#endif

/*
   Copyright 2017 Matteo Maggioni, 2022 Justin Meiners

   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

