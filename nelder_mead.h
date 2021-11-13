/*
   Do this:
      #define NELDER_MEAD_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the implementation.

    
*/

#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

// MIT Licence.
// Copyright (c) 2017 Matteo Maggioni, 2021 Justin Meiners

#ifdef __cplusplus
extern "C"
{
#endif

//-----------------------------------------------------------------------------
// Definitions
//-----------------------------------------------------------------------------

// define optimization settings
typedef struct {
  double tolx;
  double tolf;
  int max_iter;
  int max_eval;
  int verbose;
} optimset_t;


// Cost function interface.
typedef double (*multivar_real_val_func_t)(int, const double*, const void *);

//-----------------------------------------------------------------------------
// Nelder-Mead simplex algorithm
// Optimization for multivariable real-valued functions having a real value, without derivatives.
//
// Parameters:
// - dimension of the data
// - initial point (unchanged in output)
// - solution_point is the minimizer
// - cost_function is a pointer to a real valued function to optimize 
// - args are the optional arguments passed to the cost_function
// - optimset are the optimisation settings
//
// Returns: the function value at the solution point
//-----------------------------------------------------------------------------
double nelder_mead(
    int dimension,
    const double *initial_point,
    double *solution_point,
    multivar_real_val_func_t cost_func,
    const void *args,
    const optimset_t *options
);


#ifdef __cplusplus
}
#endif

#endif // NELDER_MEAD_H

#ifdef NELDER_MEAD_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RHO 1.0
#define CHI 2.0
#define GAMMA 0.5
#define SIGMA 0.5

// define a generic point containing a position (x) and a value (fx)
typedef struct {
  double *x;
  double fx;
} point_t;

// define a simplex struct containing an array of n+1 points (p)
// each having dimension (n)
typedef struct {
  point_t *p;
  int n;
} simplex_t;

//-----------------------------------------------------------------------------
// Simplex sorting
//-----------------------------------------------------------------------------

static int compare(const void *arg1, const void *arg2) {
  double fx1 = ((const point_t *)arg1)->fx;
  double fx2 = ((const point_t *)arg2)->fx;
  return (fx1 > fx2) - (fx1 < fx2);
}

static void simplex_sort(simplex_t *simplex) {
  qsort(simplex->p, simplex->n + 1, sizeof(point_t), compare);
}

//-----------------------------------------------------------------------------
// Get centroid (average position) of simplex
//-----------------------------------------------------------------------------

static void get_centroid(const simplex_t *simplex, point_t *centroid) {
  double scale = 1.0 / (double)simplex->n;
  // (last point is not included)
  for (int j = 0; j < simplex->n; j++) {
    centroid->x[j] = 0;
    for (int i = 0; i < simplex->n; i++) {
      centroid->x[j] += simplex->p[i].x[j];
    }
    centroid->x[j] *= scale; 
  }
}

//-----------------------------------------------------------------------------
// Assess if simplex satisfies the minimization requirements
//-----------------------------------------------------------------------------

static int continue_minimization(const simplex_t *simplex, int eval_count,
                          int iter_count, const optimset_t *optimset) {
  if (eval_count > optimset->max_eval || iter_count > optimset->max_iter) {
    // stop if #evals or #iters are greater than the max allowed
    return 0;
  }
  int n = simplex->n;

  double condf = simplex->p[n].fx - simplex->p[0].fx;
  
  double condx = -1.0;
  for (int i = 1; i < simplex->n + 1; i++) {
    for (int j = 0; j < simplex->n; j++) {
      const double temp = fabs(simplex->p[0].x[j] - simplex->p[i].x[j]);
      if (condx < temp) {
        condx = temp;
      }
    }
  }
  // continue if both tolx or tolf condition is not met
  return condx > optimset->tolx || condf > optimset->tolf;
}

//-----------------------------------------------------------------------------
// Update current point
//-----------------------------------------------------------------------------

static void update_point(const simplex_t *simplex, const point_t *centroid,
                  double lambda, point_t *point) {
  const int n = simplex->n;
  for (int j = 0; j < n; j++) {
    point->x[j] = (1.0 + lambda) * centroid->x[j] - lambda * simplex->p[n].x[j];
  }
}

//-----------------------------------------------------------------------------
// Simple point_t manipulation utlities
//-----------------------------------------------------------------------------

static void copy_point(int n, const point_t *src, point_t *dst) {
  memcpy(dst->x, src->x, sizeof(double) * n);
  dst->fx = src->fx;
}

void swap_points(point_t *p1, point_t *p2) {
  point_t temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

double nelder_mead(
	int n,
	const double *start,
	double *solution,
	multivar_real_val_func_t cost_func,
	const void *args,
	const optimset_t *optimset
) {

  // internal points
  point_t point_r;
  point_t point_e;
  point_t point_c;
  point_t centroid;

  // allocate memory for internal points
  double* internal_buffer = malloc(n * 4 * sizeof(double));
  point_r.x = internal_buffer + 0 * n;
  point_e.x = internal_buffer + 1 * n; 
  point_c.x = internal_buffer + 2 * n;
  centroid.x = internal_buffer + 3 * n; 

  int iter_count = 0;
  int eval_count = 0;

  // initial simplex has size n + 1 where n is the dimensionality pf the data
  simplex_t simplex;
  simplex.n = n;
  simplex.p = malloc((n + 1) * sizeof(point_t));

  double* point_buffer = malloc((n + 1) * n * sizeof(double));

  for (int i = 0; i < n + 1; i++) {
    simplex.p[i].x = point_buffer + i * n; 
    for (int j = 0; j < n; j++) {
      simplex.p[i].x[j] =
          (i - 1 == j) ? (start[j] != 0.0 ? 1.05 * start[j] : 0.00025)
                       : start[j];
    }
    simplex.p[i].fx = cost_func(n, simplex.p[i].x, args);
    eval_count++;
  }
  // sort points in the simplex so that simplex.p[0] is the point having
  // minimum fx and simplex.p[n] is the one having the maximum fx
  simplex_sort(&simplex);
  // compute the simplex centroid
  get_centroid(&simplex, &centroid);
  iter_count++;

  // continue minimization until stop conditions are met
  while (continue_minimization(&simplex, eval_count, iter_count, optimset)) {
    int shrink = 0;

    if (optimset->verbose) {
      printf("Iteration %04d     ", iter_count);
    }
    update_point(&simplex, &centroid, RHO, &point_r);
    point_r.fx = cost_func(n, point_r.x, args);

    eval_count++;
    if (point_r.fx < simplex.p[0].fx) {
      update_point(&simplex, &centroid, RHO * CHI, &point_e);
      point_e.fx = cost_func(n, point_e.x, args);
      eval_count++;
      if (point_e.fx < point_r.fx) {
        // expand
        if (optimset->verbose) {
          printf("expand          ");
        }
        copy_point(n, &point_e, simplex.p + n);
      } else {
        // reflect
        if (optimset->verbose) {
          printf("reflect         ");
        }
        copy_point(n, &point_r, simplex.p + n);
      }
    } else {
      if (point_r.fx < simplex.p[n - 1].fx) {
        // reflect
        if (optimset->verbose) {
          printf("reflect         ");
        }
        copy_point(n, &point_r, simplex.p + n);
      } else {
        if (point_r.fx < simplex.p[n].fx) {
          update_point(&simplex, &centroid, RHO * GAMMA, &point_c);
          point_c.fx = cost_func(n, point_c.x, args);
          eval_count++;
          if (point_c.fx <= point_r.fx) {
            // contract outside
            if (optimset->verbose) {
              printf("contract out    ");
            }
            copy_point(n, &point_c, simplex.p + n);
          } else {
            // shrink
            if (optimset->verbose) {
              printf("shrink         ");
            }
            shrink = 1;
          }
        } else {
          update_point(&simplex, &centroid, -GAMMA, &point_c);
          point_c.fx = cost_func(n, point_c.x, args);
          eval_count++;
          if (point_c.fx <= simplex.p[n].fx) {
            // contract inside
            if (optimset->verbose) {
              printf("contract in     ");
            }
            copy_point(n, &point_c, simplex.p + n);
          } else {
            // shrink
            if (optimset->verbose) {
              printf("shrink         ");
            }
            shrink = 1;
          }
        }
      }
    }
    if (shrink) {
      for (int i = 1; i < n + 1; i++) {
        for (int j = 0; j < n; j++) {
          simplex.p[i].x[j] = simplex.p[0].x[j] +
                              SIGMA * (simplex.p[i].x[j] - simplex.p[0].x[j]);
        }

        simplex.p[i].fx = cost_func(n, simplex.p[i].x, args);
        eval_count++;
      }
      simplex_sort(&simplex);
    } else {
      for (int i = n - 1; i >= 0 && simplex.p[i + 1].fx < simplex.p[i].fx; i--) {
        swap_points(simplex.p + (i + 1), simplex.p + i);
      }
    }
    get_centroid(&simplex, &centroid);
    iter_count++;
    if (optimset->verbose) {
      // print current minimum
      printf("[ ");
      for (int i = 0; i < n; i++) {
        printf("%.2f ", simplex.p[0].x[i]);
      }
      printf("]    %.2f \n", simplex.p[0].fx);
    }
  }

  // save solution in output argument
  double min_fx = simplex.p[0].fx;
  if (solution) {
    memcpy(solution, simplex.p[0].x, n * sizeof(double));
  }

  // free memory
  free(internal_buffer);
  free(point_buffer);
  free(simplex.p);

  return min_fx;
}

#endif
