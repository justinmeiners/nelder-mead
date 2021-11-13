/*
   MIT Licence.
   Copyright (c) 2017 Matteo Maggioni, 2021 Justin Meiners

   Do this:
      #define NELDER_MEAD_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the implementation.


*/

#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#ifdef __cplusplus
extern "C"
{
#endif


//-----------------------------------------------------------------------------
// Definitions
//-----------------------------------------------------------------------------

// Cost function interface.
// Multivariable real-valued function.

// Parameters:
// - dimension of the data
// - input point (array of n values)
// - args is user provided data 
typedef double (*nm_multivar_real_func_t)(int, const double *, void *);

typedef struct {
  double tolx;
  double tolf;
  int max_iter;
  int verbose;
} nm_optimset_t;

//-----------------------------------------------------------------------------
// Nelder-Mead simplex algorithm
// Multivariable optimization without derivatives.

// See: http://www.scholarpedia.org/article/Nelder-Mead_algorithm
//
// Parameters:
// - dimension of the data
// - start is the initial point to search around (array of n values)
// - out is the point which minimizes function (array of n values)
// - out_val is the cost_func
// - cost_fuc is a pointer to a function to optimize 
// - args are the optional arguments passed to the cost_function
// - optimset is tolerance and iteration control. 
//
// Returns:
// - 1: found solution within tolerances.
// - 0: reached iteration limit.
//-----------------------------------------------------------------------------
int nm_multivar_optimize(
    int dimension,
    const double *start,
    double *out,
    double *out_val,
    nm_multivar_real_func_t cost_func,
    void *args,
    const nm_optimset_t *optimset
);

// Same as above but with additional options.

// Parameters
// - simplex_scale is how large the initial simplex is (array of n-values)

int nm_multivar_optimize_opt(
    int dimension,
    const double *start,
    double *out,
    double *out_val,
    nm_multivar_real_func_t cost_func,
    void *args,
    const nm_optimset_t *optimset,
    const double* simplex_scale
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

static int continue_minimization(const simplex_t *simplex, 
                          int iter_count, const nm_optimset_t *optimset) {
  if (iter_count > optimset->max_iter) {
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

static void swap_points(point_t *p1, point_t *p2) {
  point_t temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

int nm_multivar_optimize(
	int n,
	const double *start,
	double *out,
    double *out_val,
	nm_multivar_real_func_t cost_func,
	void *args,
	const nm_optimset_t *optimset
) {
    double* simplex_scale = malloc(n * sizeof(double));

    int result = nm_multivar_optimize_opt(n, start, out, out_val, cost_func, args, optimset, simplex_scale);

    free(simplex_scale);
    return result;
}



#define RHO 1.0
#define CHI 2.0
#define GAMMA 0.5
#define SIGMA 0.5

int nm_multivar_optimize_opt(
    int n,
    const double *start,
    double *out,
    double *out_val,
    nm_multivar_real_func_t cost_func,
    void *args,
    const nm_optimset_t *optimset,
    const double* simplex_scale
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
  }
  // sort points in the simplex so that simplex.p[0] is the point having
  // minimum fx and simplex.p[n] is the one having the maximum fx
  simplex_sort(&simplex);
  // compute the simplex centroid
  get_centroid(&simplex, &centroid);
  iter_count++;

  // continue minimization until stop conditions are met
  while (continue_minimization(&simplex, iter_count, optimset)) {
    int shrink = 0;

    if (optimset->verbose) {
      printf("Iteration %04d     ", iter_count);
    }
    update_point(&simplex, &centroid, RHO, &point_r);
    point_r.fx = cost_func(n, point_r.x, args);

    if (point_r.fx < simplex.p[0].fx) {
      update_point(&simplex, &centroid, RHO * CHI, &point_e);
      point_e.fx = cost_func(n, point_e.x, args);
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
  if (out) {
    memcpy(out, simplex.p[0].x, n * sizeof(double));
  }
  if (out_val) {
     *out_val = simplex.p[0].fx;
  }

  // free memory
  free(internal_buffer);
  free(point_buffer);
  free(simplex.p);

  return iter_count < optimset->max_iter;
}

#undef RHO
#undef CHI
#undef GAMMA
#undef SIGMA


#endif
