#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#ifdef __cplusplus
extern "C"
{
#endif

//-----------------------------------------------------------------------------
// Definitions
//-----------------------------------------------------------------------------

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

// define optimization settings
typedef struct {
  double tolx;
  double tolf;
  int max_iter;
  int max_eval;
  int verbose;
} optimset_t;

//-----------------------------------------------------------------------------
// Cost function interface
//-----------------------------------------------------------------------------

typedef double (*multivar_real_val_func_t)(int, const double*, const void *);

//-----------------------------------------------------------------------------
// Nelder-Mead algorithm
// - dimension of the data
// - initial point (unchanged in output)
// - solution_point is the minimizer
// - cost_function is a pointer to a real valued function to optimize 
// - args are the optional arguments passed to the cost_function
// - optimset are the optimisation settings
//-----------------------------------------------------------------------------
void nelder_mead(
    int dimension,
    const point_t *initial_point,
    point_t *solution_point,
    multivar_real_val_func_t cost_func,
    const void *args,
    const optimset_t *options
);

#ifdef __cplusplus
}
#endif

#endif // NELDER_MEAD_H
