#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nelder_mead.h"

#define PI 3.1415926535897932384626433832795
#define SQUARE(x) ((x) * (x))

//-----------------------------------------------------------------------------
// Implementation of a cost function f : R^n->R compatible with fun_t
// In this instance we use the Ackley Function as it allows us to demonstrate
// the use of the optional fixed arguments.
// More details on the function at http://www.sfu.ca/%7Essurjano/ackley.html
//-----------------------------------------------------------------------------

typedef struct {
  double a;
  double b;
  double c;
} ackley_param_t;

double ackley_fun(int n, const double* x, const void *arg) {
  // cast the void pointer to what we expect to find
  const ackley_param_t *params = (const ackley_param_t *)arg;

  // cost function computation for arguments of exp
  double sum_squares = 0;
  double sum_cos = 0;
  for (int i = 0; i < n; i++) {
    sum_squares += SQUARE(x[i]);
    sum_cos += cos(params->c * x[i]);
  }

  // final result
  return -params->a * exp(-params->b * sqrt(sum_squares / n)) -
              exp(sum_cos / n) + params->a + exp(1.0);
}

//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------

int main(int argc, const char *argv[]) {
  if (argc == 1) {
    printf("%s: error: not enough inputs \n", argv[0]);
    return 0;
  }

  // reading initial point from command line
  const int n = argc - 1;
  point_t start; // initial point
  start.x = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    start.x[i] = atof(argv[i + 1]);
  }

  // optimisation settings
  optimset_t optimset;
  optimset.tolx = 0.001;    // tolerance on the simplex solutions coordinates
  optimset.tolf = 0.001;    // tolerance on the function value
  optimset.max_iter = 1000; // maximum number of allowed iterations
  optimset.max_eval = 1000; // maximum number of allowed function evaluations
  optimset.verbose = 0;     // toggle verbose output during minimization

  // cost function parameters
  ackley_param_t ackley_params; // parameters of our cost function
  ackley_params.a = 20.0;
  ackley_params.b = 0.2;
  ackley_params.c = 2.0 * PI;

  // call optimization methods
  point_t solution;    // container for the solution of the minimisation
  nelder_mead(n, &start, &solution, &ackley_fun, &ackley_params, &optimset);

  // evaluate and print starting point
  printf("Initial point\n");
  double start_fx = ackley_fun(n, start.x, &ackley_params);
  printf("x = [ ");
  for (int i = 0; i < n; i++) {
    printf("%.8f ", start.x[i]);
  }
  printf("], fx = %.8f \n", start_fx);
  // print solution
  printf("Solution\n");
  printf("x = [ ");
  for (int i = 0; i < n; i++) {
    printf("%.8f ", solution.x[i]);
  }
  printf("], fx = %.8f \n", solution.fx);

  // free memory
  free(start.x);
  free(solution.x);

  return 0;
}
