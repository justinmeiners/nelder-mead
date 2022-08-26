#define NM_OPTIMIZER_IMPLEMENTATION
#define NM_NO_DEBUG_LOG

#include "nelder_mead.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define PI 3.1415926535897932384626433832795
#define SQUARE(x) ((x) * (x))

void print_point(double* x, double fx, int n) {
    printf("x = [ ");
    for (int i = 0; i < n; i++) printf("%.8f ", x[i]);
    printf("], fx = %.8f \n", fx);
}

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

double ackley_fun(int n, const double* x, void *arg) {
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

void test_ackley() {
  int n = 3;

  // cost function parameters
  ackley_param_t ackley_params; 
  ackley_params.a = 20.0;
  ackley_params.b = 0.2;
  ackley_params.c = 2.0 * PI;

  // optimisation settings
  nm_params_t params;
  nm_params_init_default(&params, n);
  params.debug_log = 1;

  double start[3] = { -2.1 -3.04, 4.5 };
  double range[3] = { 5.0, 5.0, 5.0 };
  double solution[3];

  nm_result_t result = nm_multivar_optimize(n, start, range, &ackley_fun, &ackley_params, &params, solution);

  printf("Initial point\n");
  double start_fx = ackley_fun(n, start, &ackley_params);
  print_point(start, start_fx, n);

  printf("Solution\n");
  print_point(solution, result.min_fx, n);

  if (result.tol_satisfied) printf("Tolerance critera met\n");
}

int main(int argc, const char *argv[]) {
    test_ackley();
    return 0;
}
