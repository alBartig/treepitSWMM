#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rosenbrock.h"

#define MAX_ITER 100
#define SAFE 0.9
#define MIN_STEP 1e-12

// Structure to hold the solver's state
typedef struct {
    int nmax;
    double eps_global;
    double h_global;

    double *k1, *k2, *k3;
    double *y_temp;
    double **jacobian;
    double **a_matrix;
    int *ipiv;
} RosenbrockSolverState;

// Static variable to hold the solver's state
static RosenbrockSolverState *solver_state = NULL;

// Function to compute the Jacobian matrix numerically
static void compute_jacobian(
        int n,
        double t,
        double y[],
        void (*derivs)(double, double*, double*),
        RosenbrockSolverState *state)
{
    double *f0 = (double *)malloc(n * sizeof(double));
    double *y_temp = (double *)malloc(n * sizeof(double));
    derivs(t, y, f0);

    double epsilon = 1e-8;
    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i];
    }

    for (int j = 0; j < n; j++) {
        double temp = y[j];
        y_temp[j] = temp + epsilon;
        derivs(t, y_temp, state->jacobian[j]);
        for (int i = 0; i < n; i++) {
            state->jacobian[j][i] = (state->jacobian[j][i] - f0[i]) / epsilon;
        }
        y_temp[j] = temp;
    }

    free(f0);
    free(y_temp);
}

// Function to perform LU decomposition using Crout's method
static int lu_decompose(double **a, int n, int *ipiv)
{
    double *scales = (double *)malloc(n * sizeof(double));
    if (!scales) {
        return 0; // Memory allocation failed
    }

    for (int i = 0; i < n; i++) {
        double max = 0.0;
        for (int j = 0; j < n; j++) {
            double tmp = fabs(a[i][j]);
            if (tmp > max) {
                max = tmp;
            }
        }
        if (max == 0.0) {
            free(scales);
            return 0; // Singular matrix
        }
        scales[i] = 1.0 / max;
    }

    for (int k = 0; k < n; k++) {
        double max = 0.0;
        int imax = k;
        for (int i = k; i < n; i++) {
            double tmp = scales[i] * fabs(a[i][k]);
            if (tmp > max) {
                max = tmp;
                imax = i;
            }
        }

        if (imax != k) {
            // Swap rows
            double *temp_row = a[k];
            a[k] = a[imax];
            a[imax] = temp_row;
            int temp_piv = ipiv[k];
            ipiv[k] = ipiv[imax];
            ipiv[imax] = temp_piv;
            scales[imax] = scales[k];
        }

        if (a[k][k] == 0.0) {
            free(scales);
            return 0; // Singular matrix
        }

        for (int i = k + 1; i < n; i++) {
            a[i][k] /= a[k][k];
            for (int j = k + 1; j < n; j++) {
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }

    free(scales);
    return 1; // Success
}

// Function to solve linear system using LU decomposition
static void lu_solve(double **a, int n, int *ipiv, double *b, double *x)
{
    for (int i = 0; i < n; i++) {
        x[i] = b[i];
    }

    // Forward substitution
    for (int i = 0; i < n; i++) {
        int pivot = ipiv[i];
        double sum = x[pivot];
        x[pivot] = x[i];
        if (i > 0) {
            for (int j = 0; j < i; j++) {
                sum -= a[i][j] * x[j];
            }
        }
        x[i] = sum;
    }

    // Backward substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = x[i];
        if (i < n - 1) {
            for (int j = i + 1; j < n; j++) {
                sum -= a[i][j] * x[j];
            }
        }
        x[i] = sum / a[i][i];
    }
}

// Open the Rosenbrock solver and allocate resources
int rosenbrock_open(int n)
{
    solver_state = (RosenbrockSolverState *)malloc(sizeof(RosenbrockSolverState));
    if (!solver_state) {
        return 0; // Allocation failed
    }

    solver_state->nmax = n;
    solver_state->eps_global = 0.0;
    solver_state->h_global = 0.0;

    solver_state->k1 = (double *)calloc(n, sizeof(double));
    solver_state->k2 = (double *)calloc(n, sizeof(double));
    solver_state->k3 = (double *)calloc(n, sizeof(double));
    solver_state->y_temp = (double *)calloc(n, sizeof(double));

    solver_state->jacobian = (double **)calloc(n, sizeof(double *));
    solver_state->a_matrix = (double **)calloc(n, sizeof(double *));
    solver_state->ipiv = (int *)calloc(n, sizeof(int));

    if (!solver_state->k1 || !solver_state->k2 || !solver_state->k3 ||
        !solver_state->y_temp || !solver_state->jacobian ||
        !solver_state->a_matrix || !solver_state->ipiv) {
        rosenbrock_close();
        return 0; // Allocation failed
    }

    for (int i = 0; i < n; i++) {
        solver_state->jacobian[i] = (double *)calloc(n, sizeof(double));
        solver_state->a_matrix[i] = (double *)calloc(n, sizeof(double));
        if (!solver_state->jacobian[i] || !solver_state->a_matrix[i]) {
            rosenbrock_close();
            return 0; // Allocation failed
        }
    }

    return 1; // Success
}

// Integrate the ODEs using the Rosenbrock method
int rosenbrock_integrate(
        double ystart[], int n, double x1, double x2,
        double eps, double h1, void (*derivs)(double, double*, double*))
{
    if (!solver_state || n > solver_state->nmax) {
        printf("Error: Solver not initialized or n exceeds nmax\n");
        return 1; // Error code
    }

    int i;
    double t = x1;
    double h = h1;
    double t_end = x2;

    solver_state->eps_global = eps;
    solver_state->h_global = h1;

    // Copy initial values to y
    double *y = (double *)malloc(n * sizeof(double));
    if (!y) {
        printf("Memory allocation failed.\n");
        return 1; // Error code
    }
    for (i = 0; i < n; i++) {
        y[i] = ystart[i];
    }

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }

        // Compute Jacobian matrix at current time and state
        compute_jacobian(n, t, y, derivs, solver_state);

        // Compute function value at current time and state
        derivs(t, y, solver_state->k1);

        // Build the matrix A = I - gamma * h * J
        double gamma = 1.0 / (2.0 + sqrt(2.0));
        for (i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                solver_state->a_matrix[i][j] = -gamma * h * solver_state->jacobian[j][i];
            }
            solver_state->a_matrix[i][i] += 1.0; // Add identity matrix
        }

        // Perform LU decomposition of A
        for (i = 0; i < n; i++) {
            solver_state->ipiv[i] = i;
        }
        if (!lu_decompose(solver_state->a_matrix, n, solver_state->ipiv)) {
            printf("Error: Singular matrix in LU decomposition at t = %f\n", t);
            free(y);
            return 1; // Error code
        }

        // Compute k2
        for (i = 0; i < n; i++) {
            solver_state->k2[i] = solver_state->k1[i];
        }
        lu_solve(solver_state->a_matrix, n, solver_state->ipiv, solver_state->k2, solver_state->k2);

        // Update the solution
        for (i = 0; i < n; i++) {
            solver_state->y_temp[i] = y[i] + h * solver_state->k2[i];
        }

        // Estimate error
        double err = 0.0;
        for (i = 0; i < n; i++) {
            double y_err = h * (solver_state->k2[i] - solver_state->k1[i]);
            double denom = fabs(solver_state->y_temp[i]) + fabs(y[i]) + 1e-6;
            err += (y_err / denom) * (y_err / denom);
        }
        err = sqrt(err / n);

        // Adjust step size
        if (err > solver_state->eps_global) {
            // Reject step and reduce h
            h *= SAFE * pow(solver_state->eps_global / err, 0.25);
            if (h < MIN_STEP) {
                printf("Error: Step size too small at t = %f\n", t);
                free(y);
                return 1; // Error code
            }
            solver_state->h_global = h;
            continue; // Retry the step with new h
        } else {
            // Accept step and possibly increase h
            h *= SAFE * pow(solver_state->eps_global / err, 0.2);
            if (h > 5 * solver_state->h_global) {
                h = 5 * solver_state->h_global;
            }
            solver_state->h_global = h;
        }

        // Update time and solution
        t += h;
        for (i = 0; i < n; i++) {
            y[i] = solver_state->y_temp[i];
        }
    }

    // Copy final values back to ystart
    for (i = 0; i < n; i++) {
        ystart[i] = y[i];
    }

    free(y);
    return 0; // Success
}

// Close the Rosenbrock solver and free resources
void rosenbrock_close()
{
    if (!solver_state) {
        return;
    }

    if (solver_state->k1) free(solver_state->k1);
    if (solver_state->k2) free(solver_state->k2);
    if (solver_state->k3) free(solver_state->k3);
    if (solver_state->y_temp) free(solver_state->y_temp);
    if (solver_state->ipiv) free(solver_state->ipiv);

    if (solver_state->jacobian) {
        for (int i = 0; i < solver_state->nmax; i++) {
            if (solver_state->jacobian[i]) free(solver_state->jacobian[i]);
        }
        free(solver_state->jacobian);
    }

    if (solver_state->a_matrix) {
        for (int i = 0; i < solver_state->nmax; i++) {
            if (solver_state->a_matrix[i]) free(solver_state->a_matrix[i]);
        }
        free(solver_state->a_matrix);
    }

    free(solver_state);
    solver_state = NULL;
}
