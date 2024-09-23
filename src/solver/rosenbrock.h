#ifndef ROSENBROCK_H
#define ROSENBROCK_H

// Function to open the Rosenbrock solver (allocate resources)
int rosenbrock_open(int n);

// Function to close the Rosenbrock solver (free resources)
void rosenbrock_close();

// Function to perform the integration
int rosenbrock_integrate(
        double ystart[], int n, double x1, double x2,
        double eps, double h1, void (*derivs)(double, double*, double*));

#endif // ROSENBROCK_H
