#ifndef SPLINE_H
#define SPLINE_H 

#include <stddef.h>
#include <stdio.h>

struct SplineImpl;
typedef struct Spline {
    size_t n;
    double * x;
    double * y;
    struct SplineImpl * impl;
} Spline;

Spline * CreateSpline(double * x,
                      double * y,
                      size_t n);
void FreeSpline(Spline * spline);

Spline * CopySpline(Spline * spline);
Spline * ReadSpline(FILE * stream);

// Returns interpolation at x
double EvalSpline(Spline * spline,
                  double x);
double EvalSplineDeriv(Spline * spline,
                  double x);

// Call after changing the vectors
void RearmSpline(Spline * spline);

#endif /* SPLINE_H */
