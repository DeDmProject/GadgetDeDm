#ifndef SPLINEIMPL_H
#define SPLINEIMPL_H 

struct Spline;
struct SplineImpl;
struct SplineImpl * CreateSplineImpl(struct Spline * spline);
void FreeSplineImpl(struct SplineImpl * impl);
double EvalSplineImpl(struct Spline * spline, double x);
double EvalSplineDerivImpl(struct Spline * spline, double x);

#endif /* SPLINEIMPL_H */
