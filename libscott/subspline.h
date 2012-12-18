#ifndef SUBSPLINE_H
#define SUBSPLINE_H

#ifdef __cplusplus
#extern "C" {
#endif

#include  <gsl/gsl_interp.h>

  struct SubSpline;
  typedef struct SubSpline SubSpline;
  
  struct SplineGroup {
    size_t length;
    double * x;
  };
  typedef struct SplineGroup SplineGroup;  
  
  struct SubSpline {
    double * y;
    const SplineGroup * parent;		/* Non-owning */
    gsl_interp * interp;
    gsl_interp_accel * accel;
  };

  SplineGroup * alloc_SplineGroup(const size_t n);
  SubSpline * alloc_SubSpline(const SplineGroup * group);

  double SubSplineInterp(SubSpline * spline,const double x);
  double SubSplineInterpDeriv(SubSpline * spline,const double x);
  double SubSplineInterpInteg(SubSpline * spline,const double a,const double b);

  void arm_SubSpline(SubSpline * spline);
  void free_SubSpline(SubSpline * spline);
  void free_SplineGroup(SplineGroup * group);  /* Only free the group after all children */

#ifdef __cplusplus
}
#endif

#endif
