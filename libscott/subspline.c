#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

#include "subspline.h"
#include "../allvars.h"


SplineGroup * alloc_SplineGroup(const size_t n)
{
  SplineGroup * group = (SplineGroup*) malloc(sizeof(SplineGroup));
  group->length = n;
  group->x = (double*) malloc(sizeof(double)*group->length);
  return group;
}

SubSpline * alloc_SubSpline(const SplineGroup * group)
{
  SubSpline * spline = (SubSpline*) malloc(sizeof(SubSpline));
  spline->parent = group;
  spline->y = (double*) malloc(sizeof(double)*group->length);
  spline->accel=NULL;
  spline->interp=NULL;

  return spline;
}

double SubSplineInterp(SubSpline * spline,const double x)
{
  double value;
  int err = gsl_interp_eval_e(spline->interp,spline->parent->x,spline->y,x,spline->accel,&value);
  assert(err == GSL_SUCCESS);
  return value;
}

double SubSplineInterpDeriv(SubSpline * spline,const double x)
{
  double value;
  int err = gsl_interp_eval_deriv_e(spline->interp,spline->parent->x,spline->y,x,spline->accel,&value);
  assert(err == GSL_SUCCESS);
  return value;
}

double SubSplineInterpInteg(SubSpline * spline,const double a,const double b)
{
  double value;
  int err = gsl_interp_eval_integ_e(spline->interp,spline->parent->x,spline->y,a,b,spline->accel,&value);
  assert(err == GSL_SUCCESS);
  return value;
}

void arm_SubSpline(SubSpline * spline)
{
  /* Create a GSL interpolator */
  const gsl_interp_type * interp_type = gsl_interp_cspline;
  spline->interp = gsl_interp_alloc(interp_type,spline->parent->length);
  int err = gsl_interp_init(spline->interp,spline->parent->x,spline->y,spline->parent->length);
  assert(err == GSL_SUCCESS);
  spline->accel = gsl_interp_accel_alloc();
}

void free_SplineGroup(SplineGroup * group)
{
  free(group->x);
  free(group);
}

void free_SubSpline(SubSpline * spline)
{
  free(spline->y);
  if (spline->accel)  gsl_interp_accel_free(spline->accel);
  if (spline->interp) gsl_interp_free(spline->interp);
  free(spline);
}
