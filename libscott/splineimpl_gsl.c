#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "splineimpl.h"
#include "spline.h"


struct SplineImpl {
    gsl_interp * interp;
    gsl_interp_accel * accel;
};

struct SplineImpl * CreateSplineImpl(struct Spline * spline){
    struct SplineImpl * i = malloc(sizeof(*i));
    i->interp = gsl_interp_alloc(gsl_interp_cspline,spline->n);
    i->accel = gsl_interp_accel_alloc();
    gsl_interp_init(i->interp,spline->x,spline->y,spline->n);
    return i;
}

void FreeSplineImpl(struct SplineImpl * impl){
    if (impl){
        gsl_interp_accel_free(impl->accel);
        gsl_interp_free(impl->interp);
    }
    free(impl);
}

double EvalSplineImpl(struct Spline * spline, double x){
    gsl_set_error_handler_off();
    double result = NAN;
    int err = gsl_interp_eval_e(spline->impl->interp,
                                spline->x,
                                spline->y,
                                x,
                                spline->impl->accel,&result);

    if (err!=GSL_SUCCESS){
        fprintf(stderr,"Interpolation error at %e (Spline is [ %e : %e ])\n",
                x,
                spline->x[0],
                spline->x[spline->n-1]);
    }
    return result;
}
double EvalSplineDerivImpl(struct Spline * spline, double x){
    gsl_set_error_handler_off();
    double result = NAN;
    int err = gsl_interp_eval_deriv_e(spline->impl->interp,
				      spline->x,
				      spline->y,
				      x,
				      spline->impl->accel,&result);

    if (err!=GSL_SUCCESS){
	fprintf(stderr,"Interpolation error at %e (Spline is [ %e : %e ])\n",
		x,
		spline->x[0],
		spline->x[spline->n-1]);
    }
    return result;
}
