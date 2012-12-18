#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "../allvars.h"
#include "interpolate.h"
#include "vdevars.h"
/*
 * interpolate.c : routine to read and interpolate values of:
 * 			- Hubble constant
 *     	   from a previously calculated table
 * */

/* Check correct initialization */
void print_hubble()
{
	int j;
	for(j=0; j<HubTable.npts ; j++) {
		fprintf(stderr, "%s %lf %s %lf \n", "a: ", HubTable.a[j], " h: ", HubTable.hubble[j]);
	}

}


double get_interpolated_value(double x_array[], double y_array[], int npts, double A)
{
	double V; 
 
		gsl_interp_accel *gia = gsl_interp_accel_alloc();	
		gsl_spline *gs = gsl_spline_alloc(gsl_interp_cspline, npts);
		gsl_spline_init(gs,x_array, y_array, npts);
		V = gsl_spline_eval(gs,A,gia);
		gsl_spline_free (gs);
		gsl_interp_accel_free (gia);

	return V;	
}


double getH_a(double a)
{
	return get_interpolated_value(HubTable.a, HubTable.hubble, HubTable.npts ,a);
}


double critical_density(double a)
{
	return 3*All.Hubble*All.Hubble*getH_a(a) / (8*All.G*PI);
}
