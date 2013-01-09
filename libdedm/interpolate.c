#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "../allvars.h"
#include "interpolate.h"
#include "dedmvars.h"

#ifdef DEDM_DRAG
#include "../libscott/read_scott_tables.h"
#endif

/*
 * interpolate.c : routine to read and interpolate values of:
 * 			- Mass
 * 			- Hubble constant
 * 			- Coupling parameter
 * 			- Phi_Dot
 *
 *     	   from a previously calculated table
 * */

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

double critical_density(double a)
{
	return 3*All.Hubble*All.Hubble*getH_a(a) / (8*All.G*PI);
}

#ifdef DEDM_HUBBLE
double getH_a(double a)
{
	double H;
	static int iter=0;

	H = get_interpolated_value(HubTable.a, HubTable.hubble, HubTable.npts ,a);
#ifdef DEDM_INFO
	if(All.Time > 0.02)
	{
		if(ThisTask==0 && iter < 5)
		{	
			fprintf(stdout, "Task=%d, Time=%lf, Hubble=%lf\n", ThisTask, All.Time, H);
		fprintf(All.outDeDmFile, "Task=%d, Time=%lf, Hubble=%lf\n", ThisTask, All.Time, H);
	}
	iter++;
	}
#endif
	return H;
}
#endif

#ifdef DEDM_COUPLING
double getB_a(double a, void *param)
{
	double B;
	/* return a constant coupling parameter */
	B = BetaTable.beta_0;

	/*TODO: Beta can be time varying*/
	// B = get_interpolated_value(BetaTable.a, BetaTable.beta, BetaTable.npts, a)
	return B;
}


double getG_a(double a)
{
	double G;
	/* G is modified by a constant amount  */
	G = GTable.G_tilde;

	/* TODO: G can be time-varying */
	// G = get_interpolated_value()
	return G;
}


#endif

#ifdef DEDM_MASS
double get_variable_mass_factor(double a)
{
	return get_interpolated_value(VariableMassTable.a, VariableMassTable.mass, VariableMassTable.npts ,a);
}
#endif


#ifdef DEDM_DRAG
double getPhiDot_a(double a, void *param){
	return sqrt(2./3.)/M_p_cmbeasy * history_phidot(a);
}
#endif

#if defined(DEDM_DRAG) && defined(DEDM_COUPLING)
	/* function used for GSL integration of a possible time dependent beta coupling */
double get_a_PhiDot_Beta_a(double a, void *p){
	return getPhiDot_a(a,p) * getB_a(a,p);
}
#endif

