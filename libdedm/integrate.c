#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <math.h>

#include "../allvars.h"

#include "integrate.h"
#include "dedmvars.h"
#include "interpolate.h"

#define WORKSIZE 100000
/*
 * integrate.c - contains methods and routines to integrate and store
 *		 functions introduced because of the dedm coupling
 * */


#ifdef DEDM_DRAG
void init_phidot_table()
{
	int i;
	double result, abserr;
	gsl_function F;
	gsl_integration_workspace *workspace;

	double logTimeBegin = log(All.TimeBegin);
	double logTimeMax = log(All.TimeMax);

	workspace = gsl_integration_workspace_alloc(WORKSIZE);

	for(i=0; i< PHIDOT_TABLE_LENGTH; i++){
		/* getPhiDot_Beta_a returns PhiDot/M_p * Beta (to allow a time varying coupling) */
	F.function = &get_a_PhiDot_Beta_a;

	gsl_integration_qag(&F, exp(logTimeBegin), 
			exp(logTimeBegin + ((logTimeMax - logTimeBegin) / PHIDOT_TABLE_LENGTH) * (i + 1)), 
			All.Hubble,  /* note: absolute error just a dummy */
			1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

	PhiDotIntegrationTable[i] = result;

	}

	gsl_integration_workspace_free(workspace);

}


/* returns the integral of the integral of phidot over M */
double get_phidot_factor(double time0, double time1)
{
	double a1, a2, df1, df2, u1, u2;
 	int i1, i2;

	double logTimeBegin = log(All.TimeBegin);
	double logTimeMax = log(All.TimeMax);

	a1 = logTimeBegin + time0 * All.Timebase_interval;
      	a2 = logTimeBegin + time1 * All.Timebase_interval;

        u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * PHIDOT_TABLE_LENGTH;
        i1 = (int) u1;
	    if(i1 >= PHIDOT_TABLE_LENGTH)
	    i1 = PHIDOT_TABLE_LENGTH - 1;

	      if(i1 <= 1)
		    df1 = u1 * PhiDotIntegrationTable[0];
	      else
 df1 = PhiDotIntegrationTable[i1 - 1] + (PhiDotIntegrationTable[i1] - PhiDotIntegrationTable[i1 - 1]) * (u1 - i1);
		    
	      u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * PHIDOT_TABLE_LENGTH;
	      i2 = (int) u2;

	        if(i2 >= PHIDOT_TABLE_LENGTH)
		    i2 = PHIDOT_TABLE_LENGTH - 1;

	      if(i2 <= 1)
	    df2 = u2 * PhiDotIntegrationTable[0];
	
	      else
 df2 = PhiDotIntegrationTable[i2 - 1] + (PhiDotIntegrationTable[i2] - PhiDotIntegrationTable[i2 - 1]) * (u2 - i2);

		    return df2 - df1;
}
#endif
