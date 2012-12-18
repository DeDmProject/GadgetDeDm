#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#define M_PI		3.14159265358979323846	/* pi */

#include "spline.h"
#include "subspline.h"
#include "read_scott_tables.h"
#include "../allvars.h"

Spline * history_phidot_spline;


Spline * SingleRead(const char * key){
	FILE * file = fopen(key,"r");
	Spline * s = ReadSpline(file);
	fclose(file);
	return s;
}


void load_phidot_spline()
{
  /* H_spline is assumed to be uninitialised when this is called */
 	if(ThisTask==0)
 	   printf("Loading history splines\n");

  {
      Spline * his = SingleRead(All.HisDeDmFile);
      Spline * quint = SingleRead(All.PhiDeDmFile);
      Spline * phi = CreateSpline(his->y,quint->y,his->n);
       FreeSpline(his);
       FreeSpline(quint);

      history_phidot_spline = CopySpline(phi);
	int i = 0; 
      for (i=0;i<phi->n;++i){
	history_phidot_spline->y[i] = EvalSplineDeriv(phi,phi->x[i]);
      }

      RearmSpline(history_phidot_spline);
      FreeSpline(phi);
  }

}


void free_phidot_spline()
{
    FreeSpline(history_phidot_spline);
}


double history_phidot(double a){
    return EvalSpline(history_phidot_spline,a);
}
