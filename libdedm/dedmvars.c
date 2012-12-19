/*
 * dedmvars.c
 * */
#include "dedmvars.h"
#include "../allvars.h"

struct table_interpolation Table;

#ifdef DEDM_MASS
struct mass_interpolation VariableMassTable;
#endif

#ifdef DEDM_HUBBLE
struct hubble_interpolation HubTable;
#endif

#ifdef DEDM_COUPLING
struct beta_interpolation BetaTable;
#endif

#if defined(DEDM_COUPLING) || defined(VARIABLE_G)
struct g_newton_interpolation GTable;

void setBetaZero()
{
	BetaTable.beta_0 = All.BetaZero;
#ifdef DEDM_INFO
	if( ThisTask == 0)
	{
		fprintf(stdout, "setG_tilde() Task=%d, BetaZero: %lf\n", ThisTask, BetaTable.beta_0);
		fprintf(All.outDeDmFile, "setG_tilde() Task=%d, BetaZero: %lf\n", ThisTask, BetaTable.beta_0);
	}
#endif
}

void setG_tilde() 
{
	GTable.G_tilde = 2*All.G*All.BetaZero*All.BetaZero;
#ifdef DEDM_INFO
	if( ThisTask == 0)
	{
		fprintf(stdout, "Task=%d. The old gravitational constant G: %lf\n is now set to: %lf for DM particles\n", 
			ThisTask, All.G, All.G+GTable.G_tilde);
		fprintf(All.outDeDmFile, 
			"Task=%d. The old gravitational constant G: %lf\n is now set to: %lf for DM particles\n",
			ThisTask, All.G, All.G+GTable.G_tilde);
	}
#endif
}
#endif

extern struct his History;
extern struct phi Phi;
extern struct phi_dot PhiDot;
