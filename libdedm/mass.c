#include "dedmvars.h"
#include "interpolate.h"
#include "../allvars.h"

/*
 * mass.c - routines and methods for variable mass particles
 * */

#ifdef DEDM_MASS
void initialize_m0()
{
	int i=0;
	for(i=0; i<6; i++)
	{
		if(ThisTask==0) 
			fprintf(stdout, "Initializing m0() for DM to:%lf \n", All.MassTable[i]);
		}

	VariableMassTable.mass_0=All.MassTable[1];
}


void modify_particles_masses() 
{
	int n=0;
	double delta_m;
	delta_m = get_variable_mass_factor(All.Time);
#ifdef DEDM_INFO
	if(ThisTask==0) 
		fprintf(stdout, "Task=%d, Calculating new particle's masses, modified by a factor of=%lf.\n", 
				ThisTask, delta_m);
		fprintf(All.outDeDmFile, "Task=%d, Calculating new particle's masses, modifieda by a factor of=%lf.\n", 
				ThisTask, delta_m);
#endif
     	for(n=0; n<NumPart; n++)
			if(P[n].Type==1) 
				P[n].Mass = All.MassTable[1]*delta_m;
}                  
#endif
