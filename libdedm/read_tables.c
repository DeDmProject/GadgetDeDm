#include "integrate.h"
#include "interpolate.h"
#include "dedmvars.h"
#include "../allvars.h"

#ifdef DEDM_DRAG
#include "../libscott/read_scott_tables.h"
#endif
/*
 *	read_tables.c
 *	Reads input tables and stores them into strucures.
 * */


int get_lines(FILE *f)
{
	int n=-1;
	char line[2048];

	if(ThisTask==0)
	if(f==NULL) 
		fprintf(stderr,"\nError in get_lines(): input file could not be opened. \n");
	
		while(!feof(f))
		{
			fgets(line,2048,f);
			n++;
		}

		rewind(f);

	return n;
}


void read_interp_table(char * URL) 
{
	if(ThisTask==0)
		fprintf(stderr, "Reading from table %s\n", URL);

	int k=0;
	char dummyline[2048];
	FILE * file=NULL;
	
		file = fopen(URL,"r");
	        Table.npts = get_lines(file);
		Table.x = (double *) calloc(Table.npts, sizeof(double));
		Table.y = (double *) calloc(Table.npts, sizeof(double));

		for(k=0; k<Table.npts; k++)
		{
			fgets(dummyline,2048,file);
			sscanf(dummyline, "%lf %lf", &Table.x[k], &Table.y[k]);
		}

	fclose(file);
}


void read_all_interpolation_tables()
{
int k=0;
double hub_0, mass_0;
#ifdef DEDM_HUBBLE
	read_interp_table(All.HubbleDeDmFile); 
	// Force HubFactor = 1.00 at z=0
	hub_0 = get_interpolated_value(Table.x, Table.y, Table.npts, 1.0);

		HubTable.npts = Table.npts;
		HubTable.a = (double *) calloc(Table.npts, sizeof(double));
		HubTable.hubble = (double *) calloc(Table.npts, sizeof(double));

		for(k=0; k<Table.npts; k++) 
		{
			HubTable.a[k]=Table.x[k];
			HubTable.hubble[k]=Table.y[k]/hub_0; 
		}
#ifdef DEDM_INFO
		if(ThisTask==0)
			fprintf(stderr, "Task=%d, reading hub_0=%lf\n", ThisTask, hub_0);

		// All the tasks are being dumped to the out file for debugging pourposes
		fprintf(All.outDeDmFile, "Task=%d, reading hub_0=%lf\n", ThisTask, hub_0);
#endif
#endif


#ifdef DEDM_MASS
	read_interp_table(All.VariableMassDeDmFile); 

        VariableMassTable.npts = Table.npts;
	VariableMassTable.a = (double *) calloc(Table.npts, sizeof(double));
	VariableMassTable.mass = (double *) calloc(Table.npts, sizeof(double));

		// Ensure mass factor at z=0 is normalized to 1.00
		mass_0 = get_interpolated_value(Table.x, Table.y, Table.npts, 1.0);

		for(k=0; k<Table.npts; k++) 
		{
			VariableMassTable.a[k]=Table.x[k];
			VariableMassTable.mass[k]=Table.y[k]/mass_0;
		}
#ifdef DEDM_INFO
		if(ThisTask==0)
			fprintf(stderr, "Task=%d, reading mass_0=%lf\n", ThisTask, mass_0);

		fprintf(All.outDeDmFile, "Task=%d, reading mass_0=%lf\n", ThisTask, mass_0);
#endif
#endif
}
