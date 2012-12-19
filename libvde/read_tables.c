#include "vdevars.h"
#include "../allvars.h"


void generate_vde_hubble_table()
{
	/* Generate test tables */
	int k;
	FILE *output_vde;
	double b=0, z=0;

	output_vde = fopen("vde_hubble_vs_a.dat", "w");	
	
		for (k=0; k< HubTable.npts; k++) 
		{
			b = HubTable.hubble[k];
			z = 1./HubTable.a[k]-1;  
			fprintf(output_vde, "%g %f \n", z, b);
		}

	fclose(output_vde);
}


void invertA(int dim)
{
	int i=0;
	double *v1;

		v1 = (double *) calloc(dim, sizeof(double));

		for(i=0; i< dim; i++)
		{
			v1[i] = HubTable.a[dim-i];
		}

	for(i=0; i< dim; i++)
		HubTable.a[i] = v1[i];
  free(v1);
}


void invertH(int dim)
{
	int i=0;
	double *v1;

	v1 = (double *) calloc(dim, sizeof(double));

		for(i=0; i < dim; i++)
		{
			v1[i] = HubTable.hubble[dim-i];
		}

		for(i=0; i < dim; i++)
		{
			HubTable.hubble[i] = v1[i];
		}

	free(v1);
}


void printVectors(double vec1[], double vec2[], int dim)
{
	int i=0;
	for(i=0; i<dim; i++) 
	{
		fprintf(stdout, "%d %s %lf %s %lf \n", i, " vector1: ", vec1[i], " vector2: ", vec2[i]);
	}
}


int get_lines(FILE *f)
{
	int n=-1;
	char line[2048];

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
	int npts=0, k=0;
	double hh, zz;
	char dummyline[2048];
	FILE * file=NULL;
	
	file = fopen(URL,"r");
	npts = get_lines(file);
	
		if(file==NULL)
		{
			fprintf(stderr,"%s %s \n", "\nError in read_interp_table(). Could not open file: ", URL );
			exit(0);
		}

		if(ThisTask==0) 
			fprintf(stdout, " Found VDE hubble table h(a) file, allocating memory ... \n");

		HubTable.npts   = npts;
		HubTable.a      = (double *) calloc(npts, sizeof(double));
		HubTable.hubble = (double *) calloc(npts, sizeof(double));

			for(k=0; k<npts; k++)
			{
				fgets(dummyline,2048,file);

					// VDE Table contains only z and H 
		        	sscanf(dummyline, "%lf %lf", &zz, &hh);
				HubTable.a[npts - k - 1] = 1./(1.+zz);
			
				//NOTE: we need to tabulate H^2 in order 
				//      to be consistent with the calls to getH_a()
			HubTable.hubble[npts - k - 1] = hh*hh; 
		}

	fclose(file);
}


void read_all_interpolation_tables()
{
	char *urls = All.HubbleVDEFile;
	read_interp_table(urls);
}
