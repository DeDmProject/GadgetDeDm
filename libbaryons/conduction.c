#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef CONDUCTION


void conduction_smoothed_temperature(void)
{
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n;
  int ndone;
  long long ntot, ntotleft;
  int maxfill, source;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  MPI_Status status;


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
	}
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeDensity - NTask; i++)
	if(P[i].Type == 0)
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      for(j = 0; j < NTask; j++)
		Exportflag[j] = 0;

	      conduction_smoothed_evaluate(i, 0);

	      for(j = 0; j < NTask; j++)
		{
		  if(Exportflag[j])
		    {
		      DensDataIn[nexport].Pos[0] = P[i].Pos[0];
		      DensDataIn[nexport].Pos[1] = P[i].Pos[1];
		      DensDataIn[nexport].Pos[2] = P[i].Pos[2];

		      DensDataIn[nexport].Hsml = PPP[i].Hsml;

		      DensDataIn[nexport].Index = i;
		      DensDataIn[nexport].Task = j;
		      nexport++;
		      nsend_local[j]++;
		    }
		}
	    }


      qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeDensity)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_CONDUCT_A,
				   &DensDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				   MPI_BYTE, recvTask, TAG_CONDUCT_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }


	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    conduction_smoothed_evaluate(j, 1);


	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeDensity)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_CONDUCT_B,
				   &DensDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_CONDUCT_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = DensDataIn[source].Index;
#ifdef CONDUCTION
			  SphP[place].SmoothedEntr += DensDataPartialResult[source].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
			  SphP[place].GradEntr[0] += DensDataPartialResult[source].GradEntr[0];
			  SphP[place].GradEntr[1] += DensDataPartialResult[source].GradEntr[1];
			  SphP[place].GradEntr[2] += DensDataPartialResult[source].GradEntr[2];
#endif
#endif
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }


  /* do final operations on results */

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    SphP[i].SmoothedEntr /= pow(SphP[i].Density, GAMMA);
#ifdef CONDUCTION_SATURATION
	    SphP[i].GradEntr[0] /= pow(SphP[i].Density, GAMMA);
	    SphP[i].GradEntr[1] /= pow(SphP[i].Density, GAMMA);
	    SphP[i].GradEntr[2] /= pow(SphP[i].Density, GAMMA);
#endif
	  }
    }


  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void conduction_smoothed_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3, hinv4;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  FLOAT *pos;
  double smoothentr;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
#endif
#ifdef PERIODIC
  double boxSize, boxHalf;

#ifdef LONG_X
  double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
  double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
  double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif

  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif


  smoothentr = 0;

#ifdef CONDUCTION_SATURATION
  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      h = DensDataGet[target].Hsml;
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;



  startnode = All.MaxPart;

  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

	      smoothentr += mass_j * wk * pow(SphP[j].Density, GAMMA_MINUS1) * SphP[j].Entropy;

#ifdef CONDUCTION_SATURATION
	      if(r > 0)
		{
		  gradentr[0] += mass_j * dwk * dx / r * SphP[j].Entropy * pow(SphP[j].Density, GAMMA_MINUS1);
		  gradentr[1] += mass_j * dwk * dy / r * SphP[j].Entropy * pow(SphP[j].Density, GAMMA_MINUS1);
		  gradentr[2] += mass_j * dwk * dz / r * SphP[j].Entropy * pow(SphP[j].Density, GAMMA_MINUS1);
		}
#endif
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {
      SphP[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      SphP[target].GradEntr[0] = gradentr[0];
      SphP[target].GradEntr[1] = gradentr[1];
      SphP[target].GradEntr[2] = gradentr[2];
#endif
    }
  else
    {
      DensDataResult[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      DensDataResult[target].GradEntr[0] = gradentr[0];
      DensDataResult[target].GradEntr[1] = gradentr[1];
      DensDataResult[target].GradEntr[2] = gradentr[2];
#endif
    }
}


#endif
