#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SFR_METALS
#include "libbaryons/c_metals.h"
#endif

#ifdef COSMIC_RAYS
#include "libbaryons/cosmic_rays.h"
#endif

#ifdef VDE
#include "libvde/interpolate.h"
#endif

#ifdef DEDM_HUBBLE
#include "libdedm/interpolate.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */

#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  long long ntot, ntotleft;
  int ndone;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, NNN;
  int npleft;
  int maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  double dt_entr, dmax1, dmax2, fac;
  double tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  double timengb, sumtimengb;
  double desnumngb;

#ifdef NAVIERSTOKES
  int k;
  double dvel[3][3];
  double rotx, roty, rotz;
#endif
#ifdef SFR_PROMOTION
  double xhyd, yhel, ne, mu, energy, temp;
#endif

  MPI_Status status;

#if defined(SOFTEREQS) || defined(SFR_PROMOTION)
  double a3inv;
#endif
#ifdef REDUCEVISC
  double cs_h, f, hubble_a = 1.0;
#endif

#ifdef PERIODIC
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

#ifdef REDUCEVISC
  if(All.ComovingIntegrationOn)
    {
#if defined(DEDM_HUBBLE) || defined(VDE)
	    hubble_a=getH_a(All.Time);
#else
      	    /* Factors for comoving integration of hydro */
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif
#endif
    }
#endif

#if defined(SOFTEREQS) || defined(SFR_PROMOTION)
  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

#if defined(SFR_METALS) || defined(BLACK_HOLES)
  NNN = NumPart;
#else
  NNN = N_gas;
#endif

  for(n = 0, NumSphUpdate = 0; n < NNN; n++)
    {
#ifdef BLACK_HOLES
      if(P[n].Type == 0 || P[n].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
      if(P[n].Type == 0 || P[n].Type == 6 || P[n].Type == 7)
#else
      if(P[n].Type == 0)
#endif
#endif
#endif
	{
	  PPP[n].Left = PPP[n].Right = 0;

	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
	}
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


  desnumngb = All.DesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < NNN && nexport < All.BunchSizeDensity - NTask; i++)
#ifdef BLACK_HOLES
	    if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
	    if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
	    if(P[i].Type == 0)
#endif
#endif
#endif
	      if(P[i].Ti_endstep == All.Ti_Current)
		{
		  ndone++;

		  for(j = 0; j < NTask; j++)
		    Exportflag[j] = 0;

		  density_evaluate(i, 0);

		  for(j = 0; j < NTask; j++)
		    {
		      if(Exportflag[j])
			{
			  DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			  DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			  DensDataIn[nexport].Pos[2] = P[i].Pos[2];

#if defined(SFR_METALS) || defined(BLACK_HOLES)
			  if(P[i].Type == 0)
#endif
			    {
#ifdef SFR_DECOUPLING
			      DensDataIn[nexport].DensityOld = SphP[i].DensityOld;
			      DensDataIn[nexport].Entropy = SphP[i].Entropy;
#endif
			      DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			      DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			      DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
#ifdef WINDS
			      DensDataIn[nexport].DelayTime = SphP[i].DelayTime;
#endif
			    }
#ifdef SFR_DECOUPLING
			  else
			    {
			      DensDataIn[nexport].DensityOld = 0;
			    }
#endif

			  DensDataIn[nexport].Hsml = PPP[i].Hsml;
			  DensDataIn[nexport].Index = i;
			  DensDataIn[nexport].Task = j;
			  nexport++;
			  nsend_local[j]++;
			}
		    }
		}
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

	  for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
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
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		{
		  density_evaluate(j, 1);
		}
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
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
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;

			      PPP[place].NumNgb += DensDataPartialResult[source].Ngb;

#if defined(SFR_METALS) || defined(BLACK_HOLES)
			      if(P[place].Type == 0)
#endif
				{
				  SphP[place].Density += DensDataPartialResult[source].Rho;
#ifndef NOGRADHSML
				  SphP[place].DhsmlDensityFactor +=
				    DensDataPartialResult[source].DhsmlDensity;
#endif
#ifndef NAVIERSTOKES
				  SphP[place].u.s.DivVel += DensDataPartialResult[source].Div;
				  SphP[place].u.s.Rot[0] += DensDataPartialResult[source].Rot[0];
				  SphP[place].u.s.Rot[1] += DensDataPartialResult[source].Rot[1];
				  SphP[place].u.s.Rot[2] += DensDataPartialResult[source].Rot[2];
#else
				  for(k = 0; k < 3; k++)
				    {
				      SphP[place].u.DV[k][0] += DensDataPartialResult[source].DV[k][0];
				      SphP[place].u.DV[k][1] += DensDataPartialResult[source].DV[k][1];
				      SphP[place].u.DV[k][2] += DensDataPartialResult[source].DV[k][2];
				    }
#endif
#ifdef CONDUCTION
				  SphP[place].SmoothedEntr += DensDataPartialResult[source].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
				  SphP[place].GradEntr[0] += DensDataPartialResult[source].GradEntr[0];
				  SphP[place].GradEntr[1] += DensDataPartialResult[source].GradEntr[1];
				  SphP[place].GradEntr[2] += DensDataPartialResult[source].GradEntr[2];
#endif
#endif
				}
#ifdef BLACK_HOLES
			      if(P[place].Type == 5)
				{
				  P[place].BH_Density += DensDataPartialResult[source].Rho;
				  P[place].BH_Entropy += DensDataPartialResult[source].SmoothedEntr;
				  P[place].BH_SurroundingGasVel[0] += DensDataPartialResult[source].GasVel[0];
				  P[place].BH_SurroundingGasVel[1] += DensDataPartialResult[source].GasVel[1];
				  P[place].BH_SurroundingGasVel[2] += DensDataPartialResult[source].GasVel[2];
				}
#endif

			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}



      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < NNN; i++)
	{
#ifdef BLACK_HOLES
	  if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
	  if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
	  if(P[i].Type == 0)
#endif
#endif
#endif
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
#if defined(SFR_METALS) || defined(BLACK_HOLES)
		if(P[i].Type == 0)
#endif
		  {

		    if(SphP[i].Density > 0)
		      {
#ifndef NOGRADHSML
                        SphP[i].DhsmlDensityFactor *= PPP[i].Hsml / (NUMDIMS * SphP[i].Density);
                        if(SphP[i].DhsmlDensityFactor > -0.9) /* note: this would be -1 if only a single particle at zero lag is found */ 
                          SphP[i].DhsmlDensityFactor = 1 / (1 + SphP[i].DhsmlDensityFactor);
                        else
                          SphP[i].DhsmlDensityFactor = 1;

                          /*
			SphP[i].DhsmlDensityFactor =
			  1 / (1 + PPP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));
                          */
#endif

#ifndef NAVIERSTOKES
			SphP[i].u.s.CurlVel = sqrt(SphP[i].u.s.Rot[0] * SphP[i].u.s.Rot[0] +
						   SphP[i].u.s.Rot[1] * SphP[i].u.s.Rot[1] +
						   SphP[i].u.s.Rot[2] * SphP[i].u.s.Rot[2]) / SphP[i].Density;

			SphP[i].u.s.DivVel /= SphP[i].Density;
#else
			for(k = 0; k < 3; k++)
			  {
			    dvel[k][0] = SphP[i].u.DV[k][0] / SphP[i].Density;
			    dvel[k][1] = SphP[i].u.DV[k][1] / SphP[i].Density;
			    dvel[k][2] = SphP[i].u.DV[k][2] / SphP[i].Density;
			  }
			SphP[i].u.s.DivVel = dvel[0][0] + dvel[1][1] + dvel[2][2];

			SphP[i].u.s.StressDiag[0] = 2 * dvel[0][0] - 2.0 / 3 * SphP[i].u.s.DivVel;
			SphP[i].u.s.StressDiag[1] = 2 * dvel[1][1] - 2.0 / 3 * SphP[i].u.s.DivVel;
			SphP[i].u.s.StressDiag[2] = 2 * dvel[2][2] - 2.0 / 3 * SphP[i].u.s.DivVel;

			SphP[i].u.s.StressOffDiag[0] = dvel[0][1] + dvel[1][0];	/* xy */
			SphP[i].u.s.StressOffDiag[1] = dvel[0][2] + dvel[2][0];	/* xz */
			SphP[i].u.s.StressOffDiag[2] = dvel[1][2] + dvel[2][1];	/* yz */

			rotx = dvel[1][2] - dvel[2][1];
			roty = dvel[2][0] - dvel[0][2];
			rotz = dvel[0][1] - dvel[1][0];
			SphP[i].u.s.CurlVel = sqrt(rotx * rotx + roty * roty + rotz * rotz);
#endif


#ifdef CONDUCTION
			SphP[i].SmoothedEntr /= SphP[i].Density;
#ifdef CONDUCTION_SATURATION
			SphP[i].GradEntr[0] /= SphP[i].Density;
			SphP[i].GradEntr[1] /= SphP[i].Density;
			SphP[i].GradEntr[2] /= SphP[i].Density;
#endif
#endif
		      }

		    dt_entr =
		      (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

#ifndef MHM
#ifndef SOFTEREQS
		    SphP[i].Pressure =
		      (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#else
		    /* use an intermediate EQS, between isothermal and the full multiphase model */
		    if(SphP[i].Density * a3inv >= All.PhysDensThresh)
		      SphP[i].Pressure = All.FactorForSofterEQS *
			(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA) +
			(1 - All.FactorForSofterEQS) * GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;
		    else
		      SphP[i].Pressure =
			(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#endif
#else
		    /* Here we use an isothermal equation of state */
		    SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Density * All.InitGasU;
		    SphP[i].Entropy = SphP[i].Pressure / pow(SphP[i].Density, GAMMA);
#endif

#ifdef COSMIC_RAYS
		    CR_Particle_Update(SphP + i);

		    SphP[i].Pressure += CR_Particle_Pressure(SphP + i);
#endif


#ifdef REDUCEVISC
		    if(SphP[i].Density > 0.0)
		      {
			/*
			   cs_h = DMAX(sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density),
			   10 * (sqrt(SphP[i].VelPred[0] * SphP[i].VelPred[0] + 
			   SphP[i].VelPred[1] * SphP[i].VelPred[1] + 
			   SphP[i].VelPred[2] * SphP[i].VelPred[2]) +
			   SphP[i].u.s.DivVel * SphP[i].Hsml)) /
			   SphP[i].Hsml;
			 */
			cs_h = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density) / SphP[i].Hsml;
			f =
			  fabs(SphP[i].u.s.DivVel) / (fabs(SphP[i].u.s.DivVel) + SphP[i].u.s.CurlVel +
						      0.0001 * cs_h);
			if(All.ComovingIntegrationOn)
			  f *= sqrt(All.Time);
			SphP[i].Dtalpha = -(SphP[i].alpha - All.AlphaMin) * All.DecayTime * cs_h
			  + f * All.ViscSource * DMAX(0.0, -SphP[i].u.s.DivVel);
			if(All.ComovingIntegrationOn)
			  SphP[i].Dtalpha /= (hubble_a * All.Time);
		      }
		    else
		      SphP[i].Dtalpha = 0.0;
#endif

		  }

#ifdef BLACK_HOLES
		if(P[i].Type == 5)
		  {
		    if(P[i].BH_Density > 0)
		      {
			P[i].BH_Entropy /= P[i].BH_Density;
			P[i].BH_SurroundingGasVel[0] /= P[i].BH_Density;
			P[i].BH_SurroundingGasVel[1] /= P[i].BH_Density;
			P[i].BH_SurroundingGasVel[2] /= P[i].BH_Density;
		      }
		  }
#endif

		/* now check whether we had enough neighbours */

#ifdef BLACK_HOLES
		if(P[i].Type == 5)
		  desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
		else
		  desnumngb = All.DesNumNgb;
#endif

		if(PPP[i].NumNgb < (desnumngb - All.MaxNumNgbDeviation) ||
		   (PPP[i].NumNgb > (desnumngb + All.MaxNumNgbDeviation)
		    && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
		  {
		    /* need to redo this particle */
		    npleft++;

		    if(PPP[i].Left > 0 && PPP[i].Right > 0)
		      if((PPP[i].Right - PPP[i].Left) < 1.0e-3 * PPP[i].Left)
			{
			  /* this one should be ok */
			  npleft--;
			  P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			  continue;
			}

		    if(PPP[i].NumNgb < (desnumngb - All.MaxNumNgbDeviation))
		      PPP[i].Left = DMAX(PPP[i].Hsml, PPP[i].Left);
		    else
		      {
			if(PPP[i].Right != 0)
			  {
			    if(PPP[i].Hsml < PPP[i].Right)
			      PPP[i].Right = PPP[i].Hsml;
			  }
			else
			  PPP[i].Right = PPP[i].Hsml;
		      }

		    if(iter >= MAXITER - 10)
		      {
			printf
			  ("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			   i, ThisTask, (int) P[i].ID, PPP[i].Hsml, PPP[i].Left, PPP[i].Right,
			   (float) PPP[i].NumNgb, PPP[i].Right - PPP[i].Left, P[i].Pos[0], P[i].Pos[1],
			   P[i].Pos[2]);
			fflush(stdout);
		      }

		    if(PPP[i].Right > 0 && PPP[i].Left > 0)
		      PPP[i].Hsml = pow(0.5 * (pow(PPP[i].Left, 3) + pow(PPP[i].Right, 3)), 1.0 / 3);
		    else
		      {
			if(PPP[i].Right == 0 && PPP[i].Left == 0)
			  endrun(8188);	/* can't occur */

			if(PPP[i].Right == 0 && PPP[i].Left > 0)
			  {
#ifndef NOGRADHSML
			    if(P[i].Type == 0 && fabs(PPP[i].NumNgb - desnumngb) < 0.5 * desnumngb)
			      {
				fac = 1 - (PPP[i].NumNgb -
					   desnumngb) / (NUMDIMS * PPP[i].NumNgb) *
				  SphP[i].DhsmlDensityFactor;

				if(fac < 1.26)
				  PPP[i].Hsml *= fac;
				else
				  PPP[i].Hsml *= 1.26;
			      }
			    else
			      PPP[i].Hsml *= 1.26;
#else
			    PPP[i].Hsml *= 1.26;
#endif
			  }

			if(PPP[i].Right > 0 && PPP[i].Left == 0)
			  {
#ifndef NOGRADHSML
			    if(P[i].Type == 0 && fabs(PPP[i].NumNgb - desnumngb) < 0.5 * desnumngb)
			      {
				fac = 1 - (PPP[i].NumNgb -
					   desnumngb) / (NUMDIMS * PPP[i].NumNgb) *
				  SphP[i].DhsmlDensityFactor;

				if(fac > 1 / 1.26)
				  PPP[i].Hsml *= fac;
				else
				  PPP[i].Hsml /= 1.26;
			      }
			    else
			      PPP[i].Hsml /= 1.26;
#else
			    PPP[i].Hsml /= 1.26;
#endif
			  }
		      }

		    if(PPP[i].Hsml < All.MinGasHsml)
		      PPP[i].Hsml = All.MinGasHsml;
		  }
		else
		  P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	      }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  for(i = 0; i < NNN; i++)
    {
#if defined(SFR_METALS) && defined(SFR_PROMOTION)
      if(P[i].Type == 0 && P[i].EnergySN < 0)	/* only gas particles should enter here */
	{
	  P[i].EnergySN = 0;
	  SphP[i].Entropy *= pow(SphP[i].DensityAvg / SphP[i].Density, GAMMA_MINUS1);

	  printf("Entropy Conservation: Part=%d factor=%g EnergySN=%g\n", P[i].ID,
		 pow(SphP[i].DensityAvg / SphP[i].Density, GAMMA_MINUS1), P[i].EnergySN);
	  fflush(stdout);
	}

      if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
	{
	  xhyd = P[i].Zm[6] / P[i].Mass;
	  yhel = (1 - xhyd) / (4. * xhyd);
	  ne = SphP[i].Ne;
	  mu = (1 + 4 * yhel) / (1 + yhel + ne);
	  energy = SphP[i].Entropy * P[i].Mass / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
	  temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[i].Mass * PROTONMASS * mu;
	  temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

	  fprintf(FdPromotion, "%g %d %g %g %g %g %g %g\n", All.Time, P[i].ID, SphP[i].DensPromotion,
		  SphP[i].TempPromotion, SphP[i].Density, temp, SphP[i].DensityAvg, SphP[i].EntropyAvg);
	  fflush(FdPromotion);


	  SphP[i].TempPromotion = 0;
	  SphP[i].DensPromotion = 0;
	}
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR_METALS
      if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
      if(P[i].Type == 0)
#endif
#endif
	{
	  /* mark as active again */
	  if(P[i].Ti_endstep < 0)
	    P[i].Ti_endstep = -P[i].Ti_endstep - 1;
	}
    }
  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */

  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz;

#ifndef NAVIERSTOKES
  double divv, rotv[3];
#else
  int k;
  double dvel[3][3];
#endif

  FLOAT *pos, *vel;

#if defined(SFR_METALS) || defined(BLACK_HOLES)
  static FLOAT veldummy[3] = { 0, 0, 0 };
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
  double smoothentr;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
#endif
#endif
#ifdef WINDS
  double delaytime;
#endif
#ifndef NOFIXEDMASSINKERNEL
  double weighted_numngb;
#endif
#ifndef NOGRADHSML
  double dhsmlrho;
#endif
#ifdef SFR_DECOUPLING
  double densityold, entropy;
#endif
#ifdef BLACK_HOLES
  double gasvel[3];
#endif

  rho = 0;

#ifndef NAVIERSTOKES
  divv = rotv[0] = rotv[1] = rotv[2] = 0;
#else
  for(k = 0; k < 3; k++)
    dvel[k][0] = dvel[k][1] = dvel[k][2] = 0;
#endif


#ifndef NOFIXEDMASSINKERNEL
  weighted_numngb = 0;
#endif
#ifndef NOGRADHSML
  dhsmlrho = 0;
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
  smoothentr = 0;
#ifdef CONDUCTION_SATURATION
  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif
#endif

#ifdef BLACK_HOLES
  gasvel[0] = gasvel[1] = gasvel[2] = 0;
#endif


  if(mode == 0)
    {
      pos = P[target].Pos;
#if defined(SFR_METALS) || defined(BLACK_HOLES)
      vel = veldummy;
      if(P[target].Type == 0)
#endif
	{
	  vel = SphP[target].VelPred;
#ifdef WINDS
	  delaytime = SphP[target].DelayTime;
#endif
	}
      h = PPP[target].Hsml;
#ifdef SFR_DECOUPLING
      if(P[target].Type == 0)
	{
	  densityold = SphP[target].DensityOld;
	  entropy = SphP[target].Entropy;
	}
      else
	densityold = 0;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
#ifdef WINDS
      delaytime = DensDataGet[target].DelayTime;
#endif
#ifdef SFR_DECOUPLING
      densityold = DensDataGet[target].DensityOld;
      entropy = DensDataGet[target].Entropy;
#endif
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
  numngb = 0;
  do
    {
#ifdef SFR_DECOUPLING
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, densityold, entropy, &vel[0]);
#else
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);
#endif

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];
#ifdef WINDS
	  if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
	    if(!(delaytime > 0))	/* if I'm not wind, then ignore the wind particle */
	      continue;
#endif
#ifdef BLACK_HOLES
	  if(P[j].Mass == 0)
	    continue;
#endif
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
	      numngb++;

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

	      rho += mass_j * wk;

#if defined(CONDUCTION) || defined(BLACK_HOLES)
	      smoothentr += mass_j * wk * SphP[j].Entropy;
#ifdef CONDUCTION_SATURATION
	      if(r > 0)
		{
		  gradentr[0] += mass_j * dwk * dx / r * SphP[j].Entropy;
		  gradentr[1] += mass_j * dwk * dy / r * SphP[j].Entropy;
		  gradentr[2] += mass_j * dwk * dz / r * SphP[j].Entropy;
		}
#endif
#endif

#ifdef BLACK_HOLES
	      gasvel[0] += mass_j * wk * SphP[j].VelPred[0];
	      gasvel[1] += mass_j * wk * SphP[j].VelPred[1];
	      gasvel[2] += mass_j * wk * SphP[j].VelPred[2];
#endif


#ifndef NOFIXEDMASSINKERNEL
	      weighted_numngb += NORM_COEFF * wk / hinv3;	/* 4.0/3 * PI = 4.188790204786 */
#endif

#ifndef NOGRADHSML
	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);
#endif
	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
#ifndef NAVIERSTOKES
		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);
#else
		  dvel[0][0] -= fac * dx * dvx;
		  dvel[0][1] -= fac * dx * dvy;
		  dvel[0][2] -= fac * dx * dvz;
		  dvel[1][0] -= fac * dy * dvx;
		  dvel[1][1] -= fac * dy * dvy;
		  dvel[1][2] -= fac * dy * dvz;
		  dvel[2][0] -= fac * dz * dvx;
		  dvel[2][1] -= fac * dz * dvy;
		  dvel[2][2] -= fac * dz * dvz;
#endif
		}
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {
#ifndef NOFIXEDMASSINKERNEL
      PPP[target].NumNgb = weighted_numngb;
#else
      PPP[target].NumNgb = numngb;
#endif

#if defined(SFR_METALS) || defined(BLACK_HOLES)
      if(P[target].Type == 0)
#endif
	{
	  SphP[target].Density = rho;
#ifdef CONDUCTION
	  SphP[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
	  SphP[target].GradEntr[0] = gradentr[0];
	  SphP[target].GradEntr[1] = gradentr[1];
	  SphP[target].GradEntr[2] = gradentr[2];
#endif
#endif


#ifndef NOGRADHSML
	  SphP[target].DhsmlDensityFactor = dhsmlrho;
#endif

#ifndef NAVIERSTOKES
	  SphP[target].u.s.DivVel = divv;
	  SphP[target].u.s.Rot[0] = rotv[0];
	  SphP[target].u.s.Rot[1] = rotv[1];
	  SphP[target].u.s.Rot[2] = rotv[2];
#else
	  for(k = 0; k < 3; k++)
	    {
	      SphP[target].u.DV[k][0] = dvel[k][0];
	      SphP[target].u.DV[k][1] = dvel[k][1];
	      SphP[target].u.DV[k][2] = dvel[k][2];
	    }
#endif
	}

#ifdef BLACK_HOLES
      P[target].BH_Density = rho;
      P[target].BH_Entropy = smoothentr;
      P[target].BH_SurroundingGasVel[0] = gasvel[0];
      P[target].BH_SurroundingGasVel[1] = gasvel[1];
      P[target].BH_SurroundingGasVel[2] = gasvel[2];
#endif
    }
  else
    {
      DensDataResult[target].Rho = rho;
#ifndef NAVIERSTOKES
      DensDataResult[target].Div = divv;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
#else
      for(k = 0; k < 3; k++)
	{
	  DensDataResult[target].DV[k][0] = dvel[k][0];
	  DensDataResult[target].DV[k][1] = dvel[k][1];
	  DensDataResult[target].DV[k][2] = dvel[k][2];
	}
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
      DensDataResult[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      DensDataResult[target].GradEntr[0] = gradentr[0];
      DensDataResult[target].GradEntr[1] = gradentr[1];
      DensDataResult[target].GradEntr[2] = gradentr[2];
#endif
#endif

#ifndef NOFIXEDMASSINKERNEL
      DensDataResult[target].Ngb = weighted_numngb;
#else
      DensDataResult[target].Ngb = numngb;
#endif

#ifndef NOGRADHSML
      DensDataResult[target].DhsmlDensity = dhsmlrho;
#endif

#ifdef BLACK_HOLES
      DensDataResult[target].GasVel[0] = gasvel[0];
      DensDataResult[target].GasVel[1] = gasvel[1];
      DensDataResult[target].GasVel[2] = gasvel[2];
#endif
    }
}





/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}


#ifdef NAVIERSTOKES
double get_shear_viscosity(int i)
{
  return All.NavierStokes_KinematicViscosity * SphP[i].Density;
}
#endif
