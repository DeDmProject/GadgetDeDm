#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef DEDM_TREE
#include "libdedm/dedmvars.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif


/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

/*! This function computes the gravitational forces for all active particles.
 *  If needed, a new tree is constructed, otherwise the dynamically updated
 *  tree is used.  Particles are only exported to other processors when really
 *  needed, thereby allowing a good use of the communication buffer.
 */
void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

  int k;

#ifdef DEDM_TREE
  double ax_dm, ay_dm, az_dm;
  double fac_dm;

  if(ThisTask==0)
 	 fprintf(stdout, "Using modified tree algorithm to compute short range dark energy-dark matter interaction. \n");
#endif

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif
#ifdef STATICNFW
  double r, m;
#endif
#ifdef STATICHQ
  double r, m, a;
#endif
#ifdef EVALPOTENTIAL
  double r2;
#endif


  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

#if defined(SFR) || defined(BLACK_HOLES)
      rearrange_particle_sequence();
#endif
      force_treebuild();

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;


#if defined(BLACK_HOLES)
  /* in this case, the number of active particles may have changed 
   * since NumForceUpdate was computed in run.c due to elimination 
   * of particles (black hole mergers) in rearrange_particle_sequence() 
   * Need to recompute NumForceUpdate for this reason.
   */
  for(i = 0, NumForceUpdate = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      NumForceUpdate++;
#endif

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  if(ntot == 0)
    return;

#ifndef NOGRAVITY
#ifndef ISOTHERM


  /* Note: 'NumForceUpdate' has been determined in find_next_sync_point_and_drift() */

  if(ThisTask == 0)
    printf("Begin tree force.\n");

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */


  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
#endif
	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++) {
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];	
#ifdef DEDM_TREE
		      if(P[i].Type == 1){
	GravDataGet[nexport].u_dm.Pos_dm[k] = P[i].Pos[k];
		      }
#endif		      
		    }
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].v.Type = P[i].Type;
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;

#ifdef DEDM_TREE 
		    if(P[i].Type == 1) { // TODO Add other particle types
	GravDataGet[nexport].v.Type = P[i].Type;
	//GravDataGet[nexport].w_dm.OldAcc_dm = P[i].OldAcc_dm;
		//     } else {
 	// do something? 
		     }
#endif
		    

		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;

		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
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
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++) {
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

#ifdef DEDM_TREE 
			P[place].DeDmAccel[k] += GravDataOut[j + noffset[recvTask]].u_dm.Acc_dm[k];
#endif
			}
			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
#ifdef EVALPOTENTIAL
			  P[place].Potential += GravDataOut[j + noffset[recvTask]].v.Potential;
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

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++) 
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    
    }
#endif
#endif

  for(i = 0; i < NumPart; i++) {
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;

#if defined(DEDM_TREE) || defined(DEDM_PM) || defined (DEDM_PMb)	
if(P[i].Type == 1){ //TODO add other particle types

#ifdef DEDM_TREE /* add dedm interaction components from the tree algorithm part*/
	ax_dm = P[i].DeDmAccel[0];
	ay_dm = P[i].DeDmAccel[1];
	az_dm = P[i].DeDmAccel[2];	

#if defined(DEDM_PM) || defined(DEDM_PMb)
//	if(GTable.G_tilde != 0){ // avoid division by zero
	ax_dm += P[i].DeDmPM[0] / GTable.G_tilde;
	ay_dm += P[i].DeDmPM[1] / GTable.G_tilde;
	az_dm += P[i].DeDmPM[2] / GTable.G_tilde;
//	}
#endif // dedm2
#endif // dedm5
	/*ax += ax_dm;
	ay += ay_dm;
	az += az_dm;*/
} else {
#ifdef DEDM_TREE
	ax_dm = 0;
	ay_dm = 0; 
	az_dm = 0;
#endif // dedm5
}
#endif // dedm
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif 
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);

      }
}

if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)  {
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	for(j = 0; j < 3; j++)
	{  
		P[i].GravAccel[j] *= All.G;

#ifdef DEDM_TREE
	if(P[i].Type == 1)
	{ 
		P[i].DeDmAccel[j] *= GTable.G_tilde; 
	}
#endif
	}

#ifdef EVALPOTENTIAL
	/* remove self-potential */
	P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];

	if(All.ComovingIntegrationOn)
	  if(All.PeriodicBoundariesOn)
	    P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	      pow(All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G), 1.0 / 3);

	P[i].Potential *= All.G;

#ifdef PMGRID
	P[i].Potential += P[i].PM_Potential;	/* add in long-range potential */
#endif

	if(All.ComovingIntegrationOn)
	  {
#ifndef PERIODIC
	    fac = -0.5 * All.Omega0 * All.Hubble * All.Hubble;

	    for(i = 0; i < NumPart; i++)
	      {
		for(k = 0, r2 = 0; k < 3; k++)
		  r2 += P[i].Pos[k] * P[i].Pos[k];

		P[i].Potential += fac * r2;
	      }
#endif
	  }
	else
	  {
	    fac = -0.5 * All.OmegaLambda * All.Hubble * All.Hubble;
	    if(fac != 0)
	      {
		for(i = 0; i < NumPart; i++)
		  {
		    for(k = 0, r2 = 0; k < 3; k++)
		      r2 += P[i].Pos[k] * P[i].Pos[k];

		    P[i].Potential += fac * r2;
		  }
	      }
	  }
#endif // endif EVALPOTENTIAL
      }}
  

  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;
#ifdef DARKENERGY
      fac *= -0.5 * (1 + 3 * DarkEnergy_t(All.Time));
#endif
      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* beginning of ISOTHERM-stuff */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	fac = 1 / (P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
	for(j = 0; j < 3; j++)
	  P[i].GravAccel[j] = -ISOTHERM * ISOTHERM * P[i].Pos[j] * fac;
      }
#endif


#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif





#ifdef STATICNFW
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
	m = enclosed_mass(r);
#ifdef NFW_DARKFRACTION
	m *= NFW_DARKFRACTION;
#endif
	if(r > 0)
	  {
	    for(k = 0; k < 3; k++)
	      P[i].GravAccel[k] += -All.G * m * P[i].Pos[k] / (r * r * r);
	  }
      }
#endif



#ifdef STATICHQ
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);

	a = pow(All.G * HQ_M200 / (100 * All.Hubble * All.Hubble), 1.0 / 3) / HQ_C *
	  sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));

	m = HQ_M200 * pow(r / (r + a), 2);

#ifdef HQ_DARKFRACTION
	m *= HQ_DARKFRACTION;
#endif
	if(r > 0)
	  {
	    for(k = 0; k < 3; k++)
	      P[i].GravAccel[k] += -All.G * m * P[i].Pos[k] / (r * r * r);
	  }
      }
#endif



  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}

#ifndef NOGRAVITY
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);
#endif

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*!  This function sets the (comoving) softening length of all particle types
 *  in the table All.SofteningTable[...].  We check that the physical softening
 *  length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
    All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
  else
    All.SofteningTable[0] = All.SofteningGas;

  if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
    All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
  else
    All.SofteningTable[1] = All.SofteningHalo;

  if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
    All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
  else
    All.SofteningTable[2] = All.SofteningDisk;

  if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
    All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
  else
    All.SofteningTable[3] = All.SofteningBulge;

  if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
    All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
  else
    All.SofteningTable[4] = All.SofteningStars;

  if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
    All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
  else
    All.SofteningTable[5] = All.SofteningBndry;

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}

/*! This function is used as a comparison kernel in a sort routine, which is
 *  the method used to group particles that are going to the same CPU in the
 *  communication buffer.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
