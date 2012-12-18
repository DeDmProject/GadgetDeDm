#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"

#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

#ifdef VDE
#include "libvde/interpolate.h"
#endif

#ifdef DEDM_HUBBLE
#include "libdedm/interpolate.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */



static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;

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

#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
static double u1;
#endif

/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void hydro_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask;
  int nexport, place;
  double soundspeed_i;
  double tstart, tend;
  double sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;

#ifdef NAVIERSTOKES
  double fac;
#endif
#ifdef CONDUCTION
  double rEntropyFlow, dt, fac;
#endif

#if defined(CR_SHOCK)
  double rShockEnergy;
  double rNonRethermalizedEnergy;
#endif

#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
  /* mean particle placing */
  u1 = pow(4. * M_PI / 3 / All.DesNumNgb, 1. / 3.);
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

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
#if defined(DEDM_HUBBLE) || defined(VDE)
	hubble_a = getH_a(All.Time);
#else	
      	    hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif // ifdef Darkenergy
#endif // ifdef DEDM

      hubble_a = All.Hubble * sqrt(hubble_a);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;


  /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
#ifdef CONDUCTION
      SphP[n].CondEnergyChange = 0.0;
#endif

#ifdef SFR
      if(P[n].Type == 0)
#endif
	if(P[n].Ti_endstep == All.Ti_Current)
	  NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);


  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
#ifdef SFR
	if(P[i].Type == 0)
#endif
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      for(j = 0; j < NTask; j++)
		Exportflag[j] = 0;

	      hydro_evaluate(i, 0);

	      for(j = 0; j < NTask; j++)
		{
		  if(Exportflag[j])
		    {
		      for(k = 0; k < 3; k++)
			{
			  HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			  HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
			}
		      HydroDataIn[nexport].Hsml = PPP[i].Hsml;
		      HydroDataIn[nexport].Mass = P[i].Mass;
#ifndef NOGRADHSML
		      HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
#endif
		      HydroDataIn[nexport].Density = SphP[i].Density;
		      HydroDataIn[nexport].Pressure = SphP[i].Pressure;
		      HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;

		      /* calculation of F1 */
		      soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
#ifndef ALTVISCOSITY
		      HydroDataIn[nexport].F1 = fabs(SphP[i].u.s.DivVel) /
			(fabs(SphP[i].u.s.DivVel) + SphP[i].u.s.CurlVel +
			 0.0001 * soundspeed_i / PPP[i].Hsml / fac_mu);
#else
		      HydroDataIn[nexport].F1 = SphP[i].u.s.DivVel;
#endif

#ifdef MAGNETIC
		      for(k = 0; k < 3; k++)
			HydroDataIn[nexport].BPred[k] = SphP[i].BPred[k];
#endif

#ifdef CONDUCTION
		      HydroDataIn[nexport].SmoothedEntr = SphP[i].SmoothedEntr;
		      HydroDataIn[nexport].Entropy = SphP[i].Entropy;
#ifdef CONDUCTION_SATURATION
		      HydroDataIn[nexport].GradEntr[0] = SphP[i].GradEntr[0];
		      HydroDataIn[nexport].GradEntr[1] = SphP[i].GradEntr[1];
		      HydroDataIn[nexport].GradEntr[2] = SphP[i].GradEntr[2];
#endif
#endif

#ifdef REDUCEVISC
		      HydroDataIn[nexport].alpha = SphP[i].alpha;
#endif

#ifdef SFR_DECOUPLING
/*		      HydroDataIn[nexport].DensityOld = SphP[i].DensityOld; */
		      HydroDataIn[nexport].DensityOld = SphP[i].Density;
		      HydroDataIn[nexport].Entropy = SphP[i].Entropy;
#endif


#ifdef PARTICLE_DEBUG
		      HydroDataIn[nexport].ID = P[i].ID;
#endif

#ifdef CR_DIFFUSION
		      HydroDataIn[nexport].CR_Diff_E = SphP[i].CR_Diff_E;
		      HydroDataIn[nexport].CR_Diff_N = SphP[i].CR_Diff_N;
#endif

#ifdef NAVIERSTOKES
		      for(k = 0; k < 3; k++)
			{
			  HydroDataIn[nexport].stressdiag[k] = SphP[i].u.s.StressDiag[k];
			  HydroDataIn[nexport].stressoffdiag[k] = SphP[i].u.s.StressOffDiag[k];
			}
		      HydroDataIn[nexport].shear_viscosity = get_shear_viscosity(i);
#endif


		      HydroDataIn[nexport].Index = i;
		      HydroDataIn[nexport].Task = j;
		      nexport++;
		      nsend_local[j]++;
		    }
		}
	    }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

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
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  /* now do the imported particles */
	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
	      hydro_evaluate(j, 1);
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
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = HydroDataIn[source].Index;

			  for(k = 0; k < 3; k++)
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];

			  SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;

			  if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
			    {
			      SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
			    }

#ifdef CONDUCTION
			  SphP[place].CondEnergyChange += HydroDataPartialResult[source].CondEnergyChange;
#ifdef OUTPUTCOOLRATE
			  SphP[place].CondRate += HydroDataPartialResult[source].CondRate;
#endif
#endif
#ifdef MAGNETIC
			  for(k = 0; k < 3; k++)
			    SphP[place].DtB[k] += HydroDataPartialResult[source].DtB[k];
#ifdef TRACEDIVB
			  SphP[place].divB += HydroDataPartialResult[source].divB;
#endif
#endif
#if defined ( CR_DIFFUSION )
			  CR_Particle_Inject(SphP + place,
					     HydroDataPartialResult[source].CR_EnergyChange,
					     HydroDataPartialResult[source].CR_BaryonFractionChange);
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



  /* do final operations on results */
  tstart = second();

  for(i = 0; i < N_gas; i++)
#ifdef SFR
    if(P[i].Type == 0)
#endif
      if(P[i].Ti_endstep == All.Ti_Current)
	{

#ifdef CR_SHOCK
	  /* state right here:
	   *
	   * _c denotes comoving quantities
	   * _p denotes physical quantities
	   *
	   *
	   * Delta u_p = rho_p^(gamma-1)/(gamma-1) Delta A
	   *
	   * Delta A = dA/dloga * Delta loga
	   *
	   * dA/dloga = DtE * (gamma-1) / ( H(a) a^2 rho_c^(gamma-1)
	   *
	   * => Delta u_p = DtE * dloga / ( H(a) a^(3gamma-1) )
	   */

	  if(SphP[i].DtEntropy > 0.0)
	    {
	      rShockEnergy = SphP[i].DtEntropy *
		(P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval / hubble_a2 / fac_egy;

	      /* Feed fraction "All.CR_ShockEfficiency" into CR and see what
	       * amount of energy instantly gets rethermalized
	       *
	       * for this, we need the physical time step, which is
	       * Delta t_p = Delta t_c / hubble_a
	       */

	      rNonRethermalizedEnergy =
		CR_Particle_ShockInject(SphP + i,
					rShockEnergy * All.CR_ShockEfficiency,
					(P[i].Ti_endstep - P[i].Ti_begstep) *
					All.Timebase_interval / hubble_a);


	      /* Fraction of total energy that went and remained in CR is
	       * rNonRethermalizedEnergy / rShockEnergy,
	       * hence, we conserve energy if we do:
	       */
	      SphP[i].DtEntropy *= (1.0 - rNonRethermalizedEnergy / rShockEnergy);

	      assert(rNonRethermalizedEnergy > 0.0);

	      assert(rNonRethermalizedEnergy < (rShockEnergy * All.CR_ShockEfficiency));

	    }
#endif /* CR_SHOCK */

	  SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));

#ifdef MACHNUM
	  /* Estimates the Mach number of particle i */
	  SphP[i].MachNumber = MachNumber(SphP + i);

	  SphP[i].DissipatedEnergy = DissipatedEnergy(SphP + i, P + i);

#endif


#ifdef NAVIERSTOKES
	  for(k = 0, fac = 0; k < 3; k++)
	    {
	      fac += SphP[i].u.s.StressDiag[k] * SphP[i].u.s.StressDiag[k] +
		2 * SphP[i].u.s.StressOffDiag[k] * SphP[i].u.s.StressOffDiag[k];
	    }

	  SphP[i].DtEntropy += 0.5 * GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1)) *
	    get_shear_viscosity(i) / SphP[i].Density * fac;
#endif



#ifdef MAGNETIC
	  /* take care of cosmological dilution */
	  if(All.ComovingIntegrationOn)
	    for(k = 0; k < 3; k++)
	      SphP[i].DtB[k] -= 2.0 * SphP[i].BPred[k];
#endif

#ifdef WINDS
	  /* if we have winds, we decouple particles briefly if delaytime>0 */

	  if(SphP[i].DelayTime > 0)
	    {
	      for(k = 0; k < 3; k++)
		SphP[i].HydroAccel[k] = 0;

	      SphP[i].DtEntropy = 0;
	      SphP[i].MaxSignalVel = 2 * sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
	    }
#endif

#ifdef SPH_BND_PARTICLES
	  if(P[i].ID == 0)
	    {
	      SphP[i].DtEntropy = 0;
	      for(k = 0; k < 3; k++)
		SphP[i].HydroAccel[k] = 0;
	    }
#endif
	}


#ifdef CONDUCTION		/* note: this needs to be done for all particles, not just active ones */
  for(i = 0; i < N_gas; i++)
#ifdef SFR
    if(P[i].Type == 0)
#endif
      {
	if(All.ComovingIntegrationOn)
	  SphP[i].CondEnergyChange *= All.Time / hubble_a;

	if(P[i].Ti_endstep > 0)
	  {
	    rEntropyFlow =
	      SphP[i].CondEnergyChange / P[i].Mass * GAMMA_MINUS1 * pow(SphP[i].Density * a3inv,
									-GAMMA_MINUS1);
	    if(fabs(rEntropyFlow) > 0.5 * SphP[i].Entropy)
	      {
		printf
		  ("NOTE: step=%d task=%d large entropy change for particle i=%d (ID=%d) due to conduction: A=%g  dA=%g rho=%g dt=%g\n",
		   All.NumCurrentTiStep, ThisTask, i, (int) P[i].ID, SphP[i].Entropy, rEntropyFlow,
		   SphP[i].Density, (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval);
		fflush(stdout);

		fac = 0.5 * SphP[i].Entropy / fabs(rEntropyFlow);

		rEntropyFlow *= fac;
		SphP[i].CondEnergyChange *= fac;
	      }

	    SphP[i].Entropy += rEntropyFlow;

	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
		SphP[i].CondEnergyChange /= (P[i].Mass * dt);
	      }
	  }
	else
	  {
	    SphP[i].CondEnergyChange /= P[i].Mass;
	  }
      }
#endif


  tend = second();
  timecomp += timediff(tstart, tend);

  /* collect some timing information */

  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
    }
}


/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void hydro_evaluate(int target, int mode)
{
  int startnode, numngb;
  int j, k, n, timestep;
  FLOAT *pos, *vel;
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
  double acc[3], dtEntropy, maxSignalVel;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j;
  double r, r2, u;
  double hfc_visc;
  double dmin1, dmin2;
  int imax1, imax2;

#ifdef NAVIERSTOKES
  double faci, facj;
  FLOAT *stressdiag;
  FLOAT *stressoffdiag;
  FLOAT shear_viscosity;
#endif
#ifdef CONDUCTION
  double condEnergyChange;
  double rEnergy_i, rEnergy_j, rEnergy_i_ns, rEnergy_j_ns, Entropy;
  double rKappa_i, rKappa_j;	/* Thermal Conduction Coefficients */
  double rEnergyFlow;
  double smoothentr;
  double dtcond;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
  double electron_free_path, temp_scale_length;
#endif
#endif


#if defined(CONDUCTION)
  double rEntropy2Energy_i;
#endif
#if  defined(CR_DIFFUSION)
  double rCR_EnergyChange_i = 0.0;
  double rCR_BaryonFractionChange_i = 0.0;
#endif

#ifdef CR_DIFFUSION
  double rCR_Diff_N_i, rCR_Diff_E_i;

  double rCR_DiffusionTerm;
  double rCR_EnergyFlow;
  double rCR_MassFlow;
#endif


#ifdef REDUCEVISC
  FLOAT alpha;
  double BulkVisc_ij;
#endif
#ifdef ALTVISCOSITY
  double mu_i, mu_j;
#endif
#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

#ifdef MAGNETIC
  FLOAT *bpred;
  double dtB[3];
  double magfac, magfac_i, magfac_j, magfac_i_base;

#ifdef MAGFORCE
  double mm_i[3][3], mm_j[3][3];
  double b2;
  int l;

#ifdef ARTBPRES
  double wk_i, wk_j, hinv3;
  double w1_i, w1_j, R_abp_i, R_abp_j;
  double mm_abp_i[3][3], mm_abp_j[3][3];
#endif
#endif

#ifdef TRACEDIVB
  double divB;
#endif

#endif

#ifdef SFR_DECOUPLING
  FLOAT densityold, entropy;
#endif


#ifdef PARTICLE_DEBUG
  int4byte ID;
#endif

#ifdef CONVENTIONAL_VISCOSITY
  double c_ij, h_ij;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
#ifndef NOGRADHSML
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
#endif
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#ifndef ALTVISCOSITY
      f1 = fabs(SphP[target].u.s.DivVel) /
	(fabs(SphP[target].u.s.DivVel) + SphP[target].u.s.CurlVel +
	 0.0001 * soundspeed_i / PPP[target].Hsml / fac_mu);
#else
      f1 = SphP[target].u.s.DivVel;
#endif
#ifdef MAGNETIC
      bpred = SphP[target].BPred;
#endif
#ifdef REDUCEVISC
      alpha = SphP[target].alpha;
#endif
#ifdef CONDUCTION
      smoothentr = SphP[target].SmoothedEntr;
      Entropy = SphP[target].Entropy;
#ifdef CONDUCTION_SATURATION
      gradentr[0] = SphP[target].GradEntr[0];
      gradentr[1] = SphP[target].GradEntr[1];
      gradentr[2] = SphP[target].GradEntr[2];
#endif
#endif


#ifdef SFR_DECOUPLING
/*      densityold = SphP[target].DensityOld; */
      densityold = SphP[target].Density;
      entropy = SphP[target].Entropy;
#endif

#ifdef PARTICLE_DEBUG
      ID = SphP[target].ID;
#endif

#if CR_DIFFUSION
      rCR_Diff_E_i = SphP[target].CR_Diff_E;
      rCR_Diff_N_i = SphP[target].CR_Diff_N;
#endif
#ifdef NAVIERSTOKES
      stressdiag = SphP[target].u.s.StressDiag;
      stressoffdiag = SphP[target].u.s.StressOffDiag;
      shear_viscosity = get_shear_viscosity(target);
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
#ifndef NOGRADHSML
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
#endif
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = HydroDataGet[target].F1;
#ifdef MAGNETIC
      bpred = HydroDataGet[target].BPred;
#endif
#ifdef REDUCEVISC
      alpha = HydroDataGet[target].alpha;
#endif
#ifdef CONDUCTION
      smoothentr = HydroDataGet[target].SmoothedEntr;
      Entropy = HydroDataGet[target].Entropy;
#ifdef CONDUCTION_SATURATION
      gradentr[0] = HydroDataGet[target].GradEntr[0];
      gradentr[1] = HydroDataGet[target].GradEntr[1];
      gradentr[2] = HydroDataGet[target].GradEntr[2];
#endif
#endif

#ifdef SFR_DECOUPLING
      densityold = HydroDataGet[target].DensityOld;
      entropy = HydroDataGet[target].Entropy;
#endif

#ifdef PARTICLE_DEBUG
      ID = HydroDataGet[target].ID;
#endif

#if CR_DIFFUSION
      rCR_Diff_E_i = HydroDataGet[target].CR_Diff_E;
      rCR_Diff_N_i = HydroDataGet[target].CR_Diff_N;
#endif
#ifdef NAVIERSTOKES
      stressdiag = HydroDataGet[target].stressdiag;
      stressoffdiag = HydroDataGet[target].stressoffdiag;
      shear_viscosity = HydroDataGet[target].shear_viscosity;
#endif
    }


  /* initialize variables before SPH loop is started */

  acc[0] = acc[1] = acc[2] = dtEntropy = 0;

#ifdef CONDUCTION
  condEnergyChange = 0;
  dtcond = timestep * All.Timebase_interval;
  if(timestep == 0)
    dtcond = 1;

#endif

#ifdef MAGNETIC
  for(k = 0; k < 3; k++)
    dtB[k] = 0;
#ifdef TRACEDIVB
  divB = 0;
#endif
#ifdef MAGFORCE
  magfac_i_base = 1 / (rho * rho);
#ifndef BRIOWU
  magfac_i_base /= (4 * M_PI);
#endif
#ifdef CORRECTBFRC
  magfac_i_base *= dhsmlDensityFactor;
#endif
  for(k = 0, b2 = 0; k < 3; k++)
    {
      b2 += bpred[k] * bpred[k];
      for(l = 0; l < 3; l++)
	mm_i[k][l] = bpred[k] * bpred[l];
    }

  for(k = 0; k < 3; k++)
    mm_i[k][k] -= 0.5 * b2;
#endif
#endif /* end of MAGNETIC */
  p_over_rho2_i = pressure / (rho * rho);
#ifndef NOGRADHSML
  p_over_rho2_i *= dhsmlDensityFactor;
#endif
  h_i2 = h_i * h_i;

#if defined(CONDUCTION)
  rEntropy2Energy_i = pow(rho * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif

#ifdef CONDUCTION
  /* this gives the thermal energy per unit mass for particle i */
  rEnergy_i = smoothentr * rEntropy2Energy_i;
  rEnergy_i_ns = Entropy * rEntropy2Energy_i;

#  ifdef CONDUCTION_CONSTANT
  rKappa_i = All.ConductionCoeff;
#  else
  rKappa_i = All.ConductionCoeff * pow(rEnergy_i_ns, 2.5);
#    ifdef CONDUCTION_SATURATION
  electron_free_path = All.ElectronFreePathFactor * rEnergy_i * rEnergy_i / (rho * a3inv);
  temp_scale_length =
    atime * fabs(smoothentr) / sqrt(gradentr[0] * gradentr[0] +
				    gradentr[1] * gradentr[1] + gradentr[2] * gradentr[2]);
  rKappa_i /= (1 + 4.2 * electron_free_path / temp_scale_length);
#    endif
#  endif
#  ifdef SFR
  if(rho * a3inv >= All.PhysDensThresh)
    rKappa_i = 0;
#  endif
#endif


  maxSignalVel = 0;

  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
#ifdef SFR_DECOUPLING
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode, densityold, entropy, &vel[0]);
#else
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);
#endif
      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

#ifdef BLACK_HOLES
	  if(P[j].Mass == 0)
	    continue;
#endif

#ifdef WINDS
	  if(P[j].Type == 0)
	    if(SphP[j].DelayTime > 0)	/* ignore the wind particles */
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
	  h_j = PPP[j].Hsml;
	  if(r2 < h_i2 || r2 < h_j * h_j)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{
		  p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
		  soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;

		  if(All.ComovingIntegrationOn)
		    vdotr2 = vdotr + hubble_a2 * r2;
		  else
		    vdotr2 = vdotr;

		  if(r2 < h_i2)
		    {
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
#ifndef  TWODIMS
		      hinv3 = hinv * hinv * hinv;
#else
		      hinv3 = hinv * hinv / boxSize_Z;
#endif
		      if(u <= 0.5)
			wk_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      else
			wk_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      if(u1 <= 0.5)
			w1_i = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u1 - 1) * u1 * u1);
		      else
			w1_i = hinv3 * KERNEL_COEFF_5 * (1.0 - u1) * (1.0 - u1) * (1.0 - u1);
#endif
		    }
		  else
		    {
		      dwk_i = 0;
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
		      wk_i = 0;
		      w1_i = 1;
#endif
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
#ifndef  TWODIMS
		      hinv3 = hinv * hinv * hinv;
#else
		      hinv3 = hinv * hinv / boxSize_Z;
#endif
		      if(u <= 0.5)
			wk_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		      else
			wk_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		      if(u1 <= 0.5)
			w1_j = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u1 - 1) * u1 * u1);
		      else
			w1_j = hinv3 * KERNEL_COEFF_5 * (1.0 - u1) * (1.0 - u1) * (1.0 - u1);
#endif
		    }
		  else
		    {
		      dwk_j = 0;
#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
		      wk_j = 0;
		      w1_j = 1.;
#endif
		    }

#if defined(MAGNETIC) && defined(MAGFORCE) && defined(ARTBPRES)
		  R_abp_i = 0.3 * pow(wk_i / w1_i, 4);
		  R_abp_j = 0.3 * pow(wk_j / w1_j, 4);
#endif

		  if(soundspeed_i + soundspeed_j > maxSignalVel)
		    maxSignalVel = soundspeed_i + soundspeed_j;

		  if(vdotr2 < 0)	/* ... artificial viscosity */
		    {
#ifndef ALTVISCOSITY
#ifndef CONVENTIONAL_VISCOSITY
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
#else
		      c_ij = 0.5 * (soundspeed_i + soundspeed_j);
		      h_ij = 0.5 * (h_i + h_j);
		      mu_ij = fac_mu * h_ij * vdotr2 / (r2 + 0.0001 * h_ij * h_ij);
#endif
		      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
		      f2 =
			fabs(SphP[j].u.s.DivVel) / (fabs(SphP[j].u.s.DivVel) + SphP[j].u.s.CurlVel +
						    0.0001 * soundspeed_j / fac_mu / PPP[j].Hsml);
#ifdef REDUCEVISC
		      BulkVisc_ij = 0.5 * (alpha + SphP[j].alpha);
#ifndef CONVENTIONAL_VISCOSITY
		      visc = 0.25 * BulkVisc_ij * vsig * (-mu_ij) / rho_ij * (f1 + f2);
#else
		      visc =
			(-BulkVisc_ij * mu_ij * c_ij + 2 * BulkVisc_ij * mu_ij * mu_ij) /
			rho_ij * (f1 + f2) * 0.5;
#endif
#else
#ifndef CONVENTIONAL_VISCOSITY
		      visc = 0.25 * All.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);
#else
		      visc =
			(-All.ArtBulkViscConst * mu_ij * c_ij +
			 2 * All.ArtBulkViscConst * mu_ij * mu_ij) / rho_ij * (f1 + f2) * 0.5;
#endif
#endif

#else /* start of ALTVISCOSITY block */
		      if(f1 < 0)
			mu_i = h_i * fabs(f1);	/* f1 hold here the velocity divergence of particle i */
		      else
			mu_i = 0;
		      if(SphP[j].u.s.DivVel < 0)
			mu_j = h_j * fabs(SphP[j].u.s.DivVel);
		      else
			mu_j = 0;
		      visc = All.ArtBulkViscConst * ((soundspeed_i + mu_i) * mu_i / rho +
						     (soundspeed_j + mu_j) * mu_j / SphP[j].Density);
#endif /* end of ALTVISCOSITY block */


		      /* .... end artificial viscosity evaluation */
		      /* now make sure that viscous acceleration is not too large */

#ifndef NOVISCOSITYLIMITER
		      dt = 2 * IMAX(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		      if(dt > 0 && (dwk_i + dwk_j) < 0)
			{
#ifdef BLACK_HOLE
			  if((mass + P[j].Mass) > 0)
#endif
			    visc = DMIN(visc, 0.5 * fac_vsic_fix * vdotr2 /
					(0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			}
#endif
		    }
		  else
		    visc = 0;
#ifndef NOGRADHSML
		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
#endif
		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
		  /* Formulation derived from the Lagrangian */
		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#ifndef NOACCEL
		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;
#endif

		  dtEntropy += 0.5 * hfc_visc * vdotr2;


#ifdef NAVIERSTOKES
		  faci = mass * shear_viscosity / (rho * rho) * dwk_i / r;
		  facj = P[j].Mass * get_shear_viscosity(j) / (SphP[j].Density * SphP[j].Density) * dwk_j / r;

		  acc[0] += faci * (stressdiag[0] * dx + stressoffdiag[0] * dy + stressoffdiag[1] * dz)
		    + facj * (SphP[j].u.s.StressDiag[0] * dx + SphP[j].u.s.StressOffDiag[0] * dy +
			      SphP[j].u.s.StressOffDiag[1] * dz);

		  acc[1] += faci * (stressoffdiag[0] * dx + stressdiag[1] * dy + stressoffdiag[2] * dz)
		    + facj * (SphP[j].u.s.StressOffDiag[0] * dx + SphP[j].u.s.StressDiag[1] * dy +
			      SphP[j].u.s.StressOffDiag[2] * dz);

		  acc[2] += faci * (stressoffdiag[1] * dx + stressoffdiag[2] * dy + stressdiag[2] * dz)
		    + facj * (SphP[j].u.s.StressOffDiag[1] * dx + SphP[j].u.s.StressOffDiag[2] * dy +
			      SphP[j].u.s.StressDiag[2] * dz);
#endif







#ifdef CONDUCTION
		  rEnergy_j =
		    SphP[j].SmoothedEntr * pow(SphP[j].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
		  rEnergy_j_ns = SphP[j].Entropy * pow(SphP[j].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#ifdef CONDUCTION_CONSTANT
		  rKappa_j = All.ConductionCoeff;
#else
		  rKappa_j = All.ConductionCoeff * pow(rEnergy_j_ns, 2.5);
#ifdef CONDUCTION_SATURATION
		  electron_free_path =
		    All.ElectronFreePathFactor * rEnergy_j * rEnergy_j / (SphP[j].Density * a3inv);
		  temp_scale_length =
		    atime * fabs(smoothentr) / sqrt(SphP[j].GradEntr[0] * SphP[j].GradEntr[0] +
						    SphP[j].GradEntr[1] * SphP[j].GradEntr[1] +
						    SphP[j].GradEntr[2] * SphP[j].GradEntr[2]);
		  rKappa_j /= (1 + 4.2 * electron_free_path / temp_scale_length);
#endif
#endif
#ifdef SFR
		  if(SphP[j].Density * a3inv >= All.PhysDensThresh)
		    rKappa_j = 0;
#endif
		  if((rKappa_i + rKappa_j) > 0)
		    {
		      rEnergyFlow =
			(mass * P[j].Mass) /
			(rho * SphP[j].Density) *
			4 * (rKappa_i * rKappa_j) /
			(rKappa_i + rKappa_j) * (rEnergy_j - rEnergy_i) * 0.5 * (dwk_i + dwk_j) / r * dtcond;
		    }
		  else
		    {
		      rEnergyFlow = 0;
		    }

		  if(timestep > 0)
		    SphP[j].CondEnergyChange += 0.5 * rEnergyFlow;
		  condEnergyChange += -0.5 * rEnergyFlow;
#endif /* end of CONDUCTION */

#ifdef CR_DIFFUSION
		  /* *************** COSMIC RAY DIFFUSION ********************
		   */

		  /* In the following, the cosmic ray diffusion is computed.
		   * First, the total energy/mass transfer during particle i's
		   * timestep are evaluated, then the corresponding amount of
		   * specific CR energy and baryon fraction are adjusted in
		   * both particle i and particle j
		   */

		  /* Compute the common term element that is typical for
		   * conduction/diffusion effects. So we only have to do
		   * some of the computations once. 
		   */
		  rCR_DiffusionTerm =
		    -(mass * P[j].Mass) /
		    (rho * SphP[j].Density) *
		    0.5 * (dwk_i + dwk_j) / r * timestep * All.Timebase_interval * atime / hubble_a;

		  /* TODO: Check whether atime/hubble_a is the right
		   *       conversion from comoving time to physical time.
		   */


		  /* Compute the flow of total CR energy 
		   * (conserved quantity) 
		   */
		  rCR_EnergyFlow =
		    rCR_DiffusionTerm * 2.0 * All.CR_DiffusionCoeff * (SphP[j].CR_Diff_E - rCR_Diff_E_i);

		  /* Compute the flow of total CR mass transferred
		   * (conserved quantity) 
		   */
		  rCR_MassFlow =
		    rCR_DiffusionTerm * 2.0 * All.CR_DiffusionCoeff * (SphP[j].CR_Diff_N - rCR_Diff_N_i);

		  /* Subtract half of the flow from particle j 
		   */
		  CR_Particle_Inject(SphP + j,
				     -rCR_EnergyFlow * 0.5 / P[j].Mass, -rCR_MassFlow * 0.5 / P[j].Mass);

		  /*And add it to particle i 
		   */
		  rCR_EnergyChange_i += rCR_EnergyFlow * 0.5 / mass;
		  rCR_BaryonFractionChange_i += rCR_MassFlow * 0.5 / mass;

		  /* ************ COSMIC RAY DIFFUSION - END *****************
		   */
#endif


#ifdef MAGNETIC
		  magfac = P[j].Mass / rho * dwk_i / r;
#ifdef TRACEDIVB
		  divB += magfac *
		    ((SphP[j].BPred[0] - bpred[0]) * dx +
		     (SphP[j].BPred[1] - bpred[1]) * dy + (SphP[j].BPred[2] - bpred[2]) * dz);
#endif
		  if(All.ComovingIntegrationOn)
		    magfac *= 1. / (hubble_a * All.Time * All.Time);
		  /* last factor takes care of all cosmological prefactor */

#ifdef CORRECTDB
		  magfac *= dhsmlDensityFactor;
#endif
		  dtB[0] +=
		    magfac * ((bpred[0] * dvy - bpred[1] * dvx) * dy +
			      (bpred[0] * dvz - bpred[2] * dvx) * dz);
		  dtB[1] +=
		    magfac * ((bpred[1] * dvz - bpred[2] * dvy) * dz +
			      (bpred[1] * dvx - bpred[0] * dvy) * dx);
		  dtB[2] +=
		    magfac * ((bpred[2] * dvx - bpred[0] * dvz) * dx +
			      (bpred[2] * dvy - bpred[1] * dvz) * dy);

#ifdef MAGFORCE
		  magfac_j = 1 / (SphP[j].Density * SphP[j].Density);
#ifndef BRIOWU
		  magfac_j /= (4 * M_PI);
#endif
#ifdef CORRECTBFRC
		  magfac_j *= dwk_j * SphP[j].DhsmlDensityFactor;
		  magfac_i = dwk_i * magfac_i_base;
#else
		  magfac_i = magfac_i_base;
#endif
		  for(k = 0, b2 = 0; k < 3; k++)
		    {
		      b2 += SphP[j].BPred[k] * SphP[j].BPred[k];
		      for(l = 0; l < 3; l++)
			mm_j[k][l] = SphP[j].BPred[k] * SphP[j].BPred[l];
		    }

		  for(k = 0; k < 3; k++)
		    mm_j[k][k] -= 0.5 * b2;
#ifdef CORRECTBFRC
		  magfac = P[j].Mass / r;
#else
		  magfac = P[j].Mass * 0.5 * (dwk_i + dwk_j) / r;
#endif
		  if(All.ComovingIntegrationOn)
		    magfac *= pow(All.Time, 3 * GAMMA);
		  /* last factor takes care of all cosmological prefactor */
#ifndef BRIOWU
		  magfac *= All.UnitTime_in_s * All.UnitTime_in_s *
		    All.UnitLength_in_cm / (All.UnitMass_in_g * All.HubbleParam * All.HubbleParam);
		  /* take care of B unit conversion into GADGET units ! */
#endif
		  for(k = 0; k < 3; k++)
		    acc[k] +=
		      magfac * ((mm_i[k][0] * magfac_i + mm_j[k][0] * magfac_j) * dx +
				(mm_i[k][1] * magfac_i + mm_j[k][1] * magfac_j) * dy +
				(mm_i[k][2] * magfac_i + mm_j[k][2] * magfac_j) * dz);
#ifdef ARTBPRES
		  for(k = 0; k < 3; k++)
		    for(l = 0; l < 3; l++)
		      mm_abp_j[k][l] = SphP[j].BPred[k] * SphP[j].BPred[l] * R_abp_j / 2.0;
		  for(k = 0; k < 3; k++)
		    for(l = 0; l < 3; l++)
		      mm_abp_i[k][l] = bpred[k] * bpred[l] * R_abp_i / 2.0;
		  for(k = 0; k < 3; k++)
		    acc[k] -= magfac * ((mm_abp_i[k][0] * magfac_i +
					 mm_abp_j[k][0] * magfac_j) * dx +
					(mm_abp_i[k][1] * magfac_i +
					 mm_abp_j[k][1] * magfac_j) * dy +
					(mm_abp_i[k][2] * magfac_i + mm_abp_j[k][2] * magfac_j) * dz);
#endif
#ifdef DIVBFORCE
		  for(k = 0; k < 3; k++)
		    acc[k] -=
		      magfac * (((bpred[k] * bpred[0]) * magfac_i
				 + (bpred[k] * SphP[j].BPred[0]) * magfac_j) * dx
				+ ((bpred[k] * bpred[1]) * magfac_i
				   + (bpred[k] * SphP[j].BPred[1]) * magfac_j) * dy
				+ ((bpred[k] * bpred[2]) * magfac_i +
				   (bpred[k] * SphP[j].BPred[2]) * magfac_j) * dz);
#endif
#endif
#endif /* end of MAGNETIC */
		}
	    }
	}
    }
  while(startnode >= 0);
  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].HydroAccel[k] = acc[k];
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;
#ifdef CONDUCTION
      SphP[target].CondEnergyChange += condEnergyChange;
#ifdef OUTPUTCOOLRATE
      SphP[target].CondRate = 2 * atime * condEnergyChange / (dtcond * mass);
#endif
#endif
#ifdef MAGNETIC
      for(k = 0; k < 3; k++)
	SphP[target].DtB[k] = dtB[k];
#ifdef TRACEDIVB
      SphP[target].divB = divB;
#endif
#endif

#if defined(CR_DIFFUSION)
      CR_Particle_Inject(SphP + target, rCR_EnergyChange_i, rCR_BaryonFractionChange_i);
#endif

    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
#ifdef CONDUCTION
      HydroDataResult[target].CondEnergyChange = condEnergyChange;
#ifdef OUTPUTCOOLRATE
      HydroDataResult[target].CondRate = 2 * atime * condEnergyChange / (dtcond * mass);
#endif
#endif
#ifdef MAGNETIC
      for(k = 0; k < 3; k++)
	HydroDataResult[target].DtB[k] = dtB[k];
#ifdef TRACEDIVB
      HydroDataResult[target].divB = divB;
#endif
#endif
#if defined( CR_DIFFUSION )
      HydroDataResult[target].CR_EnergyChange = rCR_EnergyChange_i;
      HydroDataResult[target].CR_BaryonFractionChange = rCR_BaryonFractionChange_i;
#endif
    }
}




/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int hydro_compare_key(const void *a, const void *b)
{
  if(((struct hydrodata_in *) a)->Task < (((struct hydrodata_in *) b)->Task))
    return -1;
  if(((struct hydrodata_in *) a)->Task > (((struct hydrodata_in *) b)->Task))
    return +1;
  return 0;
}


#ifdef MACHNUM

/* Part used for the Mach number finder (by Christoph Pfrommer): */

/* --- function EntropieJump --- */
double EntropyJump(double gamma, double M)
{
  double f;

  f = (1 + gamma * (2 * M * M - 1)) / (gamma + 1) *
    pow(((gamma - 1) * M * M + 2) / ((gamma + 1) * M * M), gamma);
  return (f);
}

/* --- end of function EntropieJump --- */

/* --- function DEntropieJumpDM: derivative of EntropieJump with respect to to M --- */
double DEntropyJumpDM(double gamma, double M)
{
  double f;

  f = 4 * (gamma - 1) * gamma * pow(M * M - 1, 2) / (2 * (gamma + 1) * M + (pow(gamma, 2) - 1) * M * M * M) *
    pow(((gamma - 1) * M * M + 2) / ((gamma + 1) * M * M), gamma);
  return (f);
}

/* --- end of function DEntropieJumpDM --- */

/* --- function AuxMach --- */
void AuxMach(double gamma, double M, double *f, double *df)
{
  double DelS;

  DelS = EntropyJump(gamma, M);

  *f = (DelS - 1) * M;
  *df = DEntropyJumpDM(gamma, M) * M + DelS - 1;
}

/* --- end of function AuxMach --- */


/*  function MachNumber: computes the Mach number per SPH particle given 
 *  the SPH smoothing length, the local sound speed, the total entropic function A, 
 *  and the instantaneous rate of entropy injection due to shocks dA/dt.
 *  It uses the Newton-Raphson method to find the root of following equation:
 *  [f(M) - 1] * M = h / (c_sound * A) * dA/dt,
 *  where f(M) = A_2 / A_1 (ratio of post-shock to pre-shock entropic function)
 */

#ifdef COSMIC_RAYS

double MachNumber(struct sph_particle_data *Particle)
{
  double M, dM, r;
  double f, df, rhs, csnd;
  double gamma, PCR, Pgas, Ugas, TotalEntropy, TotalDtEntropy;
  int i = 0;

  /* Particle->Pressure = P_CR + P_gas (comoving) */
  /* Particle->Entropy = entropic function of the gas = exp(s / c_V) */
  PCR = CR_Particle_Pressure(Particle);	/* (comoving) */
  Pgas = Particle->Pressure - PCR;

  /* effective adiabatic index gamma */
  gamma = (GAMMA * Pgas + Particle->CR_Gamma0 * PCR) / Particle->Pressure;

  /* effective sound velocity */
  csnd = sqrt(gamma * Particle->Pressure / Particle->Density);

  /* thermal specific energy of the gas (in physical units) */
  Ugas = Particle->Entropy * pow(Particle->Density, GAMMA - 1.0) / (GAMMA - 1.0) / atime / atime;

  /* total entropic function A */
  TotalEntropy = (gamma - 1.0) / pow(Particle->Density, gamma - 1.0) *
    (CR_Particle_SpecificEnergy(Particle) + Ugas);

  /* rate of change of the total entropic function, dA/dt */
  TotalDtEntropy = (gamma - 1.0) / (GAMMA - 1.0) * Particle->DtEntropy *
    pow(Particle->Density, GAMMA - gamma) * pow(atime, -3.0 * (GAMMA - gamma));

  rhs = Particle->Hsml / (csnd * TotalEntropy) * TotalDtEntropy * hubble_a * atime * atime;

  /* initial guess */
  if(rhs > 0)
    M = 2;
  else
    M = 0.5;

  /* Newton-Raphson scheme */
  r = 1.0;
  while((i < 100) && (r > 1e-4))
    {
      AuxMach(gamma, M, &f, &df);
      dM = (f - rhs) / df;
      M -= dM;
      r = fabs(dM / M);
      i++;
    }
  if(i == 100)
    printf("Error: too many iterations in function MachNumber!\n");
  return (M);
}

/* --- end of function MachNumber --- */

#else

/* --- function MachNumber: for the thermal gas without cosmic rays --- */
double MachNumber(struct sph_particle_data *Particle)
{
  double M, dM, r;
  double f, df, rhs, csnd;
  int i = 0;

  /* sound velocity */
  csnd = sqrt(GAMMA * Particle->Pressure / Particle->Density);

  rhs = Particle->Hsml / (csnd * Particle->Entropy) * Particle->DtEntropy * hubble_a * atime * atime;

  /* initial guess */
  if(rhs > 0)
    M = 2.;
  else
    M = 0.5;

  /* Newton-Raphson scheme */
  r = 1.0;
  while((i < 100) && (r > 1e-4))
    {
      AuxMach(GAMMA, M, &f, &df);
      dM = (f - rhs) / df;
      M -= dM;
      r = fabs(dM / M);
      i++;
    }
  if(i == 100)
    printf("Error: too many iterations in function MachNumber!\n");
  return (M);
}

/* --- end of function MachNumber --- */

#endif

/* --- function: DissipatedEnergy (actually dissipated specific energy)  --- */
double DissipatedEnergy(struct sph_particle_data *Particle, struct particle_data *Properties)
{
  double u;

  u = Particle->DtEntropy *
    pow(Particle->Density, GAMMA_MINUS1) / (GAMMA_MINUS1 * atime * atime) *
    (Properties->Ti_endstep - Properties->Ti_begstep) * All.Timebase_interval;

  return (u);
}

/* --- end of function DissipatedEnergy --- */

#endif
