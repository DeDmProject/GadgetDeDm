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
#include "c_metals.h"
#endif

#ifdef DEDM_HUBBLE
#include "libdedm/interpolate.h"
#endif

#ifdef VDE
#include "libvde/interpolate.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fully usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  double tstart, tend;

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
  double hubble_a;
#endif

#ifdef SFR_METALS
  int i;
  double u_i, a3inv;

  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  TotalEnergy = 0;
  DEnergy_spawned = 0;
  DEnergy_converted = 0;
  DEnergy_radiation = 0;
  DEnergy_promotion = 0;
  DEnergy_feedback = 0;


/*Not only for active particles! */
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
	u_i = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
	TotalEnergy += u_i * P[i].Mass;
      }

#endif

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      tstart = second();
      long_range_force();
      tend = second();
      All.CPU_PM += timediff(tstart, tend);
    }
#endif

#ifndef ONLY_PM

  tstart = second();		/* measure the time for the full force computation */

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
  tend = second();
  All.CPU_Gravity += timediff(tstart, tend);

#endif


#ifdef FORCETEST
  gravity_forcetest();
#endif

#ifndef HPM

  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

      tstart = second();

#ifdef SFR_METALS
#ifdef SFR_ENRICH		/* equivalent to SNI OR SNII */
      flag_SN_starparticles();	/* mark SNI star particles */
#endif
#endif

#if defined(SFR_METALS) && defined(SFR_DECOUPLING)
      copy_densities();
      find_low_density_tail();
#endif
      density();		/* computes density, and pressure */


#if defined(MAGNETIC) && defined(BSMOOTH)
      bsmooth();
#endif


#if defined(SFR_METALS) && defined(SFR_DECOUPLING) && defined(SFR_PROMOTION)
      find_hot_neighbours();
      promote_particles();
      copy_densities();		/* optional ? */
      density();
#endif

#ifdef CONDUCTION
      conduction_smoothed_temperature();
#endif

      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

      tstart = second();
      force_update_hmax();
      tend = second();

      All.CPU_Predict += timediff(tstart, tend);


      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

      tstart = second();

      hydro_force();		/* adds hydrodynamical accelerations 
				   and computes du/dt  */
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

#ifdef MHM
      kinetic_feedback_mhm();
#endif


#ifdef BLACK_HOLES
#ifdef FOF
      /* this will find new black hole seed halos */
      if(All.Time >= All.TimeNextBlackHoleCheck)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    All.TimeNextBlackHoleCheck *= All.TimeBetBlackHoleSearch;
	  else
	    All.TimeNextBlackHoleCheck += All.TimeBetBlackHoleSearch;
	}
#endif
      blackhole_accretion();
#endif



#ifdef COOLING
      tstart = second();

      cooling_and_starformation();	/* do radiative cooling and star formation */

#ifdef SFR_METALS
#ifdef SFR_ENRICH
#ifndef SFR_FEEDBACK
      if(ThisTask == 0)
	{
	  printf("... updating weights ...\n");
	  fflush(stdout);
	}
      update_weights();

      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation...\n");
	  fflush(stdout);
	}
      enrichment();		/* distribute metals within neighbourhood */

#else
/* HOT phase Flag_phase =1*/
      Flag_phase = 1;
      update_weights();
      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation phase I...\n");
	  fflush(stdout);
	}

      enrichment();		/* distribute metals within neighbourhood */

/* COLD phase Flag_phase =2*/
      Flag_phase = 2;
      update_weights();
      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation phase II..Cleaning...\n");
	  fflush(stdout);
	}

      enrichment();		/* distribute metals within neighbourhood */

      Flag_phase = 0;

#endif
#endif

#ifdef SFR_FEEDBACK
      /*      phase_mass(); */
      energy_test();
#endif

#endif
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);
      All.CPU_SfrCool += timediff(tstart, tend);
#endif



#ifdef BUBBLES
      if(All.Time >= All.TimeOfNextBubble)
	{
#ifdef FOF
	  fof_fof(-1);
	  bubble();
#else
	  bubble();
#endif
	  if(All.ComovingIntegrationOn)
	    {
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
#endif // ifdef Dedm
	      hubble_a = All.Hubble * sqrt(hubble_a);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif


#if defined(MULTI_BUBBLES) && defined(FOF)
      if(All.Time >= All.TimeOfNextBubble)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    {
#if defined(DEDM_HUBBLE) || defined(VDE)
		    hubble_a = getH_a(All.Time);
#else
	      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
		+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
		+ DarkEnergy_a(All.Time);
#else
		+ All.OmegaLambda;
#endif
#endif
	      hubble_a = All.Hubble * sqrt(hubble_a);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;

	  printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif



    }

#endif /* end of HPM */

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
