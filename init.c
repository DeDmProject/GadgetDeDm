#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef SFR_METALS
#include "c_metals.h"
#endif

/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;

#ifdef CHEMISTRY
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
    case 2:
      read_ic(All.InitCondFile);
      break;
    case 4:
      if(RestartFlag == 2)
	read_ic(All.InitCondFile);
#ifdef IC_CLUSTER
      else
	read_ic_cluster(All.InitCondFile);
#endif
      break;
    case 5:			/* this is a special format for bepis' cluster ICs */
      if(RestartFlag == 2)
	read_ic(All.InitCondFile);
#ifdef IC_CLUSTER
      else
	read_ic_cluster_gas(All.InitCondFile);
#endif
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;

#ifdef COOLING
#ifdef SFR_METALS
  if(RestartFlag == 0)
    XH = HYDROGEN_MASSFRAC;
#endif
  IonizeParams();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
  All.MainTimestepCounts = 0;
#endif

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#ifdef BLACK_HOLES
  All.TimeNextBlackHoleCheck = All.TimeBegin;
#endif



#ifdef BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + 3.);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  printf("Initial time: %g and first bubble time %g \n", 1. / All.TimeBegin - 1,
	 1. / All.TimeOfNextBubble - 1);
#endif

#ifdef MULTI_BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + 3.);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  printf("Initial time: %g and time of the first bubbles %g \n", 1. / All.TimeBegin - 1, 
	 1. / All.TimeOfNextBubble-1);
#endif

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
#ifdef STELLARAGE
      if(RestartFlag == 0)
	P[i].StellarAge = 0;
#endif
#ifdef METALS
#ifndef SFR_METALS
      if(RestartFlag == 0)
	P[i].Metallicity = 0;
#else
      if(RestartFlag == 0)
	{
	  for(j = 0; j < 12; j++)
	    {
	      P[i].Zm[j] = 0;
	      P[i].ZmReservoir[j] = 0;
	    }

	  P[i].Zm[6] = HYDROGEN_MASSFRAC * (P[i].Mass);
	  P[i].Zm[0] = (1 - HYDROGEN_MASSFRAC) * (P[i].Mass);
	  /*     Order of chemical elements:   He, Carbon,Mg,O,Fe,Si,H,N,Ne,S,Ca,Zn */

	  PPP[i].Hsml = 0;
#ifdef SFR_FEEDBACK
	  PPP[i].EnergySN = 0;
	  PPP[i].EnergySNCold = 0;
#endif
	}
#endif

#endif

#ifdef BLACK_HOLES
      if(RestartFlag == 0 && P[i].Type == 5)
	P[i].BH_Mass = All.SeedBlackHoleMass;
#ifdef BH_KINETICFEEDBACK
      if(P[i].Type == 5)
	{
	  P[i].ActiveTime = 0;
	  P[i].ActiveEnergy = 0;
	}
#endif
#endif
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif


#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

      SphP[i].DtEntropy = 0;

#ifdef CHEMISTRY
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif

#ifdef SFR_DECOUPLING
      SphP[i].DensityOld = 0;
#endif
#ifdef SFR_PROMOTION
      SphP[i].DensityAvg = 0;
      SphP[i].EntropyAvg = 0;
      SphP[i].DensPromotion = 0;
      SphP[i].TempPromotion = 0;

#endif

      if(RestartFlag == 0)
	{
	  PPP[i].Hsml = 0;
	  SphP[i].Density = -1;
#ifdef COOLING
	  SphP[i].Ne = 1.0;
#endif
	}
#ifdef WINDS
      SphP[i].DelayTime = 0;
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#ifdef MAGNETIC
      for(j = 0; j < 3; j++)
	{
	  SphP[i].DtB[j] = 0;
	  SphP[i].BPred[j] = SphP[i].B[j];
	}
#ifdef BINISET
      SphP[i].B[0] = All.BiniX;
      SphP[i].B[1] = All.BiniY;
      SphP[i].B[2] = All.BiniZ;
      SphP[i].BPred[0] = All.BiniX;
      SphP[i].BPred[1] = All.BiniY;
      SphP[i].BPred[2] = All.BiniZ;
#endif
#endif
#ifdef REDUCEVISC
#ifdef HIGHVISCSTART
      if(HIGHVISCSTART == 0)
	SphP[i].alpha = All.ArtBulkViscConst;
      if(HIGHVISCSTART > 0)
	if(P[i].Pos[0] > HIGHVISCSTART)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
      if(HIGHVISCSTART < 0)
	if(P[i].Pos[0] < -HIGHVISCSTART)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
      SphP[i].Injected_BH_Energy = 0;
#endif
    }

#ifdef HPM
  All.HPM_entr0 = All.HPM_entr1 = All.HPM_ne0 = All.HPM_ne1 = 0;
#endif

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  DomainDecomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;


  setup_smoothinglengths();

  /* at this point, the entropy variable actually contains the 
   * internal energy, read in from the initial conditions file. 
   * Once the density has been computed, we can convert to entropy.
   */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      if(header.flag_entropy_instead_u == 0)
	{
	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
	}
      SphP[i].DtEntropy = 0;

#ifdef CR_IC_PHYSICAL
      /* Scale CR variables so that values from IC file are now the
       * physical values, not the adiabatic invariants
       */

      SphP[i].CR_C0 *= pow(SphP[i].Density, (1.0 - All.CR_Alpha) / 3.0);
      SphP[i].CR_q0 *= pow(SphP[i].Density, -1.0 / 3.0);
#endif
    }


  /* savepositions(0); 
     endrun(0); */
#ifdef CHEMISTRY

  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif



}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, j, N;
  double xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];
  double vol, meanspacing, meanspacingtot, hsmlguess;


  if(RestartFlag == 0)
    {

      for(j = 0; j < 3; j++)	/* find enclosing rectangle */
	xmin[j] = xmax[j] = P[0].Pos[j];

      for(i = 0; i < N_gas; i++)
	for(j = 0; j < 3; j++)
	  {
	    if(P[i].Pos[j] > xmax[j])
	      xmax[j] = P[i].Pos[j];
	    if(P[i].Pos[j] < xmin[j])
	      xmin[j] = P[i].Pos[j];
	  }

      MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if(All.TotN_gas)
	{
	  for(j = 0, vol = 1; j < NUMDIMS; j++)	/* find enclosing rectangle */
	    vol *= xmax_glob[j] - xmin_glob[j];

	  meanspacing = pow(vol / All.TotN_gas, 1.0 / NUMDIMS);
	}
      else
	meanspacing = 25.0 * All.SofteningGas;

      MPI_Allreduce(&meanspacing, &meanspacingtot, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

#ifdef BLACK_HOLES
      N = NumPart;
#else
      N = N_gas;
#endif

      if(All.ComovingIntegrationOn)
	{
	  for(i = 0; i < N; i++)
	    PPP[i].Hsml =
	      pow(2 * All.G * All.DesNumNgb * P[i].Mass / (All.OmegaBaryon * All.Hubble * All.Hubble),
		  1.0 / 3);
	}
      else
	{
#ifndef TWODIMS
	  hsmlguess = 0.2 * pow(3 * All.DesNumNgb / (4 * M_PI), 1.0 / 3) * meanspacingtot;
#else
	  hsmlguess = 0.2 * pow(All.DesNumNgb / M_PI, 1.0 / 2) * meanspacingtot;
#endif
	  for(i = 0; i < N; i++)
	    PPP[i].Hsml = hsmlguess;
	}
    }

#ifdef BLACK_HOLES
  if(RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif


  density();
#if defined(MAGNETIC) && defined(BSMOOTH)
  bsmooth();
#endif
}
