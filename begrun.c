#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

#ifdef COSMIC_RAYS
#include "libbaryons/cosmic_rays.h"
#endif

#ifdef VDE
#include "libvde/vdevars.h"
#include "libvde/read_tables.h"
#include "libvde/interpolate.h"
#endif

#if defined(DEDM_HUBBLE) || defined(DEDM_MASS) || defined(DEDM_DRAG)
#include "libdedm/interpolate.h"
#include "libdedm/read_tables.h"
#include "libdedm/dedmvars.h"
#endif

#ifdef DEDM_DRAG
#include "libscott/read_scott_tables.h"
#include "libdedm/integrate.h"
#endif 

#ifdef DEDM_MASS
#include "libdedm/mass.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif

/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

  if(ThisTask == 0)
    {
      printf("\nThis is P-Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

#ifdef DEDM_INFO
	char outDedm[100];
	char *time = system("date");
	sprintf(outDedm, "%s%s", All.OutputDir,"dedm_info.txt");
		
		// The file is being opened here in write mode, then closed and reopened in append mode
	All.outDeDmFile = fopen(outDedm, "w");

	if(All.outDeDmFile == NULL) {
		  if(ThisTask == 0)
			fprintf(stdout, "\nError. Unable to open file:%s to store dedm log info.\n", outDedm);
	} else {
  		if(ThisTask == 0)
			fprintf(stdout, "\nSaving dedm log info to%s\n", outDedm);
	}

	fprintf(All.outDeDmFile,"%s\n\n",time);
	fclose(All.outDeDmFile);

	All.outDeDmFile = fopen(outDedm, "a");
#endif

#ifdef DEDM_COUPLING
      setBetaZero();
#endif

#if defined(DEDM_HUBBLE) || defined(DEDM_MASS) || defined(VDE)
  read_all_interpolation_tables();
#ifdef DEDM_DRAG
  load_phidot_spline();
#endif
#endif

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

  allocate_commbuffers();	/* ... allocate buffer-memory for particle 
				   exchange during force computation */
  set_units();

#ifdef DEDM_COUPLING
setG_tilde();
#endif

#ifdef COOLING
  All.Time = All.TimeBegin;
  InitCool();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

#ifdef SFR
#ifndef SFR_METALS
  init_clouds();
#else
  Flag_phase = 0;
  Flag_promotion = 0;
#endif
  Stars_converted = 0;
#endif

#ifdef PERIODIC
  ewald_init();
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef REDUCEVISC
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
  All.DecayTime = 1 / All.DecayLength * sqrt((GAMMA - 1) / 2 * GAMMA);
#endif

  open_outputfiles();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef PMGRID
  long_range_init();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2)
    {
      set_random_numbers();

      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets 
				   all variables in the struct `All'. 
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to 
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.  
				 */

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.BunchSizeForce = all.BunchSizeForce;
      All.BunchSizeDensity = all.BunchSizeDensity;
      All.BunchSizeHydro = all.BunchSizeHydro;
      All.BunchSizeDomain = all.BunchSizeDomain;
#ifdef SFR_METALS
      All.BunchSizeMetal = all.BunchSizeMetal;
#ifdef SFR_FEEDBACK
      All.BunchSizeHotNgbs = all.BunchSizeHotNgbs;
#endif
#endif
#ifdef BLACK_HOLES
      All.BunchSizeBlackhole = all.BunchSizeBlackhole;
#endif
#ifdef FOF
      All.BunchSizeFoF = all.BunchSizeFoF;
#endif
#ifdef MHM
      All.BunchSizeKinetic = all.BunchSizeKinetic;
#endif
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);

#ifdef REDUCEVISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);

#ifdef NO_TREEDATA_IN_RESTART
      /* if this is not activated, the tree was stored in the restart-files,
         which also allocated the storage for it already */

      ngb_treeallocate(MAX_NGB);
      force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

      /* ensures that domain reconstruction will be done and new tree will be constructed */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
#endif
    }

#ifdef PMGRID
  long_range_init_regionsize();
#endif


#ifdef COSMIC_RAYS
  CR_initialize_beta_tabs(All.CR_Alpha);
#endif

  if(All.ComovingIntegrationOn) {
    init_drift_table();

#ifdef DEDM_DRAG
	init_phidot_table();
#endif
  }

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

  All.TimeLastRestartFile = CPUThisRun;

}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT
  double coulomb_log;
#endif
#endif
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

#ifdef SFR
  set_units_sfr();
#endif


#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT

#define cm (1.0/All.UnitLength_in_cm)
#define g  (1.0/All.UnitMass_in_g)
#define s  (1.0/All.UnitTime_in_s)
#define erg (1.0/All.UnitEnergy_in_cgs)
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)


  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  /* assuming full ionization */

  coulomb_log = 37.8;
  /* accordin1g to Sarazin's book */

  All.ConductionCoeff *=
    (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003 
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with 
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of  ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif

#endif /* CONDUCTION_CONSTANT */
#endif /* CONDUCTION */


#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
    (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs)))
     / ((NFW_Eps - 1) * (NFW_Eps - 1)) +
     (Rs * Rs * Rs *
      (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) +
       NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef SFR_PROMOTION
  sprintf(buf, "%spromotion_%d.txt", All.OutputDir, ThisTask);
  if(!(FdPromotion = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif


#ifdef SFR_METALS
  sprintf(buf, "%s%s", All.OutputDir, "energy_test.txt");
  if(!(FdMphase = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "SNenergy.txt");
  if(!(FdSNE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#if defined(SFR_SNI) || defined(SFR_SNII)
  sprintf(buf, "%s%s", All.OutputDir, "SN.txt");
  if(!(FdSN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif
#ifdef SFR_PROMOTION
  sprintf(buf, "%spromotion_%d.txt", All.OutputDir, ThisTask);
  if(!(FdPromotion = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif


#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif

#ifdef XXLINFO
  sprintf(buf, "%s%s", All.OutputDir, "xxl.txt");
  if(!(FdXXL = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdXXL, "nstep time ");
#ifdef MAGNETIC
	  fprintf(FdXXL, "<|B|> ");
#ifdef TRACEDIVB
	  fprintf(FdXXL, "max(divB) ");
#endif
#endif
#ifdef REDUCEVISC
	  fprintf(FdXXL, "<alpha> ");
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
    }
#endif

}

/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

#ifdef SFR_PROMOTION
  fclose(FdPromotion);
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
#ifdef SFR
  fclose(FdSfr);
#endif
#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif
#ifdef XXLINFO
  fclose(FdXXL);
#endif
#ifdef SFR_METALS
  fclose(FdMphase);
  fclose(FdSNE);
#if defined(SFR_SNI) || defined(SFR_SNII)
  fclose(FdSN);
#endif
#ifdef SFR_PROMOTION
  fclose(FdPromotion);
#endif
#endif
}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

  All.StarformationOn = 0;	/* defaults */


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = INT;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;
#ifdef COOLING
      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;
#endif

#ifdef SFR
      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;
#endif 

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = DOUBLE;

	/* Parameters for non standard cosmologies */
#ifdef DEDM_HUBBLE
        strcpy(tag[nt], "HubbleDeDmFile");
        addr[nt] = All.HubbleDeDmFile;
        id[nt++] = STRING;
#endif

#ifdef VDE
        strcpy(tag[nt], "HubbleVDEFile");
        addr[nt] = All.HubbleVDEFile;
        id[nt++] = STRING;
#endif

#ifdef DEDM_MASS
	strcpy(tag[nt], "VariableMassDeDmFile");
	addr[nt] = All.VariableMassDeDmFile;
	id[nt++] = STRING; 
#endif

#ifdef DEDM_DRAG
	strcpy(tag[nt], "PhiDeDmFile");
	addr[nt] = All.PhiDeDmFile;
	id[nt++] = STRING;
	strcpy(tag[nt], "HisDeDmFile");
	addr[nt] = All.HisDeDmFile;
	id[nt++] = STRING;
#endif

#ifdef DEDM_COUPLING
	strcpy(tag[nt], "BetaZero");
	addr[nt] = &All.BetaZero;
	id[nt++] = DOUBLE;
#endif



#ifdef NAVIERSTOKES
      strcpy(tag[nt], "NavierStokes_KinematicViscosity");
      addr[nt] = &All.NavierStokes_KinematicViscosity;
      id[nt++] = DOUBLE;
#endif


#ifdef CHEMISTRY
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = DOUBLE;
#endif


#ifdef CONDUCTION
      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionCoeff;
      id[nt++] = DOUBLE;
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BubbleTimeInterval");
      addr[nt] = &All.BubbleTimeInterval;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = DOUBLE;
#endif

#ifdef MULTI_BUBBLES
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ClusterMass200");
      addr[nt] = &All.ClusterMass200;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = DOUBLE;

#endif


#ifdef COSMIC_RAYS
      strcpy(tag[nt], "CR_SpectralIndex");
      addr[nt] = &All.CR_Alpha;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CR_SupernovaEfficiency");
      addr[nt] = &All.CR_SNEff;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CR_SupernovaSpectralIndex");
      addr[nt] = &All.CR_SNAlpha;
      id[nt++] = DOUBLE;

#ifdef CR_DIFFUSION
      strcpy(tag[nt], "CR_DiffusionCoeff");
      addr[nt] = &All.CR_DiffusionCoeff;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CR_DiffusionProton");
      addr[nt] = &All.CR_Diffusion_Proton;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CR_DiffusionGamma");
      addr[nt] = &All.CR_Diffusion_Gamma;
      id[nt++] = DOUBLE;
#endif /* CR_DIFFUSION */

#ifdef CR_SHOCK
      strcpy(tag[nt], "CR_ShockEfficiency");
      addr[nt] = &All.CR_ShockEfficiency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CR_ShockSpectralIndex");
      addr[nt] = &All.CR_ShockAlpha;
      id[nt++] = DOUBLE;
#endif /* CR_SHOCK */


#endif /* COSMIC_RAYS */


#ifdef BLACK_HOLES
      strcpy(tag[nt], "TimeBetBlackHoleSearch");
      addr[nt] = &All.TimeBetBlackHoleSearch;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = DOUBLE;


      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BlackHoleNgbFactor");
      addr[nt] = &All.BlackHoleNgbFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BlackHoleActiveTime");
      addr[nt] = &All.BlackHoleActiveTime;
      id[nt++] = DOUBLE;


#endif




#ifdef MOREPARAMS

#ifdef SFR_METALS
      strcpy(tag[nt], "FactorSFR");
      addr[nt] = &All.FactorSFR;
      id[nt++] = DOUBLE;


      strcpy(tag[nt], "TlifeSNII");
      addr[nt] = &All.TlifeSNII;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinTlifeSNI");
      addr[nt] = &All.MinTlifeSNI;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxTlifeSNI");
      addr[nt] = &All.MaxTlifeSNI;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "RateSNI");
      addr[nt] = &All.RateSNI;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "FactorSN_Phase");
      addr[nt] = &All.FactorSN_Phase;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "Tcrit_Phase");
      addr[nt] = &All.Tcrit_Phase;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DensFrac_Phase");
      addr[nt] = &All.DensFrac_Phase;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "FracEnergySN_Phase");
      addr[nt] = &All.FracEnergySN_Phase;
      id[nt++] = DOUBLE;


#endif


      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "WindFreeTravelLength");
      addr[nt] = &All.WindFreeTravelLength;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = DOUBLE;

      /*      strcpy(tag[nt], "FactorForSofterEQS");
         addr[nt] = &All.FactorForSofterEQS;
         id[nt++] = DOUBLE; */

#ifdef DARKENERGY
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = DOUBLE;
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = DOUBLE;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef REDUCEVISC
      strcpy(tag[nt], "ViscositySourceScaling");
      addr[nt] = &All.ViscSource0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ViscosityDecayLength");
      addr[nt] = &All.DecayLength;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ViscosityAlphaMin");
      addr[nt] = &All.AlphaMin;
      id[nt++] = DOUBLE;
#endif



#ifdef MAGNETIC
#ifdef BINISET
      strcpy(tag[nt], "BiniX");
      addr[nt] = &All.BiniX;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BiniY");
      addr[nt] = &All.BiniY;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BiniZ");
      addr[nt] = &All.BiniZ;
      id[nt++] = DOUBLE;
#endif

#ifdef BSMOOTH
      strcpy(tag[nt], "BSmoothInt");
      addr[nt] = &All.BSmoothInt;
      id[nt++] = INT;

      strcpy(tag[nt], "BSmoothFrac");
      addr[nt] = &All.BSmoothFrac;
      id[nt++] = DOUBLE;
#endif
#endif

#endif

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
#ifdef ALLOWEXTRAPARAMS
		      fprintf(stdout, "WARNING from file %s:   Tag '%s' ignored !\n", fname, buf1);
#else
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
#endif
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef SFR

#ifndef MOREPARAMS
  if(ThisTask == 0)
    {
      printf("Code was compiled with SFR, but not with MOREPARAMS.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#ifdef METALS
#ifndef SFR
  if(ThisTask == 0)
    {
      printf("Code was compiled with METALS, but not with SFR.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#ifdef SFR_METALS
#ifdef SFR_SNI
#ifndef SFR_ENRICH
  if(ThisTask == 0)
    {
      printf("Code was compiled with SNI, but not with ENRICH.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif
#ifdef SFR_SNII
#ifndef SFR_ENRICH
  if(ThisTask == 0)
    {
      printf("Code was compiled with SNII, but not with ENRICH.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif
#endif
#endif


#ifndef MOREPARAMS
#ifdef REDUCEVISC
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with REDUCEVISC, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DARKENERGY, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef TIMEDEPDE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MAGNETIC
#ifdef TRACEDIVB
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TRACEDIVB, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DBOUTPUT
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DBOUTPUT, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MAGFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGFORCE, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MAGFORCE
#ifdef ARTBPRES
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with ARTBPRES, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DIVBFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS


#ifdef COSMIC_RAYS
  if(ThisTask == 0)
    {
      printf("CR SN Efficiency: %g\n", All.CR_SNEff);
    }
#endif
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = log(TimeMax_new / All.TimeBegin) / All.Timebase_interval;
  else
    ti_end = (TimeMax_new - All.TimeBegin) / All.Timebase_interval;

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_endstep /= 2;
	}
    }

  All.TimeMax = TimeMax_new;
}
