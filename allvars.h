/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "tags.h"
#include "assert.h"

#ifdef CHEMISTRY
#include "libbaryons/chemistry.h"
#endif

#define  GADGETVERSION   "2.0"   /*!< code version string */

#define  GENERATIONS     2       /*!< Number of star particles that may be created per gas particle
                                  */
 

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                  *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                  *   to 2^29
                                  */

#ifndef LARGE_SIMU
#define MAXTOPNODES    200000
#else
#define MAXTOPNODES    80000
#endif

typedef  long long  peanokey;

#define  BITS_PER_DIMENSION 18	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))




#define  GAMMA         (5.0/3)    /*!< adiabatic index of simulated gas */
#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76   /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  METAL_YIELD       0.02   /*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 128

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  PI		3.14159

#ifdef CHEMISTRY
#define  T_CMB0      2.728	/* present-day CMB temperature */
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25		
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5		
#endif

#ifndef HPM_SMTH
#define HPM_SMTH ASMTH
#endif

#define COND_TIMESTEP_PARAMETER 0.25

#define MAX_NGB  20000		/*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST 999	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000  /*!< length of the lookup table used to hold the drift and kick factors */ 

#ifdef DEDM_DRAG
#define PHIDOT_TABLE_LENGTH 1000
#endif

#define MAXITER 150

#define LINKLENGTH 0.2
#define GROUP_MIN_LEN 32

#define MINRESTFAC 0.05

#ifdef DOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif

#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                           /*!< Coefficient for kernel normalization. */
#endif


#if defined(SFR_METALS) || defined (BLACK_HOLES)
#define PPP P
#ifdef SFR_FEEDBACK
#define EgySNcgs 1e52
#endif
extern int Flag_promotion;   
extern int Flag_phase; 
#else
#define PPP SphP
#endif


#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)



extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;               /*!< number of processors */
extern int PTask;	        /*!< note: NTask = 2^PTask */

extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
extern int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

extern int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                  initial conditions, 1 is resuming a run from a set of restart files, while 2
                                  marks a restart from a snapshot file. */

extern char *Exportflag;


extern int Flag_FullStep;       /*!< Flag used to signal that the current step involves all particles */


extern int TreeReconstructFlag;  
                             

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern long long Ntype[6];      /*!< total number of particles of each type */
extern int NtypeLocal[6];       /*!< local number of particles of each type */

extern gsl_rng *random_generator; /*!< the random number generator used */


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */

#ifdef SFR_METALS
  extern double Ucrit;
  extern double TotalEnergy;
  extern double DEnergy_spawned, DEnergy_converted;
  extern double DEnergy_radiation, DEnergy_promotion;
  extern double DEnergy_feedback, TotalReservoir;
#ifdef SFR_FEEDBACK
  extern double ESN;
  extern int nhot, ncold;
#endif  
#endif  

#endif


extern double TimeOfLastTreeConstruction; /*!< holds what it says */

extern int *Ngblist;            /*!< Buffer to hold indices of neighbours retrieved by the neighbour search
                                  routines */




extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int DomainMyStart, DomainMyLast;
extern int *DomainStartList, *DomainEndList;



extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern FLOAT *DomainTreeNodeLen;
extern FLOAT *DomainHmax;

extern struct DomainNODE
{
  FLOAT s[3];
  FLOAT vs[3];
  FLOAT mass;

#ifdef DEDM_TREE
FLOAT s_dm[3];
FLOAT vs_dm[3];
FLOAT mass_dm;
#endif

#ifdef UNEQUALSOFTENINGS
  int   bitflags;
#endif
}
*DomainMoment;


extern peanokey *Key, *KeySorted;

extern struct topnode_data
{
  int Daughter;
  peanokey Size;
  peanokey StartKey;
  long long Count;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;




extern double RndTable[RNDTABLE];


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];  /*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,   /*!< file handle for info.txt log-file. */
  *FdEnergy,           /*!< file handle for energy.txt log-file. */
  *FdTimings,          /*!< file handle for timings.txt log-file. */
  *FdCPU;              /*!< file handle for cpu.txt log-file. */

#ifdef SFR
extern FILE *FdSfr;    /*!< file handle for sfr.txt log-file. */
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;    /*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif



#ifdef SFR_METALS
extern FILE *FdMphase; 
extern FILE *FdSNE; 
#if defined(SFR_SNI) || defined(SFR_SNII)
extern FILE *FdSN;
#endif    
#ifdef SFR_PROMOTION
extern FILE *FdPromotion;
#endif
#endif

#ifdef FORCETEST
extern FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;  /*!< file handle for xxl.txt log-file. */
#ifdef MAGNETIC
extern double MeanB;
#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef REDUCEVISC
extern double MeanAlpha;
#endif 
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];  

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef DEDM_DRAG
extern double PhiDotIntegrationTable[PHIDOT_TABLE_LENGTH];
#endif

extern void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                             parallel tree-force computation, and the SPH routines. */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				     processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				     processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel; /*!< maximum number of files that may be written simultaneously when
                                   writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force
				   algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density
                                  computation */
#ifdef FOF
  int BunchSizeFoF;
#endif

#ifdef SFR_METALS
  int BunchSizeMetal;
#ifdef SFR_PROMOTION
  int BunchSizeHotNgbs;
#endif
#endif  

#ifdef BLACK_HOLES
  int BunchSizeBlackhole;
#endif

#ifdef MHM
  int BunchSizeKinetic;
#endif

  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH
                                  hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain
                                  decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  int DesNumNgb;                /*!< Desired number of SPH neighbours */
  int MaxNumNgbDeviation;       /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */



  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */

  /* system of units  */

  double UnitTime_in_s,   	/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,             	/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,   /*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,           /*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,        /*!< factor to convert internal pressure unit to cgs units (little 'h' still
                                     around!) */
    UnitDensity_in_cgs,         /*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,     /*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,          /*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,      /*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,    /*!< If set to zero in the parameterfile, the internal value of the
                                  gravitational constant is set to the Newtonian value based on the system of
                                  units specified. Otherwise the value provided is taken as internal gravity
                                  constant G. */
    G;                          /*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;  /*!< Hubble-constant in internal units */
  double Omega0,  /*!< matter density in units of the critical density (at z=0)*/
    OmegaLambda,  /*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,  /*!< baryon density in units of the critical density (at z=0)*/
    HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
		   * physical values for cooling physics
                   */

  double BoxSize; /*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
                                  criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike
                                  gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */
  int CoolingOn;                /*!< flags that cooling is enabled */
  int StarformationOn;          /*!< flags that star formation is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;     /*!< number of snapshot that is written next */
  double TimeBetSnapshot,    /*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,     /*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,   /*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,     /*!< cpu-time when last restart-file was written */
    TimeBetStatistics,       /*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;      /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;      /*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,  /*!< current time of the simulation */
    TimeBegin,  /*!< time of initial conditions of the simulation */
    TimeStep,   /*!< difference between current times of previous and current timestep */
    TimeMax;	/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval; /*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;           /*!< current time on integer timeline */ 
  int Ti_nextoutput;        /*!< next output time on integer timeline */

#ifdef FLEXSTEPS
   int PresentMinStep;
#endif


#ifdef PMGRID
  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

#ifdef CHEMISTRY
  double Epsilon;
#endif 


  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_TreeConstruction;
  double CPU_TreeWalk;
  double CPU_Gravity;
  double CPU_Potential;
  double CPU_Domain;
  double CPU_Snapshot;
  double CPU_Total;
  double CPU_CommSum;
  double CPU_Imbalance;
  double CPU_HydCompWalk;
  double CPU_HydCommSumm;
  double CPU_HydImbalance;
  double CPU_Hydro;
  double CPU_EnsureNgb;
  double CPU_Predict;
  double CPU_TimeLine;
  double CPU_PM;
  double CPU_Peano;
#ifdef COOLING
  double CPU_SfrCool;
#endif

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                  timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                  timestep determined by the timestep criteria falls below this limit. */ 
         MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional, /*!< minimum allowed SPH smoothing length in units of SPH gravitational
                                  softening length */
    MinGasHsml;                /*!< minimum allowed SPH smoothing length */


  double SofteningGas,    /*!< for type 0 */ 
    SofteningHalo,        /*!< for type 1 */ 
    SofteningDisk,        /*!< for type 2 */ 
    SofteningBulge,       /*!< for type 3 */ 
    SofteningStars,       /*!< for type 4 */ 
    SofteningBndry;       /*!< for type 5 */ 

  double SofteningGasMaxPhys,   /*!< for type 0 */ 
    SofteningHaloMaxPhys,       /*!< for type 1 */ 
    SofteningDiskMaxPhys,       /*!< for type 2 */ 
    SofteningBulgeMaxPhys,      /*!< for type 3 */ 
    SofteningStarsMaxPhys,      /*!< for type 4 */ 
    SofteningBndryMaxPhys;      /*!< for type 5 */ 

  double SofteningTable[6];  /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];  /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], 
    TimingsFile[100], 
    RestartFile[100], 
    ResubmitCommand[100], 
    OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];

  int OutputListLength; /*!< number of times stored in table of desired output times */



#ifdef MOREPARAMS		/* star formation and feedback sector */
  double OrigGasMass;
  double EgySpecCold;
  double EgySpecSN;
  double OverDensThresh;
  double PhysDensThresh;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double CritOverDensity;
  double CritPhysDensity;
  double FactorSN;
  double FactorEVP;
  double MaxSfrTimescale;
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelLength;
  double WindFreeTravelDensFac;
  double FactorForSofterEQS;

#ifdef SFR_METALS
  double FactorSFR;
  double MinTlifeSNI;
  double MaxTlifeSNI;
  double TlifeSNII;
  double RateSNI;
  double FactorSN_Phase;
  double Tcrit_Phase;
  double DensFrac_Phase;
  double FracEnergySN_Phase;
#ifdef SFR_DECOUPLING
  double DensityTailThreshold; 
#endif
#endif

#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#ifdef REDUCEVISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution*/
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#ifdef CONDUCTION
  double ConductionCoeff;         /*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;  /*!< Factor to get electron mean free path */
#endif
#endif

#ifdef BINISET                 
  double BiniX,BiniY,BiniZ;       /*!< Initial values for B */
#endif

#ifdef BSMOOTH
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#endif

#ifdef BLACK_HOLES
  double TimeNextBlackHoleCheck;
  double TimeBetBlackHoleSearch;
  double BlackHoleAccretionFactor;  /*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;   /*!< Fraction of the black luminosity feed into thermal feedback */
    double SeedBlackHoleMass;         /*!< Seed black hole mass */
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;        /*!< Factor by which the normal SPH neighbour should be increased/decreased */
  double BlackHoleActiveTime;
  double BlackHoleEddingtonFactor; /*! Factor above Eddington */
#ifdef MODIFIEDBONDI
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;
#endif
#endif

#ifdef COSMIC_RAYS
  double CR_Alpha;               /*!< Cosmic ray spectral index [2..3]*/
  double CR_SNEff;               /*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;             /*!< SN injection spectral index [2..3] */
  int bDebugFlag;                /*!< enables debug outputs after triggered */

#ifdef CR_DIFFUSION
  double CR_Diffusion_Gamma;     /*!< CR diffusion coefficient d_gamma */
  double CR_Diffusion_Proton;    /*!< CR diffusion coefficient d_p     */

  double CR_DiffusionCoeff;      /*!< (temporary) fixed value for CR diffusivity */
#endif /* CR_DIFFUSION */

#if defined(CR_SHOCK)
  double CR_ShockAlpha;          /*!< spectral index to be used in shock injection */
  double CR_ShockEfficiency;     /*!< energy fraction of shock energy fed into CR */
#endif /* CR_SHOCK */

#endif /* COSMIC_RAYS */


#ifdef HPM
  double HPM_entr0, HPM_entr1;
  double HPM_ne0, HPM_ne1;
  double HPM_rho0, HPM_rho1;
  double HPM_P0, HPM_P1, HPM_alpha;
#endif


#ifdef BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
#ifdef FOF
  int    BiggestGroupLen;
  float  BiggestGroupCM[3];
  double BiggestGroupMass;   
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double ClusterMass200;
  double massDMpart;
#endif

#ifdef NAVIERSTOKES
  double NavierStokes_KinematicViscosity;
#endif

#ifdef DEDM_HUBBLE
  char HubbleDeDmFile[100];
#endif

#ifdef VDE
  char HubbleVDEFile[100];
#endif

#ifdef DEDM_MASS
  char VariableMassDeDmFile[100];
#endif

#ifdef DEDM_DRAG
  char PhiDeDmFile[100];
  char HisDeDmFile[100];
#endif

#ifdef DEDM_INFO
	FILE *outDeDmFile;
#endif

#ifdef DEDM_COUPLING
  double BetaZero;
#endif

}
All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  FLOAT Pos[3];			/*!< particle position at its current time */
  FLOAT Mass;			/*!< particle mass */
  FLOAT Vel[3];			/*!< particle velocity at its current time */
  FLOAT GravAccel[3];		/*!< particle acceleration due to gravity */
#ifdef PMGRID

  FLOAT GravPM[3];		/*!< particle acceleration due to long-range PM gravity force */
#if defined(DEDM_PM) || defined(DEDM_PMb)
  FLOAT DeDmPM[3];		/*!< particle acceleration due to long-range PM quintessence mediated force */
#endif
#endif

#ifdef DEDM_TREE
  FLOAT DeDmAccel[3]; 		/*!< particle acceleration due to quintessence mediated force */
#endif


#ifdef FORCETEST
  FLOAT GravAccelDirect[3];	/* particle acceleration */
#endif
  FLOAT Potential;		/*!< gravitational potential */
  FLOAT OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                  criterion */

#if defined(EVALPOTENTIAL) && defined(PMGRID)
  FLOAT PM_Potential;
#endif

#ifdef STELLARAGE
  FLOAT StellarAge;		/*!< formation time of star particle */
#endif

#ifdef METALS
#ifdef SFR_METALS
  FLOAT Zm[12];
  FLOAT ZmReservoir[12];
#ifdef SFR_FEEDBACK
  FLOAT EnergySN;
  FLOAT EnergySNCold;    
#endif    
#else
  FLOAT Metallicity;		/*!< metallicity of gas or star particle */
#endif
#endif  /* closes METALS */

#if defined(SFR_METALS) || defined (BLACK_HOLES)
  FLOAT Hsml;	
  FLOAT Left,                   /*!< lower bound in iterative smoothing length search */  
    Right;                      /*!< upper bound in iterative smoothing length search */ 
#ifndef  NOFIXEDMASSINKERNEL
  FLOAT NumNgb;
#else
  int NumNgb;			/*!< (effective) number of SPH neighbours */
#endif
#endif

#ifndef LONGIDS
  unsigned int ID;			/*!< particle identifier */
#else
  unsigned long long ID;
#endif
  int Type;		/*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  int Ti_endstep;          /*!< marks start of current timestep of particle on integer timeline */ 
  int Ti_begstep;          /*!< marks end of current timestep of particle on integer timeline */
  float GravCost;		/*!< weight factor used for balancing the work-load */

#ifdef PSEUDOSYMMETRIC
  float AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
                                  prediction of the change of acceleration expected in the future, thereby
                                  allowing to guess whether a decrease/increase of the timestep should occur
                                  in the timestep that is started. */
#endif

#ifdef BLACK_HOLES
  FLOAT BH_Mass;
  FLOAT BH_Mdot;
  FLOAT BH_MdotEddington;   /* in units of the Eddington accretion rate */
  FLOAT BH_Density;
  FLOAT BH_Entropy;
  FLOAT BH_SurroundingGasVel[3];
#ifdef REPOSITION_ON_POTMIN
  FLOAT BH_MinPotPos[3];
  FLOAT BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  FLOAT ActiveTime;
  FLOAT ActiveEnergy;
#endif
#endif
}
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */
  FLOAT Density;		/*!< current baryonic mass density of particle */

#if !defined(SFR_METALS) && !defined(BLACK_HOLES)
  FLOAT Hsml;			/*!< current smoothing length */
  FLOAT Left,                   /*!< lower bound in iterative smoothing length search */  
    Right;                      /*!< upper bound in iterative smoothing length search */ 
#ifndef  NOFIXEDMASSINKERNEL
  FLOAT NumNgb;
#else
  int NumNgb;                   /*!< (effective) number of SPH neighbours */
#endif
#endif
  FLOAT Pressure;		/*!< current pressure */
  FLOAT DtEntropy;              /*!< rate of change of entropy */
  FLOAT HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
  FLOAT VelPred[3];		/*!< predicted SPH particle velocity at the current time */

  union
  {
#ifdef NAVIERSTOKES
    FLOAT DV[3][3];
#endif
    struct
    {
      FLOAT DivVel;			/*!< local velocity divergence */
      FLOAT CurlVel;	 	        /*!< local velocity curl */
#ifdef NAVIERSTOKES
      FLOAT StressDiag[3];
      FLOAT StressOffDiag[3];
#else
      FLOAT Rot[3];		        /*!< local velocity curl */
#endif
    }
    s;
  }
  u;

#ifndef NOGRADHSML
  FLOAT DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative
                                  entropy formulation of SPH */
#endif
  FLOAT MaxSignalVel;           /*!< maximum signal velocity */
#ifdef COOLING
  FLOAT Ne;			/*!< electron fraction, expressed as local electron number density normalized
                                  to the hydrogen number density. Gives indirectly ionization state and mean
                                  molecular weight. */
#ifdef SFR
  FLOAT Sfr;
#endif
#endif
#ifdef WINDS
  FLOAT DelayTime;              /*!< remaining maximum decoupling time of wind particle */
#endif

#ifdef MAGNETIC
  FLOAT B[3], BPred[3];
  FLOAT DtB[3];
#ifdef TRACEDIVB
  FLOAT divB;
#endif
#ifdef BSMOOTH
  FLOAT BSmooth[3];
  FLOAT DensityNorm;
#endif
#endif
#ifdef REDUCEVISC
  FLOAT alpha, Dtalpha;
#endif
#ifdef CONDUCTION
  FLOAT CondEnergyChange;
  FLOAT SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#ifdef OUTPUTCOOLRATE
  FLOAT CondRate;
#endif
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
  FLOAT Injected_BH_Energy;
#endif


#if defined (SFR_METALS) && defined (SFR_DECOUPLING)
  FLOAT DensityOld;     
#endif
#ifdef SFR_PROMOTION
  FLOAT DensityAvg;
  FLOAT EntropyAvg;
  FLOAT HotHsml;
  int   HotNgbNum;
  FLOAT DensPromotion;
  FLOAT TempPromotion;
#endif
#ifdef MHM
  FLOAT FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  FLOAT CR_C0;                  /*!< Cosmic ray amplitude adiabatic invariable */
  FLOAT CR_q0;                  /*!< Cosmic ray cutoff adiabatic invariable */

  FLOAT CR_Rho0;                /*!< Zero-Density of CR Interpolation values */
  FLOAT CR_P0;                  /*!< CR Pressure at Rho0 */
  FLOAT CR_Gamma0;              /*!< CR Adiabatic Index at Rho0 */
  FLOAT CR_E0;                  /*!< Specific Energy at Rho0 */
  FLOAT CR_n0;                  /*!< baryon fraction in cosmic rays */

  FLOAT CR_DeltaE;              /*!< Specific Energy growth during timestep */
  FLOAT CR_DeltaN;              /*!< baryon fraction growth during timestep */
#ifdef CR_DIFFUSION
  FLOAT CR_Diff_E;              /*!< CR diffusion term for energy */
  FLOAT CR_Diff_N;              /*!< CR diffusion term for baryon fraction */
#endif /* CR_DIFFUSION */

#endif /* COSMIC_RAYS */

#ifdef MACHNUM
  FLOAT MachNumber;              /*!< Mach number finder */
  FLOAT DissipatedEnergy;        /*!< Total dissipated energy */
#endif /* Mach number estimate */


#ifdef CHEMISTRY
  FLOAT elec;            
  FLOAT HI;
  FLOAT HII;

  FLOAT HeI;
  FLOAT HeII;
  FLOAT HeIII;

  FLOAT H2I;
  FLOAT H2II;

  FLOAT HM;

  FLOAT Gamma;
  FLOAT t_elec, t_cool;
#endif

}
 *SphP,                        	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system 
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6], 
    EnergyTotComp[6], 
    MomentumComp[6][4], 
    AngMomentumComp[6][4], 
    CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
  }
  u;

#ifdef DEDM_TREE
union 
	{
  FLOAT Pos_dm[3];
  FLOAT Acc_dm[3];
	} 
	u_dm;
#endif

  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w; 

#ifdef DEDM_TREE
union	// TODO [union not really used...]
	{	
	FLOAT OldAcc_dm;
	int Ninteractions_dm;
	} 
	w_dm;
#endif

#if (UNEQUALSOFTENINGS) || defined(DEDM_TREE) || defined(EVALPOTENTIAL) || defined(OUTPUTPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union 
  {
    FLOAT Potential;
    int Type;
  }
  v;
#endif
}
*GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;


extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;



/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
extern struct densdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
#ifdef WINDS
  FLOAT DelayTime;
#endif
  int Index;
  int Task;
#ifdef SFR_DECOUPLING
  FLOAT DensityOld;
  FLOAT Entropy;
#endif
}
 *DensDataIn, *DensDataGet;

#ifdef SFR_METALS
extern struct metaldata_in
{
  float Pos[3];
  float ZmReservoir[12];
  float Hsml;
#ifdef SFR_FEEDBACK
  float EnergySN;
  FLOAT EnergySNCold;    
#endif  
#ifndef  NOFIXEDMASSINKERNEL
  FLOAT NumNgb;
#else	
  int NumNgb;
#endif
  int Index;
  int Task;
} *MetalDataIn, *MetalDataGet;
#endif


#ifdef MHM
extern struct kindata_in
{
  FLOAT Pos[3];
  FLOAT Hsml;
  FLOAT Density;
  FLOAT Energy;
  int Index;
  int Task;
}
 *KinDataIn, *KinDataGet;
#endif


#ifdef BLACK_HOLES
extern struct blackholedata_in
{
  FLOAT Pos[3];
  FLOAT Density;
  FLOAT Mdot;
  FLOAT Dt;
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT BH_Mass;
  FLOAT Vel[3];
#ifdef BH_KINETICFEEDBACK
  FLOAT ActiveTime;
  FLOAT ActiveEnergy;
#endif
  int ID;
  int Index;
  int Task;
}
 *BlackholeDataIn, *BlackholeDataGet;

extern struct blackholedata_out
{
  FLOAT Mass;
  FLOAT BH_Mass;
  FLOAT AccretedMomentum[3];
#ifdef REPOSITION_ON_POTMIN
  FLOAT BH_MinPotPos[3];
  FLOAT BH_MinPot;
#endif
}
 *BlackholeDataResult, *BlackholeDataPartialResult;
#endif






/*! Structure for communication during the density computation. Holds data that is received from other
 * processors.
 */
extern struct densdata_out
{
  FLOAT Rho;
#ifndef NAVIERSTOKES
  FLOAT Div, Rot[3];
#else
  FLOAT DV[3][3];
#endif
#ifndef NOGRADHSML
  FLOAT DhsmlDensity;
#endif
#ifndef NOFIXEDMASSINKERNEL
  FLOAT Ngb;
#else
  int Ngb;
#endif
#if defined(MAGNETIC) && defined(BSMOOTH)
  FLOAT BSmooth[3];
  FLOAT DensityNorm;
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
  FLOAT SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#endif

#ifdef BLACK_HOLES
  FLOAT GasVel[3];
#endif
}
 *DensDataResult, *DensDataPartialResult;


#ifdef SFR_PROMOTION
extern struct hotngbs_in
{
  FLOAT Pos[3];
  FLOAT HotHsml;
  FLOAT Entropy;
  int Index;
  int Task;
}
 *HotNgbsIn, *HotNgbsGet;

extern struct hotngbs_out
{
  FLOAT DensitySum;
  FLOAT EntropySum;
  int   HotNgbNum;
}
*HotNgbsResult, *HotNgbsPartialResult;
#endif


#ifdef FOF
extern struct fofdata_in
{
  FLOAT Pos[3];
  FLOAT Hsml;
  int MinID;
  int MinIDTask;
  int Index;
  int Task;
}
 *FoFDataIn, *FoFDataGet;

extern struct fofdata_out
{
  FLOAT Distance;
  int   MinID;
  int   MinIDTask;
}
*FoFDataResult, *FoFDataPartialResult;

#endif











/* Various structures for communication during the computation of hydrodynamical forces.
 */
extern struct hydrodata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  FLOAT F1;
#ifndef NOGRADHSML
  FLOAT DhsmlDensityFactor;
#endif
#ifdef MAGNETIC
  FLOAT BPred[3];
#endif
#ifdef REDUCEVISC
  FLOAT alpha;
#endif
  int   Timestep;
  int   Task;
  int   Index;

#ifdef CONDUCTION
  FLOAT SmoothedEntr;
  FLOAT Entropy;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#endif

#ifdef SFR_DECOUPLING
  FLOAT DensityOld;
  FLOAT Entropy;
#endif

#ifdef PARTICLE_DEBUG
#ifndef LONGIDS
  unsigned int ID;			/*!< particle identifier */
#else
  unsigned long long ID;
#endif
#endif 

#ifdef CR_DIFFUSION
  FLOAT CR_Diff_E;       /*!< diffusion term for cosmic ray energy */
  FLOAT CR_Diff_N;       /*!< diffusion term for cosmic ray baryon fraction */
#endif

#ifdef NAVIERSTOKES
  FLOAT stressoffdiag[3];
  FLOAT stressdiag[3];
  FLOAT shear_viscosity;
#endif

}
 *HydroDataIn, *HydroDataGet;

extern struct hydrodata_out
{
  FLOAT Acc[3];
  FLOAT DtEntropy;
  FLOAT MaxSignalVel;
#ifdef MAGNETIC
  FLOAT DtB[3];
#ifdef TRACEDIVB
  FLOAT divB;
#endif
#endif
#ifdef CONDUCTION
  FLOAT CondEnergyChange;
#ifdef OUTPUTCOOLRATE
  FLOAT CondRate;
#endif
#endif

#if  defined(CR_SHOCK) || defined(CR_DIFFUSION) 
  FLOAT CR_EnergyChange;
  FLOAT CR_BaryonFractionChange;
#endif
}
 *HydroDataResult, *HydroDataPartialResult;



/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];           /*!< number of particles of each type in this file */
  double mass[6];              /*!< mass of particles of each type. If 0, then the masses are explicitly
                                 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                 /*!< time of snapshot file */
  double redshift;             /*!< redshift of snapshot file */
  int flag_sfr;           /*!< flags whether the simulation was including star formation */
  int flag_feedback;      /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];      /*!< total number of particles of each type in this snapshot. This can be
                                        different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;       /*!< flags whether cooling was included  */
  int num_files;          /*!< number of files in multi-file snapshot */
  double BoxSize;              /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;               /*!< matter density in units of critical density */
  double OmegaLambda;          /*!< cosmological constant parameter */
  double HubbleParam;          /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;    /*!< flags whether the file contains formation times of star particles */
  int flag_metals;        /*!< flags whether the file contains metallicity values for gas and star
                                 particles */
  unsigned int npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	       /*!< fills to 256 Bytes */
}
header;  /*!< holds header for snapshot files */




#ifndef CHEMISTRY
#define IO_NBLOCKS 32   /*!< total number of defined information blocks for snapshot files.
                         Must be equal to the number of entries in "enum iofields" */
#else
#define IO_NBLOCKS 41
#endif


enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_CR_C0,
  IO_CR_Q0,
  IO_RHO,
  IO_NE,
  IO_NH,
#ifdef CHEMISTRY
  IO_ELECT,
  IO_HI,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
#endif
  IO_HSML,
  IO_SFR,
  IO_AGE,
  IO_Z,
  IO_POT,
  IO_ACCEL,
  IO_DTENTR,
  IO_TSTP,
  IO_BFLD,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_BSMTH,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_BHMASS,
  IO_BHMDOT,
  IO_MACH,
  IO_DISSENERGY
};


extern char Tab_IO_Labels[IO_NBLOCKS][4];







/*
 * Variables for Tree
 * ------------------
 */

extern struct NODE
{
  FLOAT len;			/*!< sidelength of treenode */
  FLOAT center[3];		/*!< geometrical center of node */

  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];               /*!< center of mass of node */
      FLOAT mass;      /*!< mass of node */

#ifdef DEDM_TREE
     FLOAT s_dm[3];		/*!< center of mass of dm particles of node */
     FLOAT mass_dm;		/*!< total mass of dm particles of node */
#endif

      int bitflags;        /*!< flags certain node properties */
      int sibling;         /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;        /*!< this gives the next node in case the current node needs to be opened */
      int father;          /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 					gives the first allocated node */


extern struct extNODE
{
  FLOAT hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  FLOAT vs[3];			/*!< center-of-mass velocity */

#ifdef DEDM_TREE
  FLOAT vs_dm[3];		/*!< center-of-mass velocity for DM particles */
#endif

}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;	/*!< gives next node in tree walk  (nodes array) */
extern int *Father;	/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif


#ifdef CHEMISTRY
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
extern double  T[N_T],J0_nu[N_nu],J_nu[N_nu],nu[N_nu];
extern double  k1a[N_T],k2a[N_T],k3a[N_T],k4a[N_T],k5a[N_T],k6a[N_T],k7a[N_T],k8a[N_T],k9a[N_T],k10a[N_T],k11a[N_T];
extern double  k12a[N_T],k13a[N_T],k14a[N_T],k15a[N_T],k16a[N_T],k17a[N_T],k18a[N_T],k19a[N_T],k20a[N_T],k21a[N_T];
extern double  ciHIa[N_T],ciHeIa[N_T],ciHeIIa[N_T],ciHeISa[N_T],reHIIa[N_T],brema[N_T];
extern double  ceHIa[N_T],ceHeIa[N_T],ceHeIIa[N_T],reHeII1a[N_T],reHeII2a[N_T],reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
extern double sigma24[N_nu],sigma25[N_nu],sigma26[N_nu],sigma27[N_nu],sigma28[N_nu],sigma29[N_nu],sigma30[N_nu],sigma31[N_nu];
#endif
#endif

#endif
