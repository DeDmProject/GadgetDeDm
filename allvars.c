/*! \file allvars.c
 *  \brief creates global variables.
 *
 *  This file creates all the global variables that are declared in allvars.h
 *  
 *  To produce 'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#include "allvars.h"


int ThisTask;			/* the local processors  */
int NTask, PTask;		/* note: NTask = 2^PTask */

double CPUThisRun;

int NumForceUpdate;
int NumSphUpdate;
int TimeTreeRoot;
int RestartFlag;
int Flag_FullStep;
int TreeReconstructFlag;

int NumPart;			/* Note: this is the LOCAL processor value */
int N_gas;			/* Note: this is the LOCAL processor value */
long long Ntype[6];
int NtypeLocal[6];

#ifdef SFR
int Stars_converted;

#ifdef SFR_METALS
double TotalEnergy;
double DEnergy_spawned, DEnergy_converted;
double DEnergy_radiation, DEnergy_promotion;
double DEnergy_feedback, TotalReservoir;

#ifdef SFR_FEEDBACK
double ESN;
int nhot, ncold;
#endif
#endif
#endif

gsl_rng *random_generator;

#ifdef SFR_METALS
int Flag_phase;
int Flag_promotion;
#endif


struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;


double TimeOfLastTreeConstruction;

int *Ngblist;

peanokey *DomainKeyBuf;

peanokey *Key, *KeySorted;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int DomainMyStart, DomainMyLast;
int *DomainStartList, *DomainEndList;

double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;

int *DomainNodeIndex;

struct DomainNODE *DomainMoment;


FLOAT *DomainHmax;
FLOAT *DomainTreeNodeLen;



double RndTable[RNDTABLE];



/* variables for input/output ,  usually only used on process 0 
 */
char ParameterFile[100];
FILE *FdInfo, *FdEnergy, *FdTimings, *FdCPU;

#ifdef SFR
FILE *FdSfr;
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;
FILE *FdBlackHolesDetails;
#endif



#ifdef SFR_METALS
FILE *FdMphase;
FILE *FdSNE;

#if defined(SFR_SNI) || defined(SFR_SNII)
FILE *FdSN;
#endif
#ifdef SFR_PROMOTION
FILE *FdPromotion;
#endif
#endif


#ifdef FORCETEST
FILE *FdForceTest;
#endif

#ifdef XXLINFO
FILE *FdXXL;

#ifdef MAGNETIC
double MeanB;

#ifdef TRACEDIVB
double MaxDivB;
#endif
#endif
#ifdef REDUCEVISC
double MeanAlpha;
#endif
#endif



double DriftTable[DRIFT_TABLE_LENGTH], GravKickTable[DRIFT_TABLE_LENGTH], HydroKickTable[DRIFT_TABLE_LENGTH];

#ifdef DEDM_DRAG
double PhiDotIntegrationTable[PHIDOT_TABLE_LENGTH];
#endif

void *CommBuffer;		/* communication buffer, used at a number of places */

char *Exportflag;


/* this structure contains data which is the SAME for all 
 * tasks (mostly code parameters read from the parameter file). 
 * Holding this data in a structure is convenient for writing/reading
 * the restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce them
 * into this structure.
 */
struct global_data_all_processes All;



/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P, *DomainPartBuf;



/* the following struture holds data that is stored for each SPH particle
 * in addition to the collisionless variables.
 */
struct sph_particle_data *SphP, *DomainSphBuf;




/* global state of system 
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;






/* Various structure for communication during the gravity 
 * computation.
 */
struct gravdata_in *GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;

struct gravdata_index *GravDataIndexTable;

#ifdef SFR_METALS
struct metaldata_in *MetalDataIn, *MetalDataGet;
#endif

/* Various structure for communication during the density
 * computation.
 */
struct densdata_in *DensDataIn, *DensDataGet;

struct densdata_out *DensDataResult, *DensDataPartialResult;


#ifdef BLACK_HOLES
struct blackholedata_in *BlackholeDataIn, *BlackholeDataGet;
struct blackholedata_out *BlackholeDataResult, *BlackholeDataPartialResult;
#endif


#ifdef FOF
struct fofdata_in *FoFDataIn, *FoFDataGet;
struct fofdata_out *FoFDataResult, *FoFDataPartialResult;
#endif




/* Various structures for communication during the 
 * computation of hydrodynamical forces.
 */
struct hydrodata_in *HydroDataIn, *HydroDataGet;

struct hydrodata_out *HydroDataResult, *HydroDataPartialResult;


#ifdef SFR_PROMOTION
struct hotngbs_in *HotNgbsIn, *HotNgbsGet;
struct hotngbs_out *HotNgbsResult, *HotNgbsPartialResult;
#endif

#ifdef MHM
struct kindata_in *KinDataIn, *KinDataGet;
#endif


/* Header for the standard file format.
 */
struct io_header header;


char Tab_IO_Labels[IO_NBLOCKS][4];


/*******************
 ** Variables for Tree
 ********************
 */


struct NODE *Nodes, *Nodes_base;

struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/* maximum allowed number of internal nodes */
int Numnodestree;		/* number of (internal) nodes in each tree */


int *Nextnode;			/* gives next node in tree walk  (nodes array) */
int *Father;			/* gives parent node in tree (Prenodes array) */


#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif



#ifdef CHEMISTRY
/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T], k10a[N_T],
  k11a[N_T];
double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T], k20a[N_T],
  k21a[N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif

#endif /* chemisitry */
