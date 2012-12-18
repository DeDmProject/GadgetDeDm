#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef VDE
#include "libvde/interpolate.h"
#endif

#ifdef DEDM_HUBBLE
#include "libdedm/interpolate.h"
#endif

#ifdef DARKENERGY
#include "libdarkenergy/darkenergy.h"
#endif


/*! \file ngb.c
 *  \brief neighbour search by means of the tree
 *
 * This file contains routines for neighbour finding.  We use the
 * gravity-tree and a range-searching technique to find neighbours.
 */

#ifdef PERIODIC
static double boxHalf, boxSize;

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

/*! this macro maps a coordinate difference to the nearest periodic
 * image
 */

#ifdef LONG_X
#define NGB_PERIODIC_LONG_X(x) (((x)>boxHalf*LONG_X)?((x)-boxSize*LONG_X):(((x)<-boxHalf*LONG_X)?((x)+boxSize*LONG_X):(x)))
#else
#define NGB_PERIODIC_LONG_X(x) (((x)>boxHalf)?((x)-boxSize):(((x)<-boxHalf)?((x)+boxSize):(x)))
#endif

#ifdef LONG_Y
#define NGB_PERIODIC_LONG_Y(x) (((x)>boxHalf*LONG_Y)?((x)-boxSize*LONG_Y):(((x)<-boxHalf*LONG_Y)?((x)+boxSize*LONG_Y):(x)))
#else
#define NGB_PERIODIC_LONG_Y(x) (((x)>boxHalf)?((x)-boxSize):(((x)<-boxHalf)?((x)+boxSize):(x)))
#endif

#ifdef LONG_Z
#define NGB_PERIODIC_LONG_Z(x) (((x)>boxHalf*LONG_Z)?((x)-boxSize*LONG_Z):(((x)<-boxHalf*LONG_Z)?((x)+boxSize*LONG_Z):(x)))
#else
#define NGB_PERIODIC_LONG_Z(x) (((x)>boxHalf)?((x)-boxSize):(((x)<-boxHalf)?((x)+boxSize):(x)))
#endif



/*! This routine finds all neighbours `j' that can interact with the
 * particle `i' in the communication buffer.
 *
 * Note that an interaction can take place if 
 * \f$ r_{ij} < h_i \f$  OR if  \f$ r_{ij} < h_j \f$. 
 * 
 * In the range-search this is taken into account, i.e. it is
 * guaranteed that all particles are found that fulfil this condition,
 * including the (more difficult) second part of it. For this purpose,
 * each node knows the maximum h occuring among the particles it
 * represents.
 */

#ifdef SFR_DECOUPLING
int ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode, FLOAT densityold, FLOAT entropy,
		       FLOAT * vel)
#else
int ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
#endif
{
  int k, no, p, numngb;
  FLOAT hdiff;
  FLOAT searchmin[3], searchmax[3];
  struct NODE *this;


#ifdef SFR_DECOUPLING
  double fac_mu, hubble_a, hubble_a2;
  double dx, dy, dz, r2, dvx, dvy, dvz, vdotr, vdotr2, soundspeed_i, soundspeed_j, c_ij, h_ij, mu_ij;

  if(All.ComovingIntegrationOn)
    {
      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
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
    }
  else
    {
      fac_mu = hubble_a = hubble_a2 = 1;
    }

#endif



  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  hdiff = PPP[p].Hsml - hsml;
	  if(hdiff < 0)
	    hdiff = 0;

#ifdef SFR_DECOUPLING
	  if(densityold > 0)	/* note: stars should not get in here... */
	    {
	      dx = searchcenter[0] - P[p].Pos[0];
	      dy = searchcenter[1] - P[p].Pos[1];
	      dz = searchcenter[2] - P[p].Pos[2];
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
	      dvx = vel[0] - SphP[p].VelPred[0];
	      dvy = vel[1] - SphP[p].VelPred[1];
	      dvz = vel[2] - SphP[p].VelPred[2];
	      vdotr = dx * dvx + dy * dvy + dz * dvz;
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(All.ComovingIntegrationOn)
		vdotr2 = vdotr + hubble_a2 * r2;
	      else
		vdotr2 = vdotr;

	      soundspeed_i = sqrt(GAMMA * entropy * pow(densityold, GAMMA_MINUS1));
	      soundspeed_j = sqrt(GAMMA * SphP[p].Entropy * pow(SphP[p].Density, GAMMA_MINUS1));
	      c_ij = 0.5 * (soundspeed_i + soundspeed_j);

	      if(vdotr2 > 0)
		mu_ij = 0;
	      else
		{
		  h_ij = 0.5 * (hsml + PPP[p].Hsml);
		  mu_ij = fac_mu * h_ij * (-vdotr2) / (r2 + 0.0001 * h_ij * h_ij);
		}

	      /*              if((entropy > 10 * SphP[p].Entropy || SphP[p].Entropy > 10 * entropy) && mu_ij < c_ij) */
	      if((entropy > 50 * SphP[p].Entropy || SphP[p].Entropy > 50 * entropy) && mu_ij < c_ij)
		continue;

	      /*
	         if(10*densityold < SphP[p].DensityOld)
	         if(entropy*pow(densityold, GAMMA_MINUS1) > 10 * SphP[p].Entropy*pow(SphP[p].DensityOld, GAMMA_MINUS1))
	         continue; */
	      /* ignore it */
	    }
#endif


#ifdef PERIODIC
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) > (hsml + hdiff))
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) > (hsml + hdiff))
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) > (hsml + hdiff))
	    continue;
#else
	  if(P[p].Pos[0] < (searchmin[0] - hdiff))
	    continue;
	  if(P[p].Pos[0] > (searchmax[0] + hdiff))
	    continue;
	  if(P[p].Pos[1] < (searchmin[1] - hdiff))
	    continue;
	  if(P[p].Pos[1] > (searchmax[1] + hdiff))
	    continue;
	  if(P[p].Pos[2] < (searchmin[2] - hdiff))
	    continue;
	  if(P[p].Pos[2] > (searchmax[2] + hdiff))
	    continue;
#endif
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      printf
		("ThisTask=%d: Need to do a second neighbour loop in hydro-force for (%g|%g|%g) hsml=%g no=%d\n",
		 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
	      *startnode = no;
	      return numngb;
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];
	  hdiff = Extnodes[no].hmax - hsml;
	  if(hdiff < 0)
	    hdiff = 0;

	  no = this->u.d.sibling;	/* in case the node can be discarded */

#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < (searchmin[0] - hdiff))
	    continue;
	  if((this->center[0] - 0.5 * this->len) > (searchmax[0] + hdiff))
	    continue;
	  if((this->center[1] + 0.5 * this->len) < (searchmin[1] - hdiff))
	    continue;
	  if((this->center[1] - 0.5 * this->len) > (searchmax[1] + hdiff))
	    continue;
	  if((this->center[2] + 0.5 * this->len) < (searchmin[2] - hdiff))
	    continue;
	  if((this->center[2] - 0.5 * this->len) > (searchmax[2] + hdiff))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}



/*! This function returns neighbours with distance <= hsml and
 *  returns them in Ngblist. Actually, particles in a box of half side
 *  length hsml are returned, i.e. the reduction to a sphere still
 *  needs to be done in the calling routine.
 */

#ifdef SFR_DECOUPLING
int ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hsml, int *startnode, FLOAT densityold,
			  FLOAT entropy, FLOAT * vel)
#else
int ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
#endif
{
  int k, numngb;
  int no, p;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

#ifdef SFR_DECOUPLING
  FLOAT dx, dy, dz, r;
  double fac_mu, hubble_a, hubble_a2;
  double r2, dvx, dvy, dvz, vdotr, vdotr2, soundspeed_i, soundspeed_j, c_ij, h_ij, mu_ij;
#endif
#ifdef SFR_FEEDBACK
  double a3inv = 1;
  double ne, mu, u, temp;
  double xhyd, yhel;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1.;
  nhot = 0;
  ncold = 0;
#endif

#ifdef SFR_DECOUPLING
  if(All.ComovingIntegrationOn)
    {
      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
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
      hubble_a2 = All.Time * All.Time * hubble_a;
    }
  else
    {
      fac_mu = hubble_a = hubble_a2 = 1;
    }
#endif



  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

#ifdef SFR_METALS
#ifdef SFR_FEEDBACK
	  if(SphP[p].Density > 0)
	    {
	      xhyd = P[p].Zm[6] / P[p].Mass;
	      yhel = (1 - xhyd) / (4. * xhyd);

	      ne = SphP[p].Ne;
	      mu = (1 + 4 * yhel) / (1 + yhel + ne);
	      u = SphP[p].Entropy / GAMMA_MINUS1 * pow(SphP[p].Density * a3inv, GAMMA_MINUS1);	/* energy per mass unit */
	      temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
	      temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */
	    }

	  if(Flag_phase == 1)	/* Hot Phase */
	    {
	      if(temp < All.Tcrit_Phase && SphP[p].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase)
		{
		  ncold++;
		  continue;
		}
	      nhot++;
	    }

	  if(Flag_phase == 2)	/* Cold Phase */
	    {
	      if(!
		 (temp < All.Tcrit_Phase
		  && SphP[p].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase))
		{
		  nhot++;
		  continue;
		}
	      ncold++;
	    }
#endif
#endif


#ifdef SFR_DECOUPLING
	  if(densityold > 0)	/* note: if zero is passed, we have a star (or the first iteration upon start-up) */
	    {
	      dx = searchcenter[0] - P[p].Pos[0];
	      dy = searchcenter[1] - P[p].Pos[1];
	      dz = searchcenter[2] - P[p].Pos[2];
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
	      dvx = vel[0] - SphP[p].VelPred[0];
	      dvy = vel[1] - SphP[p].VelPred[1];
	      dvz = vel[2] - SphP[p].VelPred[2];
	      vdotr = dx * dvx + dy * dvy + dz * dvz;
	      r2 = dx * dx + dy * dy + dz * dz;

	      if(All.ComovingIntegrationOn)
		vdotr2 = vdotr + hubble_a2 * r2;
	      else
		vdotr2 = vdotr;

	      soundspeed_i = sqrt(GAMMA * entropy * pow(densityold, GAMMA_MINUS1));
	      soundspeed_j = sqrt(GAMMA * SphP[p].Entropy * pow(SphP[p].DensityOld, GAMMA_MINUS1));
	      c_ij = 0.5 * (soundspeed_i + soundspeed_j);

	      if(vdotr2 > 0)
		mu_ij = 0;
	      else
		{
		  h_ij = 0.5 * (hsml + PPP[p].Hsml);
		  mu_ij = fac_mu * h_ij * (-vdotr2) / (r2 + 0.0001 * h_ij * h_ij);
		}

	      r = sqrt(r2);

	      /*              if((entropy > 10 * SphP[p].Entropy && mu_ij < c_ij) && SphP[p].DensityOld > All.DensityTailThreshold) */
	      if((entropy > 50 * SphP[p].Entropy && mu_ij < c_ij)
		 && SphP[p].DensityOld > All.DensityTailThreshold)
		continue;
	    }
#endif


#ifdef PERIODIC
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) > hsml)
	    continue;
#else
	  if(P[p].Pos[0] < searchmin[0])
	    continue;
	  if(P[p].Pos[0] > searchmax[0])
	    continue;
	  if(P[p].Pos[1] < searchmin[1])
	    continue;
	  if(P[p].Pos[1] > searchmax[1])
	    continue;
	  if(P[p].Pos[2] < searchmin[2])
	    continue;
	  if(P[p].Pos[2] > searchmax[2])
	    continue;
#endif
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      numngb = ngb_clear_buf(searchcenter, hsml, numngb);
	      if(numngb == MAX_NGB)
		{
		  printf("ThisTask=%d: Need to do a second neighbour loop for (%g|%g|%g) hsml=%g no=%d\n",
			 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
		  *startnode = no;
		  return numngb;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];

	  no = this->u.d.sibling;	/* in case the node can be discarded */
#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < (searchmin[0]))
	    continue;
	  if((this->center[0] - 0.5 * this->len) > (searchmax[0]))
	    continue;
	  if((this->center[1] + 0.5 * this->len) < (searchmin[1]))
	    continue;
	  if((this->center[1] - 0.5 * this->len) > (searchmax[1]))
	    continue;
	  if((this->center[2] + 0.5 * this->len) < (searchmin[2]))
	    continue;
	  if((this->center[2] - 0.5 * this->len) > (searchmax[2]))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}




/*! The buffer for the neighbour list has a finite length MAX_NGB. For
 *  a large search region, this buffer can get full, in which case
 *  this routine can be called to eliminate some of the superfluous
 *  particles in the "corners" of the search box - only the ones in
 *  the inscribed sphere need to be kept.
 */
int ngb_clear_buf(FLOAT searchcenter[3], FLOAT hsml, int numngb)
{
  int i, p;
  FLOAT dx, dy, dz, r2;

  for(i = 0; i < numngb; i++)
    {
      p = Ngblist[i];
#ifdef PERIODIC
      dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
      dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
      dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
#else
      dx = P[p].Pos[0] - searchcenter[0];
      dy = P[p].Pos[1] - searchcenter[1];
      dz = P[p].Pos[2] - searchcenter[2];
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 > hsml * hsml)
	{
	  Ngblist[i] = Ngblist[numngb - 1];
	  i--;
	  numngb--;
	}
    }

  return numngb;
}



/*! Allocates memory for the neighbour list buffer.
 */
void ngb_treeallocate(int npart)
{
  double totbytes = 0;
  size_t bytes;

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

  if(!(Ngblist = malloc(bytes = npart * (long) sizeof(int))))
    {
      printf("Failed to allocate %g MB for ngblist array\n", bytes / (1024.0 * 1024.0));
      endrun(78);
    }
  totbytes += bytes;

  if(ThisTask == 0)
    printf("allocated %f Mbyte for ngb search.\n", totbytes / (1024.0 * 1024.0));
}


/*! free memory allocated for neighbour list buffer.
 */
void ngb_treefree(void)
{
  free(Ngblist);
}

/* This function constructs the neighbour tree. To this end, we
 * actually need to construct the gravitational tree, because we use
 * it now for the neighbour search.  
 */
void ngb_treebuild(void)
{
  if(ThisTask == 0)
    printf("Begin Ngb-tree construction.\n");

  force_treebuild();

  if(ThisTask == 0)
    printf("Ngb-Tree contruction finished \n");
}


int ngb_treefind_darkmatter(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
{
  int k, numngb;
  int no, p;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type != 1)
	    continue;

#ifdef PERIODIC
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) > hsml)
	    continue;
#else
	  if(P[p].Pos[0] < searchmin[0])
	    continue;
	  if(P[p].Pos[0] > searchmax[0])
	    continue;
	  if(P[p].Pos[1] < searchmin[1])
	    continue;
	  if(P[p].Pos[1] > searchmax[1])
	    continue;
	  if(P[p].Pos[2] < searchmin[2])
	    continue;
	  if(P[p].Pos[2] > searchmax[2])
	    continue;
#endif
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      numngb = ngb_clear_buf(searchcenter, hsml, numngb);
	      if(numngb == MAX_NGB)
		{
		  printf("ThisTask=%d: Need to do a second neighbour loop for (%g|%g|%g) hsml=%g no=%d\n",
			 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
		  *startnode = no;
		  return numngb;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];

	  no = this->u.d.sibling;	/* in case the node can be discarded */

#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < searchmin[0])
	    continue;
	  if((this->center[0] - 0.5 * this->len) > searchmax[0])
	    continue;
	  if((this->center[1] + 0.5 * this->len) < searchmin[1])
	    continue;
	  if((this->center[1] - 0.5 * this->len) > searchmax[1])
	    continue;
	  if((this->center[2] + 0.5 * this->len) < searchmin[2])
	    continue;
	  if((this->center[2] - 0.5 * this->len) > searchmax[2])
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}


void ngb_treefind_flagexport(FLOAT searchcenter[3], FLOAT hsml)
{
  int k, no;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  no = Nextnode[no];
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];

	  no = this->u.d.sibling;	/* in case the node can be discarded */

	  if(!(this->u.d.bitflags & 1))	/* if this is not a top-level node, we can skip this branch */
	    continue;

#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < searchmin[0])
	    continue;
	  if((this->center[0] - 0.5 * this->len) > searchmax[0])
	    continue;
	  if((this->center[1] + 0.5 * this->len) < searchmin[1])
	    continue;
	  if((this->center[1] - 0.5 * this->len) > searchmax[1])
	    continue;
	  if((this->center[2] + 0.5 * this->len) < searchmin[2])
	    continue;
	  if((this->center[2] - 0.5 * this->len) > searchmax[2])
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }
}




#ifdef BLACK_HOLES
int ngb_treefind_blackhole(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
{
  int k, numngb;
  int no, p;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

#ifndef REPOSITION_ON_POTMIN
	  if(P[p].Type != 0 && P[p].Type != 5)
	    continue;
#endif

#ifdef PERIODIC
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) > hsml)
	    continue;
#else
	  if(P[p].Pos[0] < searchmin[0])
	    continue;
	  if(P[p].Pos[0] > searchmax[0])
	    continue;
	  if(P[p].Pos[1] < searchmin[1])
	    continue;
	  if(P[p].Pos[1] > searchmax[1])
	    continue;
	  if(P[p].Pos[2] < searchmin[2])
	    continue;
	  if(P[p].Pos[2] > searchmax[2])
	    continue;
#endif
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      numngb = ngb_clear_buf(searchcenter, hsml, numngb);
	      if(numngb == MAX_NGB)
		{
		  printf("ThisTask=%d: Need to do a second neighbour loop for (%g|%g|%g) hsml=%g no=%d\n",
			 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
		  *startnode = no;
		  return numngb;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];

	  no = this->u.d.sibling;	/* in case the node can be discarded */
#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) <
	     (searchmin[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) >
	     (searchmax[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < searchmin[0])
	    continue;
	  if((this->center[0] - 0.5 * this->len) > searchmax[0])
	    continue;
	  if((this->center[1] + 0.5 * this->len) < searchmin[1])
	    continue;
	  if((this->center[1] - 0.5 * this->len) > searchmax[1])
	    continue;
	  if((this->center[2] + 0.5 * this->len) < searchmin[2])
	    continue;
	  if((this->center[2] - 0.5 * this->len) > searchmax[2])
	    continue;
#endif
	  no = this->u.d.nextnode;
	}
    }

  *startnode = -1;
  return numngb;
}
#endif


#if defined(SFR_METALS) && defined(SFR_DECOUPLING) && defined(SFR_PROMOTION)
int ngb_treefind_hotngbs(FLOAT searchcenter[3], FLOAT hsml, int *startnode, FLOAT entropy)
{
  int k, numngb;
  int no, p;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  /*  if(10 * entropy > SphP[p].Entropy) *//* if neighbour is not ignoring us, we ignore it */
	  if(50 * entropy > SphP[p].Entropy)	/* if neighbour is not ignoring us, we ignore it */
	    continue;

#ifdef PERIODIC
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]) > hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]) > hsml)
	    continue;
#else
	  if(P[p].Pos[0] < (searchmin[0]))
	    continue;
	  if(P[p].Pos[0] > (searchmax[0]))
	    continue;
	  if(P[p].Pos[1] < (searchmin[1]))
	    continue;
	  if(P[p].Pos[1] > (searchmax[1]))
	    continue;
	  if(P[p].Pos[2] < (searchmin[2]))
	    continue;
	  if(P[p].Pos[2] > (searchmax[2]))
	    continue;
#endif
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      numngb = ngb_clear_buf(searchcenter, hsml, numngb);
	      if(numngb == MAX_NGB)
		{
		  printf("ThisTask=%d: Need to do a second neighbour loop for (%g|%g|%g) hguess=%g no=%d\n",
			 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
		  *startnode = no;
		  return numngb;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  this = &Nodes[no];

	  no = this->u.d.sibling;

#ifdef PERIODIC
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_LONG_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < (searchmin[0]))
	    continue;
	  if((this->center[0] - 0.5 * this->len) > (searchmax[0]))
	    continue;
	  if((this->center[1] + 0.5 * this->len) < (searchmin[1]))
	    continue;
	  if((this->center[1] - 0.5 * this->len) > (searchmax[1]))
	    continue;
	  if((this->center[2] + 0.5 * this->len) < (searchmin[2]))
	    continue;
	  if((this->center[2] - 0.5 * this->len) > (searchmax[2]))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}
#endif
