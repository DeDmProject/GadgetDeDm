#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../allvars.h"
#include "../proto.h"
#include "../forcetree.h"

#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifdef VDE
#include "../libvde/interpolate.h"
#endif

#ifdef DEDM_HUBBLE
#include "../libdedm/interpolate.h"
#endif

#ifdef COOLING

/*
 * This routine does cooling and star formation for
 * the effective multi-phase model.
 */

#ifndef SFR_METALS
#ifndef MHM

#ifndef SFR			/* normal cooling routine when star formation is disabled */
void cooling_and_starformation(void)
{
  int i;
  double dt, dtime, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, dmax1, dmax2;

#ifdef COSMIC_RAYS
  double rCR_dE, rCR_dN;
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
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

      time_hubble_a = All.Time * hubble_a;
    }
  else
    {
      a3inv = time_hubble_a = hubble_a = 1;
    }

  for(i = 0; i < N_gas; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
	  /*  the actual time-step */

	  if(All.ComovingIntegrationOn)
	    dtime = All.Time * dt / time_hubble_a;
	  else
	    dtime = dt;

	  ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

	  unew = DoCooling(DMAX(All.MinEgySpec,
				(SphP[i].Entropy + SphP[i].DtEntropy * dt) /
				GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1)),
			   SphP[i].Density * a3inv, dtime, &ne);

	  SphP[i].Ne = ne;

	  if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
	    {
	      if(dt > 0)
		{

#ifdef COSMIC_RAYS
		  unew += CR_Particle_ThermalizeAndDissipate( SphP + i, dtime );
#endif /* COSMIC_RAYS */

		  SphP[i].DtEntropy = (unew * GAMMA_MINUS1 /
				       pow(SphP[i].Density * a3inv, GAMMA_MINUS1) - SphP[i].Entropy) / dt;

		  if(SphP[i].DtEntropy < -0.5 * SphP[i].Entropy / dt)
		    SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt;

		}
	    }
	}
    }
}

#else

void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, flag, stars_spawned, tot_spawned, stars_converted, tot_converted;
  unsigned int gen, bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew;
  double sum_sm, total_sm, sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2;

#ifdef WINDS
  int j;
  double v;
  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif
#endif
#ifdef METALS
  double w;
#endif

#ifdef COSMIC_RAYS
  double rCR_dE, rCR_dN;
#endif


#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK) 
  double temp, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
#if define(DEDM_HUBBLE) || defined(VDE)
	hubble_a = getH_a(All.Time);
#else
	  hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
#ifdef DARKENERGY
	+ DarkEnergy_a(All.Time);
#else
	+ All.OmegaLambda;
#endif

#endif // DEDM
      hubble_a = All.Hubble * sqrt(hubble_a);

      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    a3inv = ascale = time_hubble_a = 1;



  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

  for(i = 0; i < N_gas; i++)
    {
#ifdef SFR
      if(P[i].Type == 0)
#endif

	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
	    /*  the actual time-step */

	    if(All.ComovingIntegrationOn)
	      dtime = All.Time * dt / time_hubble_a;
	    else
	      dtime = dt;

	    /* check whether conditions for star formation are fulfilled.
	     *  
	     * f=1  normal cooling
	     * f=0  star formation
	     */
	    flag = 1;		/* default is normal cooling */

	    if(SphP[i].Density * a3inv >= All.PhysDensThresh)
	      flag = 0;

	    if(All.ComovingIntegrationOn)
	      if(SphP[i].Density < All.OverDensThresh)
		flag = 1;

#ifdef BLACK_HOLES
	    if(P[i].Mass == 0)
	      flag = 1;
#endif

#ifdef WINDS
	    if(SphP[i].DelayTime > 0)
	      flag = 1;		/* only normal cooling for particles in the wind */

	    if(SphP[i].DelayTime > 0)
	      SphP[i].DelayTime -= dtime;

	    if(SphP[i].DelayTime > 0)
	      if(SphP[i].Density * a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh)
		SphP[i].DelayTime = 0;

	    if(SphP[i].DelayTime < 0)
	      SphP[i].DelayTime = 0;

#endif


#ifdef QUICK_LYALPHA
	    temp = u_to_temp_fac * (SphP[i].Entropy + SphP[i].DtEntropy * dt) /
	      GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);

	    if(SphP[i].Density > All.OverDensThresh && temp < 1.0e5)
	      flag = 0;
	    else
	      flag = 1;
#endif


#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	    if(flag == 1)	/* normal implicit isochoric cooling */
#endif
	      {
		SphP[i].Sfr = 0;

		ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */


		unew = DMAX(All.MinEgySpec,
			    (SphP[i].Entropy + SphP[i].DtEntropy * dt) /
			    GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1));

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK) 
		if(SphP[i].Injected_BH_Energy > 0)
		  {
		    unew += SphP[i].Injected_BH_Energy / P[i].Mass;

		    temp = u_to_temp_fac * unew;
		    if(temp > 5.0e9)
		      unew = 5.0e9 / u_to_temp_fac;

		    SphP[i].Injected_BH_Energy = 0;
		  }
#endif
		unew = DoCooling(unew, SphP[i].Density * a3inv, dtime, &ne);
		SphP[i].Ne = ne;


		if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
		  {
		    /* note: the adiabatic rate has been already added in ! */

		    if(dt > 0)
		      {
#ifdef COSMIC_RAYS
			unew += CR_Particle_ThermalizeAndDissipate( SphP+i, dtime );
#endif /* COSMIC_RAYS */


			SphP[i].DtEntropy = (unew * GAMMA_MINUS1 /
					     pow(SphP[i].Density * a3inv,
						 GAMMA_MINUS1) - SphP[i].Entropy) / dt;

			if(SphP[i].DtEntropy < -0.5 * SphP[i].Entropy / dt)
			  SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt;
		      }
		  }
	      }

	    if(flag == 0)	/* active star formation */
	      {
#if !defined(QUICK_LYALPHA)
		tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * a3inv)) * All.MaxSfrTimescale;

		factorEVP = pow(SphP[i].Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

		egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

		ne = SphP[i].Ne;
		tcool = GetCoolingTime(egyhot, SphP[i].Density * a3inv, &ne);
		SphP[i].Ne = ne;

		y =
		  tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN -
					   (1 - All.FactorSN) * All.EgySpecCold);

		x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

		egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

		if(dt > 0)
		  {
		    if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
		      {
			trelax = tsfr * (1 - x) / x / (All.FactorSN * (1 + factorEVP));
			egycurrent =
			  SphP[i].Entropy * pow(SphP[i].Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;


#ifdef COSMIC_RAYS
			egycurrent += CR_Particle_ThermalizeAndDissipate( SphP + i, dtime );
#endif /* COSMIC_RAYS */


#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
			if(SphP[i].Injected_BH_Energy > 0)
			  {
			    egycurrent += SphP[i].Injected_BH_Energy / P[i].Mass;

			    temp = u_to_temp_fac * egycurrent;

			    if(temp > 5.0e9)
			      egycurrent = 5.0e9 / u_to_temp_fac;

			    if(egycurrent > egyeff)
			      {
				tcool = GetCoolingTime(egycurrent, SphP[i].Density * a3inv, &ne);

				if(tcool < trelax)
				  trelax = tcool;
			      }

			    SphP[i].Injected_BH_Energy = 0;
			  }
#endif



#if !defined(NOISMPRESSURE)
			SphP[i].Entropy =
			  (egyeff +
			   (egycurrent -
			    egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
			  pow(SphP[i].Density * a3inv, GAMMA_MINUS1);

			SphP[i].DtEntropy = 0;
#endif
		      }
		  }

		cloudmass = x * P[i].Mass;

		if(tsfr < dtime)
		  tsfr = dtime;

		sm = (1 - All.FactorSN) * dtime / tsfr * cloudmass;	/* amount of stars expect to form */

		p = sm / P[i].Mass;

		sum_sm += P[i].Mass * (1 - exp(-p));

		SphP[i].Sfr = (1 - All.FactorSN) * cloudmass / tsfr *
		  (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

#ifdef METALS
		w = get_random_number(P[i].ID);
		P[i].Metallicity += w * METAL_YIELD * (1 - exp(-p));
#endif

		prob = P[i].Mass / (All.OrigGasMass / GENERATIONS) * (1 - exp(-p));

#else /* belongs to ifndef(QUICK_LYALPHA) */

		prob = 2.0;	/* this will always cause a star creation event */

#endif /* ends to QUICK_LYALPHA */

		if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		  {
		    if(P[i].Mass < 1.5 * All.OrigGasMass / GENERATIONS)
		      {
			/* here we turn the gas particle itself into a star */
			Stars_converted++;
			stars_converted++;

			sum_mass_stars += P[i].Mass;

			P[i].Type = 4;
#ifdef STELLARAGE
			P[i].StellarAge = All.Time;
#endif
		      }
		    else
		      {
			/* here we spawn a new star particle */
			gen = (int) (1.1 * P[i].Mass / (All.OrigGasMass / GENERATIONS)) - 1;

			for(bits = 0; GENERATIONS > (1 << bits); bits++);

			gen <<= (32 - bits);

			if(NumPart + stars_spawned >= All.MaxPart)
			  {
			    printf
			      ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			       ThisTask, NumPart, stars_spawned, All.MaxPart);
			    fflush(stdout);
			    endrun(8888);
			  }

			P[NumPart + stars_spawned] = P[i];
			P[NumPart + stars_spawned].Type = 4;

			P[NumPart + stars_spawned].ID += gen;

			P[NumPart + stars_spawned].Mass = All.OrigGasMass / GENERATIONS;
			P[i].Mass -= P[NumPart + stars_spawned].Mass;
			sum_mass_stars += P[NumPart + stars_spawned].Mass;
#ifdef STELLARAGE
			P[NumPart + stars_spawned].StellarAge = All.Time;
#endif
			force_add_star_to_tree(i, NumPart + stars_spawned);

			stars_spawned++;
		      }
		  }

#ifdef METALS
		if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		  P[i].Metallicity += (1 - w) * METAL_YIELD * (1 - exp(-p));
#endif

#ifdef COSMIC_RAYS
		CR_Particle_SupernovaFeedback(&SphP[i], 
					      p * All.FeedbackEnergy * All.CR_SNEff,
					      dtime );
#endif


#ifdef WINDS
		/* Here comes the wind model */

		if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
		  {
		    p = All.WindEfficiency * sm / P[i].Mass;

		    prob = 1 - exp(-p);

		    if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		      {
			v =
			  sqrt(2 * All.WindEnergyFraction * All.FactorSN *
			       All.EgySpecSN / (1 - All.FactorSN) / All.WindEfficiency);
#ifdef ISOTROPICWINDS
			theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
			phi = 2 * M_PI * get_random_number(P[i].ID + 4);

			dir[0] = sin(theta) * cos(phi);
			dir[1] = sin(theta) * sin(phi);
			dir[2] = cos(theta);
#else
			dir[0] = P[i].GravAccel[1] * P[i].Vel[2] - P[i].GravAccel[2] * P[i].Vel[1];
			dir[1] = P[i].GravAccel[2] * P[i].Vel[0] - P[i].GravAccel[0] * P[i].Vel[2];
			dir[2] = P[i].GravAccel[0] * P[i].Vel[1] - P[i].GravAccel[1] * P[i].Vel[0];
#endif

			for(j = 0, norm = 0; j < 3; j++)
			  norm += dir[j] * dir[j];

			norm = sqrt(norm);
			if(get_random_number(P[i].ID + 5) < 0.5)
			  norm = -norm;

			if(norm != 0)
			  {
			    for(j = 0; j < 3; j++)
			      dir[j] /= norm;

			    for(j = 0; j < 3; j++)
			      {
				P[i].Vel[j] += v * ascale * dir[j];
				SphP[i].VelPred[j] += v * ascale * dir[j];
			      }

			    SphP[i].DelayTime = All.WindFreeTravelLength / v;
			  }
		      }
		  }
#endif
	      }
	  }

    }				/* end of main loop over active particles */


  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("\n----> spawned %d stars, converted %d gas particles into stars\n\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}


      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
      NumForceUpdate += stars_spawned;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(i = 0, sfrrate = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      sfrrate += SphP[i].Sfr;

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
    }
}


double get_starformation_rate(int i)
{
  double rateOfSF;
  double a3inv;
  int flag;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;



  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  flag = 1;			/* default is normal cooling */

  if(SphP[i].Density * a3inv >= All.PhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(SphP[i].Density * a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP;

  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].Ne;
  tcool = GetCoolingTime(egyhot, SphP[i].Density * a3inv, &ne);

  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);

  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  cloudmass = x * P[i].Mass;

  rateOfSF = (1 - All.FactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}

#endif /* closing of SFR-conditional */

#endif /* closes SFR_MHM conditional */

#endif /* closes SFR_METALS conditional */





#if defined(SFR) || defined(BLACK_HOLES)
void rearrange_particle_sequence(void)
{
  int i, j;
  struct particle_data psave;

#ifdef BLACK_HOLES
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

#ifdef SFR
  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
    }
#endif

#ifdef BLACK_HOLES
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	if(P[i].Type == 0)
	  {
	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];

	    N_gas--;

	    count_gaselim++;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];
	  }

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
#endif

}
#endif /* closing of SFR-conditional */



#ifdef SFR
void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  if(All.PhysDensThresh == 0)
    {
      A0 = All.FactorEVP;

      egyhot = All.EgySpecSN / A0;

      meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

      u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
      u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
      tcool = GetCoolingTime(egyhot, dens, &ne);

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      All.PhysDensThresh =
	x / pow(1 - x,
		2) * (All.FactorSN * All.EgySpecSN - (1 -
						      All.FactorSN) * All.EgySpecCold) /
	(All.MaxSfrTimescale * coolrate);

      if(ThisTask == 0)
	{
	  printf("\nA0= %g  \n", A0);
	  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", All.PhysDensThresh,
		 All.PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs));
	  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
	  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);
	}

      dens = All.PhysDensThresh * 10;

      do
	{
	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  fac = 1 / (log(dens * 1.025) - log(dens));
	  dens *= 1.025;

	  neff = -log(peff) * fac;

	  tsfr = sqrt(All.PhysDensThresh / (dens)) * All.MaxSfrTimescale;
	  factorEVP = pow(dens / All.PhysDensThresh, -0.8) * All.FactorEVP;
	  egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	  ne = 0.5;
	  tcool = GetCoolingTime(egyhot, dens, &ne);

	  y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	  x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));
	  egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	  peff = GAMMA_MINUS1 * dens * egyeff;

	  neff += log(peff) * fac;
	}
      while(neff > 4.0 / 3);

      thresholdStarburst = dens;

#ifdef MODIFIEDBONDI
      All.BlackHoleRefDensity = thresholdStarburst;
      All.BlackHoleRefSoundspeed = sqrt(GAMMA * GAMMA_MINUS1 * egyeff);
#endif


      if(ThisTask == 0)
	{
	  printf("Run-away sets in for dens=%g\n", thresholdStarburst);
	  printf("Dynamic range for quiescent star formation= %g\n", thresholdStarburst / All.PhysDensThresh);
	  fflush(stdout);
	}

      integrate_sfr();

      if(ThisTask == 0)
	{
	  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

	  printf("Isotherm sheet central density: %g   z0=%g\n",
		 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
		 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
	  fflush(stdout);

	}

      if(All.ComovingIntegrationOn)
	{
	  All.Time = All.TimeBegin;
	  IonizeParams();
	}

#ifdef WINDS
      if(All.WindEfficiency > 0)
	if(ThisTask == 0)
	  printf("Windspeed: %g\n",
		 sqrt(2 * All.WindEnergyFraction * All.FactorSN * All.EgySpecSN / (1 - All.FactorSN) /
		      All.WindEfficiency));
#endif
    }
}

void integrate_sfr(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2;
  FILE *fd;


  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      IonizeParams();
    }

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = All.PhysDensThresh; rho <= 1000 * All.PhysDensThresh; rho *= 1.1)
    {
      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;
      tcool = GetCoolingTime(egyhot, rho, &ne);

      y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  if(ThisTask == 0)
    fd = fopen("sfrrate.txt", "w");
  else
    fd = 0;

  for(rho0 = All.PhysDensThresh; rho0 <= 10000 * All.PhysDensThresh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > All.PhysDensThresh)
	    {
	      tsfr = sqrt(All.PhysDensThresh / rho) * All.MaxSfrTimescale;

	      factorEVP = pow(rho / All.PhysDensThresh, -0.8) * All.FactorEVP;

	      egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;
	      tcool = GetCoolingTime(egyhot, rho, &ne);

	      y =
		tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold);
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(All.PhysDensThresh / rho2) * All.MaxSfrTimescale;
	      factorEVP2 = pow(rho2 / All.PhysDensThresh, -0.8) * All.FactorEVP;
	      egyhot2 = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
	      tcool2 = GetCoolingTime(egyhot2, rho2, &ne);
	      y2 =
		tsfr2 / tcool2 * egyhot2 / (All.FactorSN * All.EgySpecSN -
					    (1 - All.FactorSN) * All.EgySpecCold);
	      x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
	      egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
	      P2 = GAMMA_MINUS1 * rho2 * egyeff2;

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;


	      sigma_u4 += rho * dz;
	    }



	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - All.FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
	}
    }


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      IonizeParams();
    }

  if(ThisTask == 0)
    fclose(fd);
}

#endif /* closing of SFR-conditional */


#ifdef SFR
void set_units_sfr(void)
{
  double meanweight, feedbackenergyinergs;

  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

#ifdef INTERNAL_CRIT_DENSITY
  All.PhysDensThresh = All.CritPhysDensity;
#else
  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;
#endif

  All.FeedbackEnergy = All.FactorSN / (1 - All.FactorSN) * All.EgySpecSN;

  feedbackenergyinergs = All.FeedbackEnergy / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);

  if(ThisTask == 0)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n", feedbackenergyinergs);
      printf("OverDensThresh= %g\nPhysDensThresh= %g (internal units)\n", All.OverDensThresh,
	     All.PhysDensThresh);
    }


#ifdef SFR_FEEDBACK
  ESN = All.FactorSN_Phase * EgySNcgs / (All.UnitEnergy_in_cgs);	/*conversion to internal 
									   energy */
  if(ThisTask == 0)
    printf("Feedback energy = %g ergs ,   %g internal units\n", EgySNcgs * All.FactorSN_Phase, ESN);
#endif
}

#endif /* closes SFR */

#endif /* closes COOLING */
