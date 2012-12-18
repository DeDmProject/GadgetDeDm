#ifdef COSMIC_RAYS

/* Sanity checks for compiler switches in conjunction with cosmic
 * rays
 */

#if defined(CR_SHOCK) && !defined(COOLING) 
#error Cannot compile with CR_SHOCK but without cooling.
#endif

typedef struct sph_particle_data SphParticle;

/* ============================================================ */
/* ============ Interface functions to GADGET ================= */
/* ============================================================ */

int extern CR_initialize_beta_tabs( double Alpha );
void extern CR_free_beta_tabs ( void );

void extern CR_Particle_Update( SphParticle* Particle );
double extern CR_Particle_Pressure( SphParticle* Particle );

double extern CR_Particle_SpecificEnergy( SphParticle* Particle );
double extern CR_Particle_SpecificNumber( SphParticle* Particle );

double extern CR_Particle_BaryonFraction( SphParticle* Particle );
double extern CR_Particle_MeanKineticEnergy( SphParticle* Particle );

#ifdef COOLING
void extern CR_Particle_GetThermalizationRate( SphParticle* Particle, double* SpecEnergyChangeRate, double* SpecNumberChangeRate );
#endif

void extern CR_Particle_GetDissipationRate( SphParticle* Particle, double* SpecEnergyChangeRate, double* BaryonFractionChangeRate );

double CR_Particle_ThermalizeAndDissipate( SphParticle *Particle,
					   double Time );


#ifdef SFR
double extern CR_Particle_SupernovaFeedback( SphParticle* Particle, double SpecificEnergyInput, double TimeScale );
#endif

#if defined(CR_SHOCK) || defined(CR_SHOCK)
double extern CR_Particle_ShockInject( SphParticle* Particle,
				       double SpecificEnergyInput,
				       double TimeScale );
#endif

void extern CR_Particle_Inject( SphParticle* Particle, double DeltaE, double DeltaN );

double extern CR_q_from_mean_kinetic_energy( double T );
double extern CR_mean_kinetic_energy( double q, double Alpha );


double extern CR_Particle_TimeScaleQ( SphParticle *Particle, double TimeScale );
#endif




