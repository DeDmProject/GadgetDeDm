double get_interpolated_value(double*,double*,int,double);

double critical_density(double);

#ifdef DEDM_HUBBLE
double getH_a(double);
#endif

#ifdef DEDM_COUPLING
double getB_a(double, void*);
#endif

#if defined(DEDM_COUPLING) || defined(VARIABLE_G)
double getG_a(double);
#endif

#ifdef DEDM_MASS
double get_variable_mass_factor(double);
#endif

#ifdef DEDM_DRAG
double getPhiDot_a(double, void*);
#endif

#if defined(DEDM_COUPLING) && defined(DEDM_DRAG)
double get_a_PhiDot_Beta_a(double, void*);
#endif
