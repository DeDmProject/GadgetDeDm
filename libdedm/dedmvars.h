#define M_p_cmbeasy 3.80917e56    // CMBEASY units Planck Mass
#define M_p         1.094e-54     // Gadget units Planck Mass (solar masses)

extern struct table_interpolation
{
	int npts;
	double *x;
	double *y;
} Table;

#ifdef DEDM_HUBBLE
extern struct hubble_interpolation
{
	int     npts;
	double *a;
	double *hubble;

} HubTable;
#endif

#ifdef DEDM_MASS
extern struct mass_interpolation 
{
	int     npts;
	double *a;
	double *mass;
	double mass_0;

} VariableMassTable;
#endif

#ifdef DEDM_COUPLING
/* Structure for a time-varying beta coupling */
extern struct beta_interpolation 
{
	int npts;
	double *a;
	double *beta;
	double beta_0;

} BetaTable;
#endif

	// TODO OPT VARIABLE_G is not implemented but is left for future possible code extensions
#if defined(DEDM_COUPLING) || defined(VARIABLE_G)
/* Time dependent effective newton constant*/
extern struct g_newton_interpolation 
{
	int npts;
	double *a;
	double *G_new;
	/* G_tilde is simply G modified by a constant factor */
	double G_tilde;
} GTable;

void setBetaZero();

void setG_tilde();
#endif

	// Scott's implementation of the scalar field friction term
extern struct his
{
	int npts;
	double *tau;
	double *a;
} History; // table tau vs. a

extern struct phi
{
	int npts; 
	double *tau;
	double *phi;
} Phi;

extern struct phi_dot
{
	int npts; 
	double *a;
	double *phi;
} PhiDot;
