/*!
  \file caa_main.c
  \brief Containing the main routines for fitting the age, lga and wgl model.
  \author Geir Storvik

  Fitting the age and lga models are performed in the ::caa_main_model1 routine.

  Fitting the wgl model is performed in the ::caa_main_model2 routine.
  
*/

/*Include Files:*/
#include <R.h>
#include <Rdefines.h>
#include "caa.h"

//#define LOG_FILE 1 !! Use option in Makefile instead


/*(Non-global) External Fuctions:*/

static int initialize_model1(int i_age_errors,double *i_A2A,
			     double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
			     double *i_pri_age_prec_par,double *i_pri_age_ar,
			     double *i_pri_lga_eff_mean,double *i_pri_lga_eff_prec,
			     double *i_pri_lga_prec_par,double *i_pri_lga_ar,
			     int i_lga_fixed_model,int i_lga_cens_model, double *i_lga_cens_par,
			     double *i_cens_mu,double *i_cens_tau,double *i_cens_pri);
static int re_initialize_model1(int i_age_errors,double *i_A2A);
static int initialize_model2(double *i_pri_wgl_eff_mean,double *i_pri_wgl_eff_prec,
			     double *i_pri_wgl_prec_par,double *i_pri_wgl_ar,
			     int i_wgl_fixed_model);
static int re_initialize_model2();
static int init_graph_model1();
static int re_init_graph_model1();
static int MCMC_model1_init();
static int MCMC_model1();
static int MCMC_model1_it(int start_h,int i_force_acc,int i_len_only,
			  int i_it,int *o_acc_h,int i_it_tot);
static int MCMC_model2();
static int MCMC_model2_it(int start_h, int it_tot);
static int MCMC_it_age(int start_h,int i_force_acc,int i_len_only,
		       int i_it,int i_it_tot,int *o_acc_h);
static int MCMC_it_lga_fixed(int start_h,int i_it);
double Qfunc_age_fix(int node, int nnode,char *arg);
double Qfunc_age_fix_new(int node, int nnode,char *arg);
double Qfunc_age_ran(int node, int nnode,char *arg);
double Qfunc_length(int node, int nnode,char *arg);
double Qfunc_length_CC(int node, int nnode,char *arg);
double Qfunc_weight(int node, int nnode,char *arg);
double Qfunc_weight_CC(int node, int nnode,char *arg);

#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif


/*(Non-global) External Variables*/

double           *s_mod1_mean_inv_lik;/*!< Mean of inverse haul-likelihood in age+lga model */
Age_struct       *s_age;             /*!< Current simulations of age-parameters */ 
Age_struct       *s_age_mean;        /*!< Mean of simulations of age-parameters */ 
double           *s_age_mean_inv_lik;/*!< Mean of inverse haul-likelihood in age model */
double           *s_ages_miss_sim;   /*!< Simulated missing ages in last iteration */
LW_struct        *s_length;          /*!< Current simulations of lga-parameters */ 
LW_struct        *s_length_mean;     /*!< Mean of simulations of lga-parameters */ 
double           *s_lga_mean_inv_lik;/*!< Mean of inverse haul-likelihood in lga model */
double           *s_lga_g_a_par_init;/*!< Initial parameters for g-function */
double           *s_lga_fixed_int;   /*!< Values of fixed intercept lga parameter */
double           *s_lga_fixed_slp;   /*!< Values of fixed slope lga parameters */
double           *s_lga_fixed_tau;   /*!< Values of fixed tau_obs lga parameters */
double           *s_lga_fixed_g_a_c;     /*!< Values of fixed g_a c-parameters */
double           *s_lga_fixed_g_a_theta; /*!< Values of fixed g_a theta-parameters */
double           *s_lga_fixed_g_a_gamma; /*!< Values of fixed g_a gamma-parameters */
LW_struct        *s_weight;          /*!< Current simulations of wgl-parameters */ 
LW_struct        *s_weight_mean;     /*!< Mean of simulations of wgl-parameters */ 
double           *s_wgl_mean_inv_lik;/*!< Mean of inverse haul-likelihood in wgl model */
double           *s_wgl_fixed_int;   /*!< Values of fixed intercept wgl parameter */
double           *s_wgl_fixed_slp;   /*!< Values of fixed slope wgl parameters */
double           *s_wgl_fixed_tau;   /*!< Values of fixed tau_obs wgl parameters */
Data_orig        *s_D_orig;          /*!< Original amigo-type data */
Data_age         *s_D_age;           /*!< Age data */
Data_age         *s_D_age_sim;       /*!< Simulated age data */
Data_lin         *s_D_lga;           /*!< Lga data */
Data_lin         *s_D_lga_sim;       /*!< Simulated lga data */
Data_l           *s_D_l;             /*!< Length-only data */
Data_lin         *s_D_wgl;           /*!< Wgl data */
Data_lin         *s_D_wgl_sim;       /*!< Simulated wgl data */
Data_g_a         *s_D_g_a;           /*!< g_a parameters and simulations */
Data_g_a         *s_D_g_a_mean;      /*!< g_a parameters and mean of simulations */
int               s_burn_in;         /*!< Burn in */
int               s_num_it_inner;    /*!< Number of iterations inside main loop */
int               s_num_it_outer;    /*!< Main loop, number of iterations stored */
int               s_constr;          /*!< Constr type, 1 for sum, 2 for treatment */
int               s_seed;            /*!< Seed number */
int              *s_acc_h;           /*!< Vector giving acceptance rates of alpha-param's */
double           *s_ppp;
double           *s_discr;
static int        DO_PVAL=0;         /*!< PVAL not completely implement yet, used for testing */

int               s_coastal_cod;     /*!< =1 if coastal cod, 0 otherwise */
Data_CC          *s_D_CC;              /*!< Parameters for coastal cod */
LW_struct        *s_length_CC;       /*!< Current simulations of lga-parameters for coastal cod */ 
LW_struct        *s_length_CC_mean;  /*!< Mean of simulations of lga-parameters for coastal cod */ 
Data_lin         *s_D_lga_CC;        /*!< Lga data for coastal cod */
Data_g_a         *s_D_g_a_CC;        /*!< g_a parameters and simulations for coastal cod */
Data_g_a         *s_D_g_a_CC_mean;
double           *s_lga_mean_inv_lik_CC;/*!< Mean of inverse haul-likelihood in lga model for coastal cod */
LW_struct        *s_weight_CC;       /*!< Current simulations of wgl-parameters for coastal cod */ 
LW_struct        *s_weight_CC_mean;  /*!< Mean of simulations of wgl-parameters for coastal cod */ 
Data_lin         *s_D_wgl_CC;        /*!< Wgl data for coastal cod */
double           *s_wgl_mean_inv_lik_CC;/*!< Mean of inverse haul-likelihood in wgl model for coastal cod */

Data_COST        *s_D_COST;          /*!< Original data for COST project */
int               s_COST;            /*!< Indicator =1 if COST project, 0 otherwise */





/*!
  \brief Fitting age and lga model through MCMC simulations.
  \author Geir Storvik

  This is a routine made to be called from R using the .Call function.
  Hence, all input parameters are of type SEXP, and the routine also
  returns an SEXP variable. Note that all input variables are read-only!
  The routine therefore starts to convert input data into approperiate
  c-structures (as defined in caa.h).

  Thereafter call the main routine ::MCMC_model1 for performing the MCMC simulations.

  The simulations are stored as a int vector \em o_mcmc with all age parameters
  stored first, then all lga parameters. Simulated variables are only
  stored for every s_num_it_inner simulation after s_burn_in burnin iterations.

  For a further description of the storing of this vector, see the
  ::write_it routine in caa_routines.c

  In addition, the mean of the loglikelihood of the age and lga models
  are returned in o_loglik_mean.

  Residuals from the lga model (from the last iteration) is stored in o_resid_lga.

  The mean of the inverse likelihood (for use of Bayesian CV) for the age and lga models 
  are stored in o_age_mean_inv_lik and o_lga_mean_inv_lik, respectively.

  o_err is equal to 1 if error and zero otherwise.
*/
SEXP caa_main_model1(SEXP i_mcmc_par, SEXP i_constr, SEXP i_seed, 
		     SEXP i_num_par, SEXP i_nBoats, SEXP i_common_par,
		     SEXP i_dataList, SEXP i_ageList, SEXP i_lgaList, SEXP i_priorList,
		     SEXP i_data_COST)
{
  SEXP      elmt = R_NilValue;
  int       n_protect = 0;
  int       i,n;
  long      time_now, time_start;

  double    loglik_age,loglik_lga;

  FILE     *fp;

  /* Variables connected to input data */
  int      *num_par,num_par_COST;
  int       nBoats,nAges;
  int       n_mcmc,n_Fish;

  int      *a_vec,*age_n_cov,*age_ispat,*age_icell;
  int      *age_int_nFac,*age_int_fix,*age_int_c_cov;
  int       age_int_nconstr_cell;
  int      *age_hsz_nFac,*age_hsz_fix,*age_hsz_c_cov;
  int      age_hsz_nconstr_cell;
  int      *age_num_adj_area,*age_adj_area;
  int       age_errors,age_hsz_quad;
  double   *age_int_constr_cell=NULL, *age_hsz_constr_cell=NULL;
  double   *age_int_Sigma_cell=NULL, *age_hsz_Sigma_cell=NULL;
  double   *A2A=NULL;

  int      *lga_n_cov,*lga_ispat,*lga_icell;
  int      *lga_int_nFac,*lga_int_fix,*lga_int_c_cov;
  int       lga_int_nconstr_cell;
  int      *lga_slp_nFac,*lga_slp_fix,*lga_slp_c_cov;
  int       lga_slp_nconstr_cell;
  int      *lga_hsz_nFac,*lga_hsz_fix,*lga_hsz_c_cov;
  int       lga_hsz_nconstr_cell;
  int      *lga_num_adj_area,*lga_adj_area;
  int       lga_g_a_model,lga_g_a_ncat,lga_g_a_nSeason,lga_fixed_model,lga_cens_model;
  int      *lga_g_a_a2Age_vec;
  double   *lga_g_a_avec;
  double   *lga_cens_par=NULL;
  double   *lga_int_constr_cell=NULL,*lga_slp_constr_cell=NULL,*lga_hsz_constr_cell=NULL;
  double   *lga_int_Sigma_cell=NULL,*lga_slp_Sigma_cell=NULL,*lga_hsz_Sigma_cell=NULL;
  double   *cens_mu,*cens_tau,*cens_pri;

  double   *pri_age_eff_mean,*pri_age_eff_prec,*pri_age_prec_par,*pri_age_ar;
  double   *pri_lga_eff_mean,*pri_lga_eff_prec,*pri_lga_prec_par,*pri_lga_ar;

  /* Output variables */
  SEXP      mcmc, mcmc_COST, loglik_mean, resid_lga;
  SEXP      mod1_mean_inv_lik, age_mean_inv_lik, lga_mean_inv_lik;
  SEXP      ppp, err;
  SEXP      resList, resList_names;
  double   *o_mcmc, *o_mcmc_COST, *o_loglik_mean, *o_resid_lga;
  double   *o_mod1_mean_inv_lik, *o_age_mean_inv_lik, *o_lga_mean_inv_lik;
  double   *o_ppp;
  int      *o_err;
  int       n_output = 9;
  char     *res_names[9] = {"mcmc", "mcmc_COST",
			     "loglik_mean", "resid_lga", 
			     "mod1_mean_inv_lik", "age_mean_inv_lik", 
			     "lga_mean_inv_lik", "ppp", "err"};

  time_start = clock();

  /* Connecting input data from R objects to variables */
  /* NOTE: Variables are read-only! */

  num_par = INTEGER_POINTER(i_num_par);
  nBoats = INTEGER_VALUE(i_nBoats);
  nAges = INTEGER_VALUE(getListElement(i_ageList, "nAges"));

  a_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "a_vec")));
  age_n_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "n_cov")));
  age_ispat = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "ispat")));
  age_icell = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "icell")));
  age_int_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "int_nFac")));
  age_int_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "int_fix")));
  age_int_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "int_c_cov")));
  age_int_Sigma_cell = NUMERIC_POINTER(getListElement(i_ageList, "int_Sigma_cell"));
  age_int_constr_cell = NUMERIC_POINTER(getListElement(i_ageList, "int_constr_cell"));
  age_int_nconstr_cell   = INTEGER_VALUE(getListElement(i_ageList, "int_nconstr_cell"));
  age_hsz_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "hsz_nFac")));
  age_hsz_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "hsz_fix")));
  age_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "hsz_c_cov")));
  age_hsz_Sigma_cell = NUMERIC_POINTER(getListElement(i_ageList, "hsz_Sigma_cell"));
  age_hsz_constr_cell = NUMERIC_POINTER(getListElement(i_ageList, "hsz_constr_cell"));
  age_hsz_nconstr_cell = INTEGER_VALUE(getListElement(i_ageList, "hsz_nconstr_cell"));
  age_num_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "num_adj_area")));
  age_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_ageList, "adj_area")));
  age_errors=INTEGER_VALUE(getListElement(i_ageList, "age_errors"));
  age_hsz_quad = INTEGER_VALUE(getListElement(i_ageList, "hsz_quad"));
  if(age_errors)
    A2A = NUMERIC_POINTER(getListElement(i_ageList, "A2A"));

  lga_n_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "n_cov")));
  lga_ispat = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "ispat")));
  lga_icell = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "icell")));
  lga_int_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "int_nFac")));
  lga_int_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "int_fix")));
  lga_int_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "int_c_cov")));
  lga_int_Sigma_cell = NUMERIC_POINTER(getListElement(i_lgaList, "int_Sigma_cell"));
  lga_int_constr_cell = NUMERIC_POINTER(getListElement(i_lgaList, "int_constr_cell"));
  lga_int_nconstr_cell = INTEGER_VALUE(getListElement(i_lgaList, "int_nconstr_cell"));
  lga_slp_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "slp_nFac")));
  lga_slp_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "slp_fix")));
  lga_slp_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "slp_c_cov")));
  lga_slp_Sigma_cell = NUMERIC_POINTER(getListElement(i_lgaList, "slp_Sigma_cell"));
  lga_slp_constr_cell = NUMERIC_POINTER(getListElement(i_lgaList, "slp_constr_cell"));
  lga_slp_nconstr_cell = INTEGER_VALUE(getListElement(i_lgaList, "slp_nconstr_cell"));
  lga_hsz_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "hsz_nFac")));
  lga_hsz_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "hsz_fix")));
  lga_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "hsz_c_cov")));
  lga_hsz_Sigma_cell = NUMERIC_POINTER(getListElement(i_lgaList, "hsz_Sigma_cell"));
  lga_hsz_constr_cell = NUMERIC_POINTER(getListElement(i_lgaList, "hsz_constr_cell"));
  lga_hsz_nconstr_cell = INTEGER_VALUE(getListElement(i_lgaList, "hsz_nconstr_cell"));
  lga_num_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "num_adj_area")));
  lga_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList, "adj_area")));
  lga_fixed_model = INTEGER_VALUE(getListElement(i_lgaList, "fixed_model"));
  lga_g_a_model = INTEGER_VALUE(getListElement(i_lgaList, "g_a_model"));
  lga_g_a_ncat = INTEGER_VALUE(getListElement(i_lgaList,"g_a_ncat"));
  lga_g_a_nSeason = INTEGER_VALUE(getListElement(i_lgaList,"g_a_nSeason"));
  lga_g_a_a2Age_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList,"g_a_a2Age_vec")));
  lga_g_a_avec = NUMERIC_POINTER(getListElement(i_lgaList,"g_a_avec"));
  lga_cens_model = INTEGER_VALUE(getListElement(i_lgaList, "cens_model"));

  if(lga_cens_model)
    lga_cens_par = NUMERIC_POINTER(getListElement(i_lgaList, "cens_par"));


  if(lga_g_a_model)
    s_lga_g_a_par_init = NUMERIC_POINTER(getListElement(i_lgaList, "g_a_par_init"));
  
  if(lga_fixed_model)
    {
      s_lga_fixed_int = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_int"));
      s_lga_fixed_slp = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_slp"));
      s_lga_fixed_tau = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_tau"));
      s_lga_fixed_g_a_c = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_g_a_c"));
      s_lga_fixed_g_a_theta = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_g_a_theta"));
      s_lga_fixed_g_a_gamma = NUMERIC_POINTER(getListElement(i_lgaList, "fixed_g_a_gamma"));
    }

  pri_age_eff_mean = NUMERIC_POINTER(getListElement(i_priorList, "age_eff_mean"));
  pri_age_eff_prec = NUMERIC_POINTER(getListElement(i_priorList, "age_eff_prec"));
  pri_age_prec_par = NUMERIC_POINTER(getListElement(i_priorList, "age_prec_par"));
  pri_age_ar = NUMERIC_POINTER(getListElement(i_priorList, "age_ar"));
  pri_lga_eff_mean = NUMERIC_POINTER(getListElement(i_priorList, "lga_eff_mean"));
  pri_lga_eff_prec = NUMERIC_POINTER(getListElement(i_priorList, "lga_eff_prec"));
  pri_lga_prec_par = NUMERIC_POINTER(getListElement(i_priorList, "lga_prec_par"));
  pri_lga_ar = NUMERIC_POINTER(getListElement(i_priorList, "lga_ar"));


  /* Connecting input variables to global variables */
  s_burn_in = INTEGER_POINTER(i_mcmc_par)[0];
  s_num_it_inner = INTEGER_POINTER(i_mcmc_par)[1];
  s_num_it_outer = INTEGER_POINTER(i_mcmc_par)[2];
  s_constr = INTEGER_VALUE(i_constr);
  s_seed = INTEGER_VALUE(i_seed);

  s_coastal_cod = INTEGER_VALUE(getListElement(i_dataList, "coastal_cod"));
  

  s_COST = INTEGER_VALUE(getListElement(i_data_COST, "COST"));
  if(s_COST)
    {
      cens_mu = NUMERIC_POINTER(getListElement(i_data_COST, "cens.mu"));
      cens_tau = NUMERIC_POINTER(getListElement(i_data_COST, "cens.tau"));
      cens_pri = NUMERIC_POINTER(getListElement(i_data_COST, "cens.pri"));
    }

  /* Allocate space for output */

  n_mcmc = (int) s_num_it_outer*(num_par[0]+num_par[1]+num_par[2]+num_par[3]+num_par[4]);
  PROTECT(mcmc = NEW_NUMERIC(n_mcmc));
  n_protect++;
  o_mcmc = NUMERIC_POINTER(mcmc);

  if(s_COST) 
    {
      num_par_COST = INTEGER_VALUE(getListElement(i_data_COST, "num_par"));
      n_mcmc = (int) s_num_it_outer*num_par_COST;
    }
  else
    n_mcmc = 0;
  PROTECT(mcmc_COST = NEW_NUMERIC(n_mcmc));
  n_protect++;
  o_mcmc_COST = NUMERIC_POINTER(mcmc_COST);

  PROTECT(loglik_mean = NEW_NUMERIC(2));
  n_protect++;
  o_loglik_mean = NUMERIC_POINTER(loglik_mean);

  n_Fish = INTEGER_VALUE(getListElement(i_common_par, "nFish"));
  PROTECT(resid_lga = NEW_NUMERIC(n_Fish));
  n_protect++;
  o_resid_lga = NUMERIC_POINTER(resid_lga);
  
  PROTECT(mod1_mean_inv_lik = NEW_NUMERIC(nBoats));
  n_protect++;
  o_mod1_mean_inv_lik = NUMERIC_POINTER(mod1_mean_inv_lik);
  PROTECT(age_mean_inv_lik = NEW_NUMERIC(nBoats));
  n_protect++;
  o_age_mean_inv_lik = NUMERIC_POINTER(age_mean_inv_lik);

  PROTECT(lga_mean_inv_lik = NEW_NUMERIC(nBoats));
  n_protect++;
  o_lga_mean_inv_lik = NUMERIC_POINTER(lga_mean_inv_lik);
  
  PROTECT(ppp = NEW_NUMERIC(10));
  n_protect++;
  o_ppp = NUMERIC_POINTER(ppp);
  PROTECT(err = NEW_INTEGER(1));
  n_protect++;
  o_err = INTEGER_POINTER(err);

  /* Make "names" attribute for output list */   
  PROTECT(resList_names = allocVector(STRSXP, n_output));
  n_protect++;
  for(i = 0; i < n_output; i++)
    SET_STRING_ELT(resList_names, i,  mkChar(res_names[i]));

  /* Make list with output variables */ 
  PROTECT(resList = allocVector(VECSXP, n_output)); 
  n_protect++;
  SET_VECTOR_ELT(resList, 0, mcmc);   
  SET_VECTOR_ELT(resList, 1, mcmc_COST);
  SET_VECTOR_ELT(resList, 2, loglik_mean);
  SET_VECTOR_ELT(resList, 3, resid_lga);
  SET_VECTOR_ELT(resList, 4, mod1_mean_inv_lik);
  SET_VECTOR_ELT(resList, 5, age_mean_inv_lik);
  SET_VECTOR_ELT(resList, 6, lga_mean_inv_lik);
  SET_VECTOR_ELT(resList, 7, ppp);
  SET_VECTOR_ELT(resList, 8, err);
  setAttrib(resList, R_NamesSymbol, resList_names); //attaching the vector names


  /* Connecting output variables to global variables */
  s_mod1_mean_inv_lik = o_mod1_mean_inv_lik;
  *s_mod1_mean_inv_lik = G_ZERO;
  s_age_mean_inv_lik = o_age_mean_inv_lik;
  *s_age_mean_inv_lik = G_ZERO;
  s_lga_mean_inv_lik = o_lga_mean_inv_lik;
  *s_lga_mean_inv_lik = G_ZERO;
  s_ppp = o_ppp;
  


  #ifdef DEBUG_INPUT
  if(s_COST==0)
    *o_err = write_input_model1_new(i_mcmc_par,i_constr,i_seed,i_num_par,i_nBoats,i_common_par,
				    i_dataList,i_ageList,i_lgaList,i_priorList);
  #endif


  #ifdef LOG_FILE
  g_caa_log = fopen("caa_logfile_model1.txt","w");
  #endif


  /* Make struct for original data */
  if(s_COST==0)
    {
      #ifdef DEBUG_PROG
      printf("Make struct for original data\n");
      #endif
      *o_err = makedata_orig(i_dataList, &s_D_orig);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_orig\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = makedata_CC(i_ageList, &s_D_CC);
	  if(*o_err)
	    {
	      write_warning("caa_main_model1:Error calling makedata_CC\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
    }
  else
    {
      #ifdef DEBUG_PROG
      printf("Make struct for original data\n");
      #endif
      *o_err = makedata_COST(i_data_COST, &s_D_orig, &s_D_COST);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_COST\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      s_D_COST->mcmc = &o_mcmc_COST[0];
      s_D_COST->num_var = num_par_COST;
      #ifdef DEBUG_INPUT
      *o_err = write_input_model1_COST(s_D_orig, s_D_COST,i_ageList,i_lgaList,i_priorList);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling write_COST_data\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      #endif
    }

  /* Make age data */
  #ifdef DEBUG_PROG
  printf("Make age data: makedata_age1\n");
  #endif
  *o_err = makedata_age1(nBoats,nAges,a_vec,
  			 age_n_cov[0],age_int_nFac,age_ispat[0],
  			 age_int_fix,age_int_c_cov,
  			 age_num_adj_area,age_adj_area,
  			 age_n_cov[1],age_hsz_nFac,age_ispat[1],
  			 age_hsz_fix,age_hsz_c_cov,
  			 age_num_adj_area,age_adj_area,
  			 &s_D_age);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_age1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  #ifdef DEBUG_PROG
  printf("Make age data: makedata_age2\n");
  #endif
  *o_err = makedata_age2(s_D_orig,s_D_age);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_age1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  #ifdef DEBUG_PROG
  printf("Make age data: make_cell_constr_age\n");
  #endif
  *o_err = make_cell_constr_age(s_D_age,age_icell,
				age_int_Sigma_cell,age_int_constr_cell,age_int_nconstr_cell,
				age_hsz_Sigma_cell,age_hsz_constr_cell,age_hsz_nconstr_cell);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling make_cell_constr_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Make g_a data */
  #ifdef DEBUG_PROG
  printf("Make g_a data: makedata_g_a\n");
  #endif
  *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			lga_g_a_model,s_D_age,&s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  #ifdef DEBUG_PROG
  printf("Make g_a data: makedata_g_a\n");
  #endif
  *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			lga_g_a_model,s_D_age,&s_D_g_a_mean);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Make extra g_a data if coastal cod */
  if(s_coastal_cod)
    {
      *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			    lga_g_a_model,s_D_age,&s_D_g_a_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_g_a\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			    lga_g_a_model,s_D_age,&s_D_g_a_CC_mean);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_g_a\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  /* Make lga data */
  #ifdef DEBUG_PROG
  printf("Make lga data: makedata_lin1\n");
  #endif
  *o_err = makedata_lin1(nBoats,
			 lga_n_cov[0],lga_int_nFac,lga_ispat[0],lga_int_fix,lga_int_c_cov,
			 lga_num_adj_area,lga_adj_area,
			 lga_n_cov[1],lga_slp_nFac,lga_ispat[1],lga_slp_fix,lga_slp_c_cov,
			 lga_num_adj_area,lga_adj_area,
			 lga_n_cov[2],lga_hsz_nFac,lga_ispat[2],lga_hsz_fix,lga_hsz_c_cov,
			 lga_num_adj_area,lga_adj_area,
			 &s_D_lga);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  s_D_lga->haulweight = s_D_orig->haulweight;
  
  #ifdef DEBUG_PROG
  printf("Make lga data: make_cell_constr_lin\n");
  #endif
  *o_err = make_cell_constr_lin(s_D_lga,lga_icell,
				lga_int_Sigma_cell,lga_int_constr_cell,lga_int_nconstr_cell,
				lga_slp_Sigma_cell,lga_slp_constr_cell,lga_slp_nconstr_cell,
				lga_hsz_Sigma_cell,lga_hsz_constr_cell,lga_hsz_nconstr_cell);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling make_cell_constr_lin\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Make extra lga data if coastal cod */
  if(s_coastal_cod)
    {
      *o_err = makedata_lin1(nBoats,
			     lga_n_cov[0],lga_int_nFac,lga_ispat[0],lga_int_fix,lga_int_c_cov,
			     lga_num_adj_area,lga_adj_area,
			     lga_n_cov[1],lga_slp_nFac,lga_ispat[1],lga_slp_fix,lga_slp_c_cov,
			     lga_num_adj_area,lga_adj_area,
			     lga_n_cov[2],lga_hsz_nFac,lga_ispat[2],lga_hsz_fix,lga_hsz_c_cov,
			     lga_num_adj_area,lga_adj_area,
			     &s_D_lga_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      s_D_lga_CC->haulweight = s_D_orig->haulweight;

      #ifdef DEBUG_PROG
      printf("Make lga data: make_cell_constr_lin\n");
      #endif
      *o_err = make_cell_constr_lin(s_D_lga_CC,lga_icell,
				    lga_int_Sigma_cell,lga_int_constr_cell,lga_int_nconstr_cell,
				    lga_slp_Sigma_cell,lga_slp_constr_cell,lga_slp_nconstr_cell,
				    lga_hsz_Sigma_cell,lga_hsz_constr_cell,lga_hsz_nconstr_cell);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling make_cell_constr_lin\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  
  #ifdef DEBUG_PROG
  printf("Make lga data: makedata_lga_suff\n");
  #endif
  if(s_coastal_cod==0)
    {
      *o_err = makedata_lga_suff(s_D_lga,s_D_orig,s_D_g_a);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_lga_suff\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  else
    {
      *o_err = makedata_lga_suff_CC(s_D_lga,s_D_lga_CC,s_D_orig,s_D_g_a,s_D_g_a_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling makedata_lga_suff_CC\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  #ifdef DEBUG_PROG
  printf("Initialize model1\n");
  #endif
  *o_err = initialize_model1(age_errors,A2A,
			     pri_age_eff_mean,pri_age_eff_prec,pri_age_prec_par,pri_age_ar,
			     pri_lga_eff_mean,pri_lga_eff_prec,pri_lga_prec_par,pri_lga_ar,
			     lga_fixed_model,lga_cens_model,lga_cens_par,
			     cens_mu,cens_tau,cens_pri);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling initialize_model1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  s_age->par->num_var = num_par[0];
  s_age->par->mcmc = &o_mcmc[0];
  s_age->hsz_quad = age_hsz_quad;
  s_length->par->num_var = num_par[1];
  s_length->par->mcmc = &o_mcmc[num_par[0]*s_num_it_outer];
  if(num_par[2] != s_D_g_a->g_a_npar)
    {
      write_warning("Number of parameters in non-linear function do not match\n");
      UNPROTECT(n_protect);
      return(resList);
    }
  s_D_g_a->g_a_mcmc = &o_mcmc[(num_par[0]+num_par[1])*s_num_it_outer];
  if(s_coastal_cod)
    {
      s_length_CC->par->num_var = num_par[3];
      s_length_CC->par->mcmc = &o_mcmc[(num_par[0]+num_par[1]+num_par[2])*s_num_it_outer];
      if(num_par[4] != s_D_g_a_CC->g_a_npar)
	{
	  write_warning("Number of parameters in non-linear function do not match\n");
	  UNPROTECT(n_protect);
	  return(resList);
	}
      s_D_g_a_CC->g_a_mcmc = &o_mcmc[(num_par[0]+num_par[1]+num_par[2]+num_par[3])*s_num_it_outer];
    }

  /* Find node numbers in fixed effect graph for age */
  #ifdef DEBUG_PROG
  printf("Find node numbers: find_node_effect\n");
  #endif
  s_D_age->glm->xcov[0]->n_cov--;
  find_node_effect(s_D_age->glm->xcov,s_D_age->glm->nxcov,s_age->gr_str_f->in_gr,
		   &(s_age->gr_str_f->node));
  s_D_age->glm->xcov[0]->n_cov++;
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling find_node_effect\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Find node numbers in graph for length */
  find_node_effect(s_D_lga->glm->xcov,s_D_lga->glm->nxcov,s_length->gr_str->in_gr,
		   &(s_length->gr_str->node));
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling find_node_effect\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      find_node_effect(s_D_lga_CC->glm->xcov,s_D_lga_CC->glm->nxcov,s_length_CC->gr_str->in_gr,
		   &(s_length_CC->gr_str->node));
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling find_node_effect\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
 
  /* Initialize GMRFLib-type graphs for age and lga models */
  #ifdef DEBUG_PROG
  printf("Initialize GMRFLib-type graphs: init_graph_model1\n");
  #endif
  *o_err = init_graph_model1();
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling init_graph_model1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Starting values of MCMC simulations */
  #ifdef DEBUG_PROG
  printf("Starting values of MCMC simulations: MCMC_model1_init\n");
  #endif
  time_now = clock();
  //fprintf(stderr,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  *o_err = MCMC_model1_init();
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling MCMC_model1_init\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }


  /* Start MCMC simulations */
  if(s_coastal_cod)
    {
      FILE *fp;
      int a,h;
      fp = fopen("N_coastal_cod.txt","w");
      fprintf(fp,"age\n");
      for(h=0;h<s_D_lga->glm->nHaul;h++)
	{
	  fprintf(fp,"%d ",h);
	  for(a=0;a<s_D_age->glm->ncat;a++)	    
	    fprintf(fp,"%d ",s_D_age->Ages[h][a]);
	  fprintf(fp,"\n");
	}
      fprintf(fp,"lga\n");
      for(h=0;h<s_D_lga->glm->nHaul;h++)
	{
	  fprintf(fp,"%d ",h);
	  for(a=0;a<s_D_age->glm->ncat;a++)	    
	    fprintf(fp,"%d ",s_D_age->Ages[h][a]);
	  fprintf(fp,"\n");
	}
      fclose(fp);
    }

  #ifdef DEBUG_PROG
  printf("Start MCMC simulations: MCMC_model1\n");
  #endif
  time_now = clock();
  //fprintf(stderr,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  *o_err = MCMC_model1();
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling MCMC_model1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  if(s_COST==0)
    {
      /* Calculate likelihood over averaged parameters */
      #ifdef DEBUG_PROG
      printf("Calculate likelihood over averaged parameters\n");
      #endif
      *o_err = calc_lik_age(s_age_mean,s_D_age,&loglik_age);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling calc_lik_age\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      *o_err = calc_lik_lin(s_length_mean,s_D_lga,&loglik_lga);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling calc_lik_lin\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = calc_lik_lin(s_length_CC_mean,s_D_lga_CC,&loglik_lga);
	  if(*o_err)
	    {
	      write_warning("caa_main_model1:Error calling calc_lik_lin\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
      o_loglik_mean[0] = loglik_age;
      o_loglik_mean[1] = loglik_lga;

      #ifdef DEBUG_PROG
      printf("caa_main.c:Calculate likelihood over averaged parameters for coastal cod?\n");
      #endif

      *o_err = calc_resid_lga(s_D_orig->totage,s_D_orig->totlength,s_D_orig->nFishBoat,
			      s_D_orig->replength,
			      s_age,s_D_age,s_length,s_D_lga,s_D_g_a,o_resid_lga);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling calc_resid_lga\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  else
    {
      #ifdef DEBUG_PROG
      printf("Calculate likelihood over averaged parameters for COST?\n");
      #endif
    }

  // Clean up
  #ifdef DEBUG_PROG
  printf("Clean up\n");
  #endif
  time_now = clock();
  //fprintf(stderr,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  /* Initialize GMRFLib-type graphs for age and lga models */
  *o_err = re_init_graph_model1();
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_init_graph_model1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  s_D_age->glm->xcov[0]->n_cov--;
  re_find_node_effect(s_D_age->glm->xcov,s_D_age->glm->nxcov,s_age->gr_str_f->in_gr,
		   &(s_age->gr_str_f->node));
  s_D_age->glm->xcov[0]->n_cov++;

  re_find_node_effect(s_D_lga->glm->xcov,s_D_lga->glm->nxcov,s_length->gr_str->in_gr,
		   &(s_length->gr_str->node));
  if(s_coastal_cod)
    re_find_node_effect(s_D_lga_CC->glm->xcov,s_D_lga_CC->glm->nxcov,s_length_CC->gr_str->in_gr,
			&(s_length_CC->gr_str->node));


  *o_err = re_initialize_model1(age_errors,A2A);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  if(s_coastal_cod==0)
    {
      *o_err = re_makedata_lga_suff(s_D_lga);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_lga_suff\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  else
    {
      *o_err = re_makedata_lga_suff_CC(s_D_lga,s_D_lga_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_lga_suff_CC\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  *o_err = re_makedata_lin1(nBoats,
			    lga_n_cov[0],lga_int_nFac,lga_ispat[0],
			    lga_int_fix,lga_int_c_cov,
			    lga_num_adj_area,lga_adj_area,
			    lga_n_cov[1],lga_slp_nFac,lga_ispat[1],
			    lga_slp_fix,lga_slp_c_cov,
			    lga_num_adj_area,lga_adj_area,
			    lga_n_cov[2],lga_hsz_nFac,lga_ispat[2],
			    lga_hsz_fix,lga_hsz_c_cov,
			    lga_num_adj_area,lga_adj_area,
			    &s_D_lga);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      *o_err = re_makedata_lin1(nBoats,
				lga_n_cov[0],lga_int_nFac,lga_ispat[0],
				lga_int_fix,lga_int_c_cov,
				lga_num_adj_area,lga_adj_area,
				lga_n_cov[1],lga_slp_nFac,lga_ispat[1],
				lga_slp_fix,lga_slp_c_cov,
				lga_num_adj_area,lga_adj_area,
				lga_n_cov[2],lga_hsz_nFac,lga_ispat[2],
				lga_hsz_fix,lga_hsz_c_cov,
				lga_num_adj_area,lga_adj_area,
				&s_D_lga_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  *o_err = re_makedata_g_a(&s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  *o_err = re_makedata_g_a(&s_D_g_a_mean);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  
  if(s_coastal_cod)
    {
      *o_err = re_makedata_g_a(&s_D_g_a_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_g_a\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      *o_err = re_makedata_g_a(&s_D_g_a_CC_mean);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_g_a\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    } 
  
  *o_err = re_makedata_age2(s_D_age);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_age2\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  *o_err = re_makedata_age1(nBoats,nAges,a_vec,
  			    age_n_cov[0],age_int_nFac,age_ispat[0],
  			    age_int_fix,age_int_c_cov,
  			    age_num_adj_area,age_adj_area,
  			    age_n_cov[1],age_hsz_nFac,age_ispat[1],
  			    age_hsz_fix,age_hsz_c_cov,
  			    age_num_adj_area,age_adj_area,
  			    &s_D_age);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_age1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  if(s_COST==0)
    {
      *o_err = re_makedata_orig(&s_D_orig);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_orig\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = re_makedata_CC(&s_D_CC);
	  if(*o_err)
	    {
	      write_warning("caa_main_model1:Error calling re_makedata_CC\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
    }
  else
    {
      *o_err = re_makedata_COST(&s_D_orig,&s_D_COST);
      if(*o_err)
	{
	  write_warning("caa_main_model1:Error calling re_makedata_COST\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }


  #ifdef LOG_FILE
  fclose(g_caa_log);
  #endif


  UNPROTECT(n_protect);
  time_now = clock();
  //fprintf(stderr,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);

  return(resList);
}		/* end of caa_main_model1 */
  


/*!
  \brief Fitting wgl model through MCMC simulations.
  \author Geir Storvik

  Since this is a routine made to be called from R or Splus, all input 
  variables are pointers to vectors. The routine therefore starts to convert 
  input data into approperiate c-structures (as defined in caa.h).
  Thereafter call the main routine ::MCMC_model2 for performing the MCMC simulations.

  The simulations are stored as a int vector \em o_mcmc with all age parameters
  stored first, then all lga parameters. Simulated variables are only
  stored for every s_num_it_inner simulation after s_burn_in burnin iterations.

  For a further description of the storing of this vector, see the
  ::write_it routine in caa_routines.c

  In addition, the mean of the loglikelihood of the wgl models is returned in o_loglik_mean.

  Residuals from the wgl model (from the last iteration) is stored in o_resid_wgl.

  The mean of the inverse likelihood (for use of Bayesian CV) for the wgl model
  is stored in o_wgl_mean_inv_lik.

  o_err is equal to 1 if error and zero otherwise.
*/
void caa_main_model2(int *i_mcmc_par,
		     int *i_constr,int *i_seed,int *i_coastal_cod,
		     int *i_wgl_nBoats,double *i_totlength,double *i_totweight,double *i_haulweight,
		     int *i_replength,int *i_tottype,int *i_nFishBoat,
		     int *i_n_cov,int *i_ispat,int *i_wgl_icell,
		     int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
		     double *i_wgl_int_Sigma_cell,double *i_wgl_int_constr_cell,
		     int *i_wgl_int_nconstr_cell,
		     int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
		     double *i_wgl_slp_Sigma_cell,double *i_wgl_slp_constr_cell,
		     int *i_wgl_slp_nconstr_cell,
		     int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
		     double *i_wgl_hsz_Sigma_cell,double *i_wgl_hsz_constr_cell,
		     int *i_wgl_hsz_nconstr_cell,
		     int *i_wgl_fixed_model,
		     double *i_wgl_fixed_int,double *i_wgl_fixed_slp,double *i_wgl_fixed_tau,
		     int *i_num_adj_area,int *i_adj_area,
		     int *i_num_par,
		     double *i_pri_wgl_eff_mean,double *i_pri_wgl_eff_prec,
		     double *i_pri_wgl_prec_par,double *i_pri_wgl_ar,
		     double *o_mcmc,double *o_loglik_mean,double *o_resid_wgl,
                     double *o_wgl_mean_inv_lik,
		     int *o_err)
{
  double      loglik_wgl;
  int         wgl_int_nconstr_cell,wgl_slp_nconstr_cell,wgl_hsz_nconstr_cell;

  #ifdef DEBUG_INPUT
  int         err;
  err = write_input_model2(i_mcmc_par,i_constr,i_seed,
			   i_wgl_nBoats,i_totlength,i_totweight,i_haulweight,
			   i_replength,i_nFishBoat,
			   i_n_cov,i_ispat,
			   i_wgl_int_nFac,i_wgl_int_fix,i_wgl_int_c_cov,
			   i_wgl_slp_nFac,i_wgl_slp_fix,i_wgl_slp_c_cov,
			   i_wgl_hsz_nFac,i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			   i_num_adj_area,i_adj_area,i_num_par);
  #endif

  #ifdef LOG_FILE
  g_caa_log = fopen("caa_logfile_model2.txt","w");
  #endif

  s_burn_in = (int) i_mcmc_par[0];
  s_num_it_inner = (int) i_mcmc_par[1];
  s_num_it_outer = (int) i_mcmc_par[2];
  s_constr = (int) *i_constr;
  s_seed = (int) *i_seed;
  s_wgl_fixed_int = i_wgl_fixed_int;
  s_wgl_fixed_slp = i_wgl_fixed_slp;
  s_wgl_fixed_tau = i_wgl_fixed_tau;
  s_wgl_mean_inv_lik = o_wgl_mean_inv_lik;

  s_coastal_cod = (int) *i_coastal_cod;
  #ifdef DEBUG_PROG
  printf("Coastal cod = %d\n",s_coastal_cod);
  #endif


  /* Read weight given length data */
  #ifdef DEBUG_PROG
  printf("Read weight given length data\n");
  #endif
  *o_err = makedata_lin1(*i_wgl_nBoats,
			 i_n_cov[0],i_wgl_int_nFac,i_ispat[0],i_wgl_int_fix,i_wgl_int_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[1],i_wgl_slp_nFac,i_ispat[1],i_wgl_slp_fix,i_wgl_slp_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[2],i_wgl_hsz_nFac,i_ispat[2],i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			 i_num_adj_area,i_adj_area,
			 &s_D_wgl);
  s_D_wgl->haulweight = i_haulweight;
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  /* Make extra wgl data if coastal cod */
  if(s_coastal_cod)
    {
      *o_err = makedata_lin1(*i_wgl_nBoats,
			     i_n_cov[0],i_wgl_int_nFac,i_ispat[0],i_wgl_int_fix,i_wgl_int_c_cov,
			     i_num_adj_area,i_adj_area,
			     i_n_cov[1],i_wgl_slp_nFac,i_ispat[1],i_wgl_slp_fix,i_wgl_slp_c_cov,
			     i_num_adj_area,i_adj_area,
			     i_n_cov[2],i_wgl_hsz_nFac,i_ispat[2],i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			     i_num_adj_area,i_adj_area,
			     &s_D_wgl_CC);
      s_D_wgl_CC->haulweight = i_haulweight;
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }
  #ifdef DEBUG_PROG
  printf("Make constraints\n");
  #endif
  wgl_int_nconstr_cell = (int) *i_wgl_int_nconstr_cell;
  wgl_slp_nconstr_cell = (int) *i_wgl_slp_nconstr_cell;
  wgl_hsz_nconstr_cell = (int) *i_wgl_hsz_nconstr_cell;
  *o_err = make_cell_constr_lin(s_D_wgl,i_wgl_icell,
				i_wgl_int_Sigma_cell,i_wgl_int_constr_cell,wgl_int_nconstr_cell,
				i_wgl_slp_Sigma_cell,i_wgl_slp_constr_cell,wgl_slp_nconstr_cell,
				i_wgl_hsz_Sigma_cell,i_wgl_hsz_constr_cell,wgl_hsz_nconstr_cell);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling make_cell_constr_lin\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  if(s_coastal_cod)
    {
      *o_err = make_cell_constr_lin(s_D_wgl_CC,i_wgl_icell,
				    i_wgl_int_Sigma_cell,i_wgl_int_constr_cell,wgl_int_nconstr_cell,
				    i_wgl_slp_Sigma_cell,i_wgl_slp_constr_cell,wgl_slp_nconstr_cell,
				    i_wgl_hsz_Sigma_cell,i_wgl_hsz_constr_cell,wgl_hsz_nconstr_cell);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling make_cell_constr_lin\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }

  #ifdef DEBUG_PROG
  printf("Make sufficient statistics\n");
  #endif
  if(s_coastal_cod==0)
    {
      *o_err = makedata_wgl_suff(s_D_wgl,i_nFishBoat,i_totlength,i_totweight,
				    i_replength,i_haulweight);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling makedata_lin2\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }
  else
    {
      *o_err = makedata_wgl_suff_CC(s_D_wgl,s_D_wgl_CC,i_nFishBoat,i_totlength,i_totweight,
				    i_replength,i_haulweight,i_tottype);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling makedata_lin2\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }

  #ifdef DEBUG_PROG
  printf("Initialize model2\n");
  #endif
  *o_err = initialize_model2(i_pri_wgl_eff_mean,i_pri_wgl_eff_prec,i_pri_wgl_prec_par,i_pri_wgl_ar,*i_wgl_fixed_model);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling initialize_model2\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  s_weight->par->num_var = i_num_par[0];
  s_weight->par->mcmc = &o_mcmc[0];
  if(s_coastal_cod)
    {
      s_weight_CC->par->num_var = i_num_par[1];
      s_weight_CC->par->mcmc = &o_mcmc[i_num_par[0]*s_num_it_outer];
    }


  /* Find node numbers in graph for weight */
  #ifdef DEBUG_PROG
  printf("Find node numbers in graph for weight\n");
  #endif
  find_node_effect(s_D_wgl->glm->xcov,s_D_wgl->glm->nxcov,s_weight->gr_str->in_gr,
		   &(s_weight->gr_str->node));
  if(s_coastal_cod)
    find_node_effect(s_D_wgl_CC->glm->xcov,s_D_wgl_CC->glm->nxcov,s_weight_CC->gr_str->in_gr,
		     &(s_weight_CC->gr_str->node));
  
  /* Run MCMC */
  #ifdef DEBUG_PROG
  printf("Run MCMC\n");
  #endif
  *o_err = MCMC_model2();
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling MCMC\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  /* Calculate likelihood over averaged parameters */
  #ifdef DEBUG_PROG
  printf("Calculate likelihood over averaged parameters\n");
  #endif
  *o_err = calc_lik_lin(s_weight_mean,s_D_wgl,&loglik_wgl);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling calc_lik_lin\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  o_loglik_mean[0] = loglik_wgl;
  
  *o_err = calc_resid_wgl(i_totlength,i_totweight,i_nFishBoat,s_weight_mean,s_D_wgl,o_resid_wgl);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling calc_resid_wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  // Clean up
  #ifdef DEBUG_PROG
  printf("Clean up\n");
  #endif
  re_find_node_effect(s_D_wgl->glm->xcov,s_D_wgl->glm->nxcov,s_weight->gr_str->in_gr,
		   &(s_weight->gr_str->node));
  if(s_coastal_cod)
    re_find_node_effect(s_D_wgl_CC->glm->xcov,s_D_wgl_CC->glm->nxcov,s_weight_CC->gr_str->in_gr,
			&(s_weight_CC->gr_str->node));
  
  *o_err = re_initialize_model2();
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling re_initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  if(s_coastal_cod==0)
    {
      *o_err = re_makedata_wgl_suff(s_D_wgl,i_nFishBoat,i_totlength,i_totweight,
				    i_replength);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling re_makedata_lin2\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }
  else
    {
      *o_err = re_makedata_wgl_suff_CC(s_D_wgl,s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling re_makedata_lin2\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }
  
  *o_err = re_makedata_lin1(*i_wgl_nBoats,
			    i_n_cov[0],i_wgl_int_nFac,i_ispat[0],
			    i_wgl_int_fix,i_wgl_int_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[1],i_wgl_slp_nFac,i_ispat[1],
			    i_wgl_slp_fix,i_wgl_slp_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[2],i_wgl_hsz_nFac,i_ispat[2],
			    i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			    i_num_adj_area,i_adj_area,
			    &s_D_wgl);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  if(s_coastal_cod)
    {
      *o_err = re_makedata_lin1(*i_wgl_nBoats,
				i_n_cov[0],i_wgl_int_nFac,i_ispat[0],
				i_wgl_int_fix,i_wgl_int_c_cov,
				i_num_adj_area,i_adj_area,
				i_n_cov[1],i_wgl_slp_nFac,i_ispat[1],
				i_wgl_slp_fix,i_wgl_slp_c_cov,
				i_num_adj_area,i_adj_area,
				i_n_cov[2],i_wgl_hsz_nFac,i_ispat[2],
				i_wgl_hsz_fix,i_wgl_hsz_c_cov,
				i_num_adj_area,i_adj_area,
				&s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_main_model2:Error calling re_makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  return;
	}
    }

  #ifdef LOG_FILE
  fclose(g_caa_log);
  #endif

  return;
}		/* end of caa_main_model2 */
  


/*!
  \brief Memory allocation and initial values for age and lga models
  \author Geir Storvik
*/
static int initialize_model1(int i_age_errors,double *i_A2A,
			     double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
			     double *i_pri_age_prec_par,double *i_pri_age_ar,
			     double *i_pri_lga_eff_mean,double *i_pri_lga_eff_prec,
			     double *i_pri_lga_prec_par,double *i_pri_lga_ar,
			     int i_lga_fixed_model,int i_lga_cens_model,double *i_lga_cens_par,
			     double *i_cens_mu,double *i_cens_tau,double *i_cens_pri)
{
  int err;
  int i;

  (*GMRFLib_uniform_init)(s_seed);
  /* Allocating memory for age structs */
  err = init_age(s_D_age,i_age_errors,i_A2A,s_coastal_cod,
		 i_pri_age_eff_mean,i_pri_age_eff_prec,i_pri_age_prec_par,i_pri_age_ar,
		 &s_age,&s_age_mean);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_age\n");
      return(err);
    }
  s_age->gr_str_f->Qfunc = Qfunc_age_fix;
  s_age->gr_str_f->Qfunc_new = Qfunc_age_fix_new;
  s_age->gr_str_r->Qfunc = Qfunc_age_ran;


  /* Allocating memory for length structs */
  err = init_lin(s_D_lga,i_pri_lga_eff_mean,i_pri_lga_eff_prec,i_pri_lga_prec_par,i_pri_lga_ar,
		 &s_length,&s_length_mean);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_lin\n");
      return(err);
    }
  s_length->gr_str->Qfunc = Qfunc_length;
  s_length->fixed_model = i_lga_fixed_model;
  /* Allocating memory for extra length structs if coastal cod */
  if(s_coastal_cod)
    {
      err = init_lin(s_D_lga_CC,i_pri_lga_eff_mean,i_pri_lga_eff_prec,
		     i_pri_lga_prec_par,i_pri_lga_ar,&s_length_CC,&s_length_CC_mean);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_lin\n");
	  return(err);
	}
      s_length_CC->gr_str->Qfunc = Qfunc_length_CC;
      s_length_CC->fixed_model = i_lga_fixed_model;
    }

  /* Initialize non-linear model */
  if(s_D_g_a->g_a_model>0)
    {
      if(s_length->fixed_model == 0)
	{
	  if(s_D_g_a->g_a_model==1) /* Schnute Richards model */
	    {
	      for(i=0;i<s_D_g_a->g_a_npar;i++)
		{
		  s_D_g_a->g_a_par[i] = s_lga_g_a_par_init[i];
		  //printf("g_a_par[%d]=%f\n",i,s_D_g_a->g_a_par[i]);
		}
	    }
	  else if(s_D_g_a->g_a_model==2) /* polyn3 model */
	    {
	      s_D_g_a->g_a_par[0] = G_ZERO;  // log(beta1)
	      s_D_g_a->g_a_par[1] = G_ZERO;  // beta3
	    }
	}
      else /* use fixed lga model (Schnute Richards model) */
	{
	  if(s_D_g_a->g_a_model!=1)
	    {
	      write_warning("initialize_model1:Error using fixed model. Unknown g-function.\n");
	      return(1);
	    }
	  s_D_g_a->g_a_par[0] = s_lga_fixed_g_a_c[0];
	  s_D_g_a->g_a_par[1] = s_lga_fixed_g_a_theta[0];
	  s_D_g_a->g_a_par[2] = s_lga_fixed_g_a_gamma[0];
	}
    }
  if(s_coastal_cod)
    {
      if(s_D_g_a_CC->g_a_model>0)
	{
	  if(s_length_CC->fixed_model == 0)
	    {
	      if(s_D_g_a_CC->g_a_model==1) /* Schnute Richards model */
		{
		  for(i=0;i<s_D_g_a_CC->g_a_npar;i++)
		    s_D_g_a_CC->g_a_par[i] = s_lga_g_a_par_init[i];
		}
	      else if(s_D_g_a_CC->g_a_model==2) /* polyn3 model */
		{
		  s_D_g_a_CC->g_a_par[0] = G_ZERO;  // log(beta1)
		  s_D_g_a_CC->g_a_par[1] = G_ZERO;  // beta3
		}
	    }
	  else /* use fixed lga model (Schnute Richards model) */
	    {
	      if(s_D_g_a_CC->g_a_model!=1)
		{
		  write_warning("initialize_model1:Error using fixed model. Unknown g-function.\n");
		  return(1);
		}
	      fprintf(stderr,"initialize_model1: correct fixed parameters for c,theta and gamma\n");
	      s_D_g_a_CC->g_a_par[0] = s_lga_fixed_g_a_c[0];
	      s_D_g_a_CC->g_a_par[1] = s_lga_fixed_g_a_theta[1];
	      s_D_g_a_CC->g_a_par[2] = s_lga_fixed_g_a_gamma[2];
	    }
	}
    }

  /* Initial values if sampling discards */
  s_length->cens_model = i_lga_cens_model;
  if(s_length->cens_model)
    {
      s_length->cens_k = i_lga_cens_par[0];
      s_length->cens_m = i_lga_cens_par[1];
      s_length->cens_r = i_lga_cens_par[2];
      s_length->cens_Nlim = i_lga_cens_par[3];
    }
  if(s_coastal_cod)
    {
      s_length_CC->cens_model = i_lga_cens_model;
      if(s_length_CC->cens_model)
	{
	  s_length_CC->cens_k = i_lga_cens_par[0];
	  s_length_CC->cens_m = i_lga_cens_par[1];
	  s_length_CC->cens_r = i_lga_cens_par[2];
	  s_length_CC->cens_Nlim = i_lga_cens_par[3];
	}
    }

  if(s_COST == 1)
    {
      err = init_cens_COST(s_D_orig,s_D_COST,i_cens_mu,i_cens_tau,i_cens_pri);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_cens_COST\n");
	  return(err);
	}
    }
  
  /* Initial values for length given age model, needed sampling ages */
  err = init_lga_par(s_length,s_D_lga,s_lga_fixed_int,s_lga_fixed_slp,s_lga_fixed_tau);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_lga_par\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = init_lga_par(s_length_CC,s_D_lga_CC,s_lga_fixed_int,s_lga_fixed_slp,s_lga_fixed_tau);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_lga_par\n");
	  return(err);
	}
    }

  err = init_evaluate(s_D_age->glm->nHaul,s_D_age->glm->ncat);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_evaluate\n");
      return(err);
    }

  for(i=0;i<2;i++)
    {
      s_ppp[i] = G_ZERO;
      //s_discr[i] = G_ZERO;
    }
  /* Struct for age simulations */
  if(DO_PVAL)
    {
      s_D_age_sim = CALLOC(1,Data_age);              // Free ok
      err = init_glm_sim(s_D_age->glm,&s_D_age_sim->glm);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_glm_sim\n");
	  return(err);
	}
      s_D_age_sim->n_h =s_D_age->n_h;

      s_D_age_sim->Ages = Mmatrix_2d(0,s_D_age->glm->nHaul-1,    // Free ok
				     0,s_D_age->glm->ncat-1,sizeof(int),1);

      /* Struct for length simulations */
      s_D_lga_sim = CALLOC(1,Data_lin);        // Free ok
      /* Free ok */
      err = init_glm_sim(s_D_lga->glm,&s_D_lga_sim->glm);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_glm_sim\n");
	  return(err);
	}

    }
  s_acc_h = CALLOC(s_D_age->glm->nHaul,int);     // Free ok

  // Initialization for routines in caa_sample_multi.c
  err = sample_multi_initialize(s_D_age->glm->ncat);
  if(err)
    {
      write_warning("initialize_model1:Error calling sampling_multi_initialize\n");
      return(err);
    }

  // Initialization for routines in caa_sample_g_a.c
  if(s_D_g_a->g_a_model>0)
    {
      err = sample_g_a_initialize(s_D_g_a->ncat,s_D_g_a->g_a_model);
      if(err)
	{
	  write_warning("initialize_model1:Error calling sampling_g_a_initialize\n");
	  return(err);
	}
    }


  return(0);
}		/* end of initialize_model1 */



/*!
  \brief Reallocate memory allocated in initialize_model1
  \author Geir Storvik
*/
static int re_initialize_model1(int i_age_errors,double *i_A2A)
{
  int  err;

  err = sample_multi_re_initialize();
  if(err)
    {
      write_warning("re_initialize_model1:Error calling sampling_multi_re_initialize\n");
      return(err);
    }
  if(s_D_g_a->g_a_model>0)
    {
      err = sample_g_a_re_initialize();
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling sampling_g_a_re_initialize\n");
	  return(err);
	}
    }

  err = re_init_age(s_D_age,i_age_errors,i_A2A,&s_age,&s_age_mean);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling re_init_age\n");
      return(err);
    }

  err = re_init_lin(s_D_lga,&s_length,&s_length_mean);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling re_init_lin\n");
      return(err);
    }

  if(s_coastal_cod)
    {
      err = re_init_lin(s_D_lga,&s_length_CC,&s_length_CC_mean);
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling re_init_lin\n");
	  return(err);
	}
    }

  err = re_init_evaluate(s_D_age->glm->ncat);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling init_evaluate\n");
      return(err);
    }

  if(DO_PVAL)
    {
      err = re_init_glm_sim(s_D_age->glm,&s_D_age_sim->glm);
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling re_init_glm_sim\n");
	  return(err);
	}
      FREE(s_D_age_sim); 
      Fmatrix_2d(&s_D_age_sim->Ages[0][0],&s_D_age_sim->Ages[0]);
      FREE(s_D_age_sim);
      err = re_init_glm_sim(s_D_lga->glm,&s_D_lga_sim->glm);
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling re_init_glm_sim\n");
	  return(err);
	}
      FREE(s_D_lga_sim);
    }

  FREE(s_acc_h);

  return(0);
}		/* end of re_initialize_model1 */




/*!
  \brief Memory allocation and initial values for wgl model
  \author Geir Storvik
*/
static int initialize_model2(double *i_pri_wgl_eff_mean,double *i_pri_wgl_eff_prec,
			     double *i_pri_wgl_prec_par,double *i_pri_wgl_ar,
			     int i_wgl_fixed_model)
{
  int err;
  
  (*GMRFLib_uniform_init)(s_seed);
  
  /* Allocating memory for weight structs */
  err = init_lin(s_D_wgl,i_pri_wgl_eff_mean,i_pri_wgl_eff_prec,
		 i_pri_wgl_prec_par,i_pri_wgl_ar,
		 &s_weight,&s_weight_mean);
  if(err)
    {
      write_warning("initialize_model2:Error calling init_lin\n");
      return(err);
    }
  
  s_weight->gr_str->Qfunc = Qfunc_weight;
  s_weight->fixed_model = i_wgl_fixed_model;
  
  /* Allocating memory for extra weight structs if coastal cod */
  if(s_coastal_cod)
    {
      err = init_lin(s_D_wgl_CC,i_pri_wgl_eff_mean,i_pri_wgl_eff_prec,
		     i_pri_wgl_prec_par,i_pri_wgl_ar,
		     &s_weight_CC,&s_weight_CC_mean);
      if(err)
	{
	  write_warning("initialize_model2:Error calling init_lin\n");
	  return(err);
	}
      
      s_weight_CC->gr_str->Qfunc = Qfunc_weight_CC;
      s_weight_CC->fixed_model = i_wgl_fixed_model;
    }
  
  if(s_weight->fixed_model == 1)
    {
      s_weight->par->tau_obs = s_wgl_fixed_tau[0];
      /* Assume first covariate is constant term */
      if(s_D_wgl->glm->xcov[0]->n_fac[0]==1)
	s_weight->par->eff[0][0][0][0] = s_wgl_fixed_int[0]; /* Intercept */
      else
  	{
  	  write_warning("initialize_model2:First covariate should be constant term\n");
  	  return(1);
  	}
      if(s_D_wgl->glm->xcov[1]->n_fac[0]==1)
  	s_weight->par->eff[0][1][0][0] = s_wgl_fixed_slp[0]; /* Slope */
      else
  	{
  	  write_warning("initialize_model2:First covariate should be constant term\n");
  	  return(1);
  	}
    }
  #ifdef DEBUG_PROG
  printf("Init_wgl_par:\n");
  printf("Int = %f\n",s_weight->par->eff[0][0][0][0]);
  printf("Slp = %f\n",s_weight->par->eff[0][1][0][0]);
  #endif
  
  err = init_evaluate(s_D_wgl->glm->nHaul,1);
  if(err)
    {
      write_warning("initialize_model2:Error calling init_evaluate\n");
      return(err);
    }
  
  if(DO_PVAL)
    {
      /* Struct for weight simulations */
      s_D_wgl_sim = CALLOC(1,Data_lin);    // Free ok
      /* Free ok */
      err = init_glm_sim(s_D_wgl->glm,&s_D_wgl_sim->glm);  
      if(err)
	{
	  write_warning("initialize_model2:Error calling init_glm_sim\n");
	  return(err);
	}
    }

  return(0);
}		/* end of initialize_model2 */




/*!
  \brief Reallocate memory allocated in initialize_model2
  \author Geir Storvik
*/
static int re_initialize_model2()
{
  int   err;

  err = re_init_lin(s_D_wgl,&s_weight,&s_weight_mean);
  if(err)
    {
      write_warning("re_initialize_model2:Error calling re_init_lin\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = re_init_lin(s_D_wgl_CC,&s_weight_CC,&s_weight_CC_mean);
      if(err)
	{
	  write_warning("re_initialize_model2:Error calling re_init_lin\n");
	  return(err);
	}
    }


  err = re_init_evaluate(1);
  if(err)
    {
      write_warning("initialize_model2:Error calling re_init_evaluate\n");
      return(err);
    }

  if(DO_PVAL)
    {
      FREE(s_D_wgl_sim);
      err = re_init_glm_sim(s_D_wgl->glm,&s_D_wgl_sim->glm);
      if(err)
	{
	  write_warning("re_initialize_model2:Error calling re_init_glm_sim\n");
	  return(err);
	}
    }
  return(0);
}		/* end of re_initialize_model2 */



/*!
  \brief Initialize graph structures for age and lga model (for GMRFLib)
  \author Geir Storvik

  For both structures, Precision matrix calculated, Constraints are defined
  and Graph structure initiated by GMRLFib_init_problem is made
*/
static int init_graph_model1()
{
  int err;

  #ifdef DEBUG_PROG
  printf("Initializing for age\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing for age\n");
  #endif

  // Remove haul effect before making graph
  s_D_age->glm->xcov[0]->n_cov--;
  /* Initialization for graph of fixed effects */
  err = make_graph_gauss(s_D_age->glm,s_age->gr_str_f,s_age->par,0,s_constr,s_age->hsz_quad);
  if(err)
    {
      write_warning("init_graph_model1:Error calling make_graph_gauss\n");
      return(err);
    }
  // Add haul effect again
  s_D_age->glm->xcov[0]->n_cov++;
  s_age->gr_str_f->x_new = CALLOC(s_age->gr_str_f->graph->n,double);  // Free ok
  s_age->gr_str_f->Q_new = Mmatrix_2d(0,s_age->gr_str_f->graph->n-1,  // Free ok
				      0,s_age->gr_str_f->graph->n-1,sizeof(double),1);
  /* Initialization for graph of random effects for use of GMRFLib */
  // This part is currently not used because other routines specialized for
  // our problem is used instead
  if(0)
    {
      err = make_graph_age_ran(s_age,s_D_age);
      if(err)
	{
	  write_warning("init_graph_model1:Error calling make_graph_age_ran\n");
	  return(err);
	}
  }

  #ifdef DEBUG_PROG
  printf("Initializing for length\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing for length\n");
  #endif
  err = make_graph_gauss(s_D_lga->glm,s_length->gr_str,s_length->par,0,s_constr,0);
  if(err)
    {
      write_warning("init_graph_model1:Error calling make_graph_gauss\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = make_graph_gauss(s_D_lga_CC->glm,s_length_CC->gr_str,s_length_CC->par,0,s_constr,0);
      if(err)
	{
	  write_warning("init_graph_model1:Error calling make_graph_gauss\n");
	  return(err);
	}
    }

  return(0);
}


/*!
  \brief Initialize graph structures for age and lga model (for GMRFLib)
  \author
*/
static int re_init_graph_model1()
{
  int err;

  // Clean up
  // Remove haul effect before making graph
  s_D_age->glm->xcov[0]->n_cov--;
  /* Initialization for graph of fixed effects */
  err = re_make_graph_gauss(s_age->gr_str_f);
  if(err)
    {
      write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
      return(err);
    }
  // Add haul effect again
  s_D_age->glm->xcov[0]->n_cov++;
  FREE(s_age->gr_str_f->x_new);
  Fmatrix_2d(&s_age->gr_str_f->Q_new[0][0],&s_age->gr_str_f->Q_new[0]);
  
  err = re_make_graph_age_ran(s_age,s_D_age);
  if(err)
    {
      write_warning("re_init_graph_model1:Error calling re_make_graph_age_ran\n");
      return(err);
    }
  
  err = re_make_graph_gauss(s_length->gr_str);
  if(err)
    {
      write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = re_make_graph_gauss(s_length_CC->gr_str);
      if(err)
	{
	  write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
	  return(err);
	}
    }
  
  return(0);
}



/*!
  \author Geir Storvik
  \brief  Initialization for MCMC simulations of age and lga model

  Start with initialization of parameters. This is mainly based on using data
  which is complete (that is amigo-type data and age-stratified by length data)
  and running a few simulations on this.

  Note 
*/
static int MCMC_model1_init()
{
  int        err,printFile;

  printFile = 0;  
  #ifdef DEBUG_PROG
  printFile = 1;
  #endif

  FILE *file;
  int a,h;
  int *N;
  if(0)
    {
      N = CALLOC(s_D_age->glm->ncat,int);      // Free ok
      for(a=0;a<s_D_age->glm->ncat;a++)
	N[a] = 0;
      file=fopen("nAges.txt","w");
      for(h=0;h<s_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<s_D_age->glm->ncat;a++)
	    N[a] += s_D_age->Ages_fix[h][a];
	}
      for(a=0;a<s_D_age->glm->ncat;a++)
	fprintf(file,"%d ",N[a]);
      fprintf(file,"\n");
    }

  /* Sample ages using only age-given length transitions */
  #ifdef DEBUG_PROG
  printf("Sample ages init\n");
  #endif
  if(s_COST)
    {
      if(s_D_COST->obs->n_trip==0 || s_D_COST->mland->n_trip==0)
	err = sample_ages_init2_COST(s_D_orig,s_age,s_D_age,s_length,s_D_lga,s_D_g_a,
				     s_D_COST,printFile);
      else
	err = sample_ages_init_COST(s_D_orig,s_age,s_D_age,s_length,s_D_lga,s_D_g_a,
				    s_D_COST,printFile);
    }
  else if(s_D_orig->coastal_cod)
    err = sample_ages_len_only_init_new_CC(s_D_orig,s_D_CC,s_age,s_D_age,
					   s_D_lga,s_D_g_a,s_D_lga_CC,s_D_g_a_CC,printFile);
  else
    err = sample_ages_len_only_init_new(s_D_orig,s_age,s_D_age,s_length,
					s_D_lga,s_D_g_a,printFile);
  if(err)
    {
      write_warning("MCMC_model1_init:Error calling sample_ages_len_only_init_new\n");
      return(err);
    }

  if(0)
    {
      for(a=0;a<s_D_age->glm->ncat;a++)
	N[a] = 0;
      for(h=0;h<s_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<s_D_age->glm->ncat;a++)
	    N[a] += s_D_age->Ages[h][a];
	}
      for(a=0;a<s_D_age->glm->ncat;a++)
	fprintf(file,"%d ",N[a]);
      fprintf(file,"\n");
      fclose(file);
    }
  FREE(N);


  /* Sample discarded */
  if(s_COST)
    {
      if(s_D_COST->mland->n_trip>0)
	{
	  if(s_D_COST->obs->n_trip==0)
	    err = sample_discard_init2_COST(s_age,s_D_age,s_length,s_D_lga,s_D_g_a,
					    s_D_orig,s_D_COST,printFile);
	  else
	    err = sample_discard_init_COST(s_age,s_D_age,s_length,s_D_lga,s_D_g_a,
					   s_D_orig,s_D_COST,printFile);
	  if(err)
	    {
	      write_warning("MCMC_model1_init:Error calling sample_discard_init_COST\n");
	      return(err);
	    }
	}
      err = sample_lambda_init_COST(s_D_age,s_D_COST);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling sample_lambda_init_COST\n");
	  return(err);
	}
    }

  // Initialization alpha-values for full data
  //Using a small precision when optimization gives alpha-values close
  //to the observations and seems to result in realisations with high
  //posterior densities. It seems however that convergence of the haul-effects
  //goes very slowly and that it is better to not do this. This choice
  //is therefore commented out, but should be investigated further later on.
  //s_age->par->tau_obs = 0.0001;
  err = age_haul_modes(0,s_D_age->glm->nHaul,s_age,s_D_age);
  //s_age->par->tau_obs = 1;

  /* Then initial fit of lga model */
  #ifdef DEBUG_PROG
  printf("Initial fit of lga model\n");
  #endif
  if(s_length->fixed_model == 0)
    {
      err = MCMC_it_lga(s_D_lga,s_D_g_a,s_length,0);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling MCMC_it_lga\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
	  err = MCMC_it_lga(s_D_lga_CC,s_D_g_a_CC,s_length_CC,0);
	  if(err)
	    {
	      write_warning("MCMC_model1_init:Error calling MCMC_it_lga\n");
	      return(err);
	    }
	}
    }
  else /* use fixed lga model */
    {
      err = MCMC_it_lga_fixed(0,0);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling MCMC_it_lga_fixed\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
	  err = MCMC_it_lga_fixed(0,0);
	  if(err)
	    {
	      write_warning("MCMC_model1_init:Error calling MCMC_it_lga_fixed\n");
	      return(err);
	    }
	}
    }

  // Sample parameters in non-linear model
  if(0)
    {
      //Test om denne er nødvendig
      fprintf(stderr,"Test of suff_g_a_init er nodvendig her\n");
      err = suff_g_a_init(s_length,s_D_age,s_D_lga,s_D_g_a,0,s_D_g_a->suff);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling suff_g_a_init\n");
	  return(err);
	}
      err = sample_g_a(s_length,s_D_age,s_D_lga,s_D_g_a,s_D_g_a->suff,0);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling sample_g_a\n");
	  return(err);
	}
    }

  return(0);
}



/*!
  \author Geir Storvik
  \brief  Perform MCMC simulations of age and lga model

  The simulations are divided into burn-in and simulations after burnin in which
  num_it_inner simulations are performed for each num_it_outer simulation.
  All inner simulations are performed through the ::MCMC_model1_it routine.

  Simulations are saved for each num_it_outer iteration using the ::write_it
  routine.
*/
static int MCMC_model1()
{
  int        err,force_acc;
  int        it,iti,ito;
  int        it_tot;
  double     ppp,discr,discr2,ppp2[2];
  long       time_now, time_start;

  time_start = clock();

  force_acc = 0;
  it_tot = 0;
  for(iti=0;iti<s_burn_in;iti++)
    {
      err = MCMC_model1_it(0,force_acc,1,iti,s_acc_h,it_tot);
      it_tot++;
      if(err)
	{
	  write_warning("MCMC_model1:Error calling MCMC_model1_it\n");
	  return(err);
	}
      force_acc = 0;
    }
  /* Start full simulation */
  it = 0;
  force_acc = 0;
  for(ito=0;ito<s_num_it_outer;ito++)
    {
      //printf("\n\nouter it=%d\n",ito);
      //time_now = clock();
      //fprintf(stderr,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
      for(iti=0;iti<s_num_it_inner;iti++)
	{
          #ifdef DEBUG_PROG
	    printf("\n\nIteration %d\n",it);
          #endif
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\n\nIteration %d\n",it);
          #endif
	  //if(s_COST && ito%10 == 0)
	  if(s_COST)
	    {
	      //printf("resample length data\n");
	      err = resample_data_COST(s_D_orig,s_D_COST);
	      if(err)
		{
		  write_warning("MCMC_model1:Error calling resample_data_COST\n");
		  return(err);
		}
	    }
	  err = MCMC_model1_it(0,force_acc,1,it,s_acc_h,it_tot);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling MCMC_model1_it\n");
	      return(err);
	    }
	  err = update_average_age(it,s_age,s_D_age,s_age_mean);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_age\n");
	      return(err);
	    }

	  err = update_average_g_a(it,s_D_g_a,s_D_g_a_mean);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_g_a\n");
	      return(err);
	    }
	  if(s_D_orig->coastal_cod)
	    {
	      err = update_average_g_a(it,s_D_g_a_CC,s_D_g_a_CC_mean);
	      if(err)
		{
		  write_warning("MCMC_model1:Error calling update_average_g_a\n");
		  return(err);
		}
	    }


	  err = update_average_lin(it,s_length,s_D_lga,s_length_mean);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_lin\n");
	      return(err);
	    }
	  if(s_D_orig->coastal_cod)
	    {
	      err = update_average_lin(it,s_length_CC,s_D_lga_CC,s_length_CC_mean);
	      if(err)
		{
		  write_warning("MCMC_model1:Error calling update_average_lin\n");
		  return(err);
		}
	    }

          err = Bayes_CV_model1(it,s_age,s_D_age,s_length,s_D_lga,
			      s_mod1_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          err = Bayes_CV_age(it,s_age,s_D_age,s_age_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          err = Bayes_CV_lin(it,s_length,s_D_lga,s_lga_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

	  if(s_D_orig->coastal_cod)
	    {
              #ifdef DEBUG_PROG
	      printf("Calculate s_lga_mean_inv_lik for coastal_cod\n");
	      #endif
	    }

          it++;
	  it_tot++;
          force_acc = 0;
	}
      if(0)
	{
	  FILE *fp;
	  int h,a;
	  double mu;
	  fp = fopen("haul_param.txt","a");
	  for(h=0;h<s_D_age->glm->nHaul;h++)
	    {
	      for(a=0;a<s_D_age->glm->ncat-1;a++)
		{
		  mu = calc_eff(s_D_age->glm->xcov[0],s_age->par->eff[a][0],h)*s_D_age->glm->suff[h][0][0];
		  fprintf(fp,"%f ",s_age->alpha[h][a]-mu);
		}
	      a=s_D_age->glm->ncat-1;
	      mu = calc_eff(s_D_age->glm->xcov[0],s_age->par->eff[a][0],h)*s_D_age->glm->suff[h][0][0];
	      fprintf(fp,"%f\n",s_age->alpha[h][a]-mu);
	    }
	  fclose(fp);
	}
      
      /* Save samples */
      s_D_age->glm->xcov[0]->n_cov--;
      err= write_it(ito,s_D_age->glm,s_age->par);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling write_it\n");
	  return(err);
	}
      s_D_age->glm->xcov[0]->n_cov++;

      err = write_it(ito,s_D_lga->glm,s_length->par);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling write_it\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
	  err = write_it(ito,s_D_lga_CC->glm,s_length_CC->par);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling write_it, lga parameters\n");
	      return(err);
	    }
	}

      err = write_it_g_a(ito,s_D_g_a);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling write_it_g_a\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
	  err = write_it_g_a(ito,s_D_g_a_CC);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling write_it_g_a\n");
	      return(err);
	    }
	}

      if(s_COST)
	{
	  err = write_it_COST(ito,s_D_COST);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling write_it_COST\n");
	      return(err);
	    }
	}

      err = calc_KS_age(s_age,s_D_age,&discr,&ppp);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling calc_KS_age\n");
	  return(err);
	}
      err = update_mean(&s_ppp[0],ppp,ito);
      //err = update_mean(&s_discr[0],discr,ito);
      //err = update_mean(&s_ppp[2],discr,ito);
      if(0)
	{
	  err = calc_D_Robins_age(s_age,s_D_age,&discr,&ppp);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling calc_D_Robins_age\n");
	      return(err);
	    }
	  err = update_mean(&s_ppp[1],ppp,ito);
	  //err = update_mean(&s_ppp[3],discr,ito);
	  //err = update_mean(&s_discr[1],discr,ito);
	}
      err = calc_entropy_age(s_age,s_D_age,&discr,&discr2,ppp2);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling calc_entr_age\n");
	  return(err);
	}
      err = update_mean(&s_ppp[1],ppp2[0],ito);
      err = update_mean(&s_ppp[2],ppp2[1],ito);

      if(DO_PVAL)
	{
	  /*
	  err = Simulate_length_weight(s_age,s_D_age,
		 s_length,s_D_lga_sim,s_weight,s_D_wgl_sim);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Simulate_length_weight\n");
	      return(err);
	    }
	  */
	  /*
	  s_discr_obs = calc_discr(s_age,s_length,s_weight,
				   s_D_age,s_D_lga,s_D_wgl);
	  s_discr_new = calc_discr(s_age,s_length,s_weight,
				   s_D_age_sim,s_D_lga_sim,s_D_wgl_sim);
	  */
	}
    }

  return(0);
}		/* end of MCMC_model1 */



/*!
  \author Geir Storvik
  \brief Main simulation steps inside each MCMC iteration for model1

  Simulations are performed by switching between the following steps:
  - Simulating age parameters
  - Simulating parameters of non-linear function in lga model
  - Simulating lga parameters

  If fixed_model, then simulate only age parameters.
*/
static int MCMC_model1_it(int start_h,int i_force_acc,int i_len_only,
			  int i_it,int *o_acc_h,int i_it_tot)
{
  int err=0;

  #ifdef DEBUG_PROG
    printf("\nSampling Age-parameters\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"\nSampling Age-parameters\n");
  #endif

  err = MCMC_it_age(start_h,i_force_acc,i_len_only,i_it,i_it_tot,o_acc_h);
  if(err)
    {
      write_warning("MCMC_model1_it:Error calling MCMC_it_age\n");
      return(err);
    }
  if(s_length->fixed_model == 0) /* Sample lga model */
    {
      #ifdef DEBUG_PROG
      printf("\nSampling Length-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling Length-parameters\n");
      #endif

      err = MCMC_it_g_a(s_D_age,s_D_lga,s_D_g_a,s_length,start_h,i_it);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_g_a\n");
          #endif
	  write_warning("MCMC_model1_it:Error calling MCMC_it_g_a\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
	  err = MCMC_it_g_a(s_D_age,s_D_lga_CC,s_D_g_a_CC,s_length_CC,start_h,i_it);
	  if(err)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_g_a\n");
              #endif
	      write_warning("MCMC_model1_it:Error calling MCMC_it_g_a\n");
	      return(err);
	    }
	}
      err = MCMC_it_lga(s_D_lga,s_D_g_a,s_length,start_h);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga\n");
          #endif
	  write_warning("MCMC_model1_it:Error calling MCMC_it_lga\n");
	  return(err);
	}
      if(s_D_orig->coastal_cod)
	{
          #ifdef DEBUG_PROG
	  printf("\n Coastal cod\n");
	  printf("tau[%d][%d]=%lf\n",0,1,s_length_CC->par->tau[0][1]);
	  printf("tau_obs=%lf\n",s_length_CC->par->tau_obs);
          #endif
	  //s_length_CC->par->tau[0][1] = s_length->par->tau[0][1];
	  //s_length_CC->par->tau_obs = s_length->par->tau_obs;
	  err = MCMC_it_lga(s_D_lga_CC,s_D_g_a_CC,s_length_CC,start_h);
	  if(err)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga\n");
              #endif
	      write_warning("MCMC_model1_it:Error calling MCMC_it_lga\n");
	      return(err);
	    }
	}
    }
  else /* Use fixed lga model */
    {
      #ifdef DEBUG_PROG
      printf("\nSampling fixed Length-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling fixed Length-parameters\n");
      #endif

      err = MCMC_it_lga_fixed(start_h,i_it_tot);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
          #endif
	  write_warning("MCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
	  return(err);
	}
    }

  return(0);
}		/* end of MCMC_model1_it */



/*!
  \author Geir Storvik
  \brief  Perform MCMC simulations of wgl model

  Start with initialization of parameters.

  The simulations are divided into burn-in and simulations after burnin in which
  num_it_inner simulations are performed for each num_it_outer simulation.
  All inner simulations are performed through the ::MCMC_model1_it routine.

  Simulations are saved for each num_it_outer iteration using the ::write_it
  routine.
*/
static int MCMC_model2()
{
  int        err;
  int        it,iti,ito,it_tot;

  #ifdef DEBUG_PROG
    printf("Initializing for weight\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing for weight\n");
  #endif
  err = make_graph_gauss(s_D_wgl->glm,s_weight->gr_str,s_weight->par,0,s_constr,0);
  if(err)
    {
      write_warning("MCMC_model2:Error calling make_graph_gauss\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = make_graph_gauss(s_D_wgl_CC->glm,s_weight_CC->gr_str,s_weight_CC->par,0,s_constr,0);
      if(err)
	{
	  write_warning("MCMC_model2:Error calling make_graph_gauss\n");
	  return(err);
	}
    }
 
  it_tot=0;
  /* Burn in */
  for(iti=0;iti<s_burn_in;iti++)
    {
      err = MCMC_model2_it(0,it_tot);
      it_tot++;
     if(err)
       {
	 write_warning("MCMC_model2:Error calling MCMC_model2_it\n");
	 return(err);
       }
    }

  /* Start full simulation */
  it = 0;
  for(ito=0;ito<s_num_it_outer;ito++)
    {
      for(iti=0;iti<s_num_it_inner;iti++)
	{
          #ifdef DEBUG_PROG
	    printf("\n\nIteration %d\n",it);
          #endif
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\n\nIteration %d\n",it);
          #endif
	  err = MCMC_model2_it(0,it_tot);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling MCMC_model2_it\n");
	      return(err);
	    }

	  err = update_average_lin(it,s_weight,s_D_wgl,s_weight_mean);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling update_average_lin\n");
	      return(err);
	    }

	  if(s_coastal_cod)
	    {
	      err = update_average_lin(it,s_weight_CC,s_D_wgl_CC,s_weight_CC_mean);
	      if(err)
		{
		  write_warning("MCMC_model2:Error calling update_average_lin\n");
		  return(err);
		}
	    }

          err = Bayes_CV_lin(it,s_weight,s_D_wgl,s_wgl_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          it++;
	  it_tot++;
	}

      /* Calculate likelihoods */
      err = calc_lik_lin(s_weight,s_D_wgl,&(s_weight->par->loglik));
      if(err)
	{
	  write_warning("MCMC_model2:Error calling calc_lik_lin\n");
	  return(err);
	}
      if(s_coastal_cod)
	{
	  err = calc_lik_lin(s_weight_CC,s_D_wgl_CC,&(s_weight_CC->par->loglik));
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling calc_lik_lin\n");
	      return(err);
	    }
	}
      #ifdef DEBUG_PROG
      printf("Loglik_wgl=%lf\n",s_weight->par->loglik);
      printf("Int=%lf\n",s_weight->par->eff[0][0][0][0]);
      printf("Slope=%lf\n",s_weight->par->eff[0][1][0][0]);
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"Loglik_wgl=%lf\n",s_weight->par->loglik);
      fprintf(g_caa_log,"Int=%lf\n",s_weight->par->eff[0][0][0][0]);
      fprintf(g_caa_log,"Slope=%lf\n",s_weight->par->eff[0][1][0][0]);
      #endif

      /* Save samples */
      err = write_it(ito,s_D_wgl->glm,s_weight->par);
      if(err)
	{
	  write_warning("MCMC_model2:Error calling write_it\n");
	  return(err);
	}

      if(s_coastal_cod)
	{
	  err = write_it(ito,s_D_wgl_CC->glm,s_weight_CC->par);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling write_it coastal cod\n");
	      return(err);
	    }
	}

      if(DO_PVAL)
	{
	  /*
	  err = Simulate_length_weight(s_age,s_D_age,
		 s_length,s_D_lga_sim,s_weight,s_D_wgl_sim);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling Simulate_length_weight\n");
	      return(err);
	    }
	  */
	  /*
	  s_discr_obs = calc_discr(s_age,s_length,s_weight,
				   s_D_age,s_D_lga,s_D_wgl);
	  s_discr_new = calc_discr(s_age,s_length,s_weight,
				   s_D_age_sim,s_D_lga_sim,s_D_wgl_sim);
	  */
	}
    }
  /*
  it2++;
  err = calc_lik_lin(s_weight_mean,s_D_wgl,&(s_weight_mean->par->loglik));
  if(err)
    {
      write_warning("MCMC_model2:Error calling calc_lik_lin\n");
      return(err);
    }
  
  err = write_it_mean(it2,s_D_wgl->glm,s_weight->par,s_weight_mean->par,1);
  if(err)
    {
      write_warning("MCMC_model2:Error calling write_it_mean\n");
      return(err);
    }
  */

  // Clean up
  err = re_make_graph_gauss(s_weight->gr_str);
  if(err)
    {
      write_warning("MCMC_model2:Error calling re_make_graph_gauss\n");
      return(err);
    }

  if(s_coastal_cod)
    {
      err = re_make_graph_gauss(s_weight_CC->gr_str);
      if(err)
	{
	  write_warning("MCMC_model2:Error calling re_make_graph_gauss\n");
	  return(err);
	}
    }
  return(0);
}		/* end of MCMC_model2 */




/*!
  \author Geir Storvik
  \brief Main simulation steps inside each MCMC iteration for model2

  Simulations are performed by the following step:
  - Simulating wgl parameters
*/
static int MCMC_model2_it(int start_h, int it_tot)
{
  int err;

  if(s_weight->fixed_model == 0) /* Sample wgl model */
    {
      #ifdef DEBUG_PROG
      printf("Sampling Weight-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"Sampling Weight-parameters\n");
      #endif
      err = MCMC_it_wgl(s_D_wgl,s_weight,start_h);
      if(err)
	{
	  write_warning("MCMC_model2_it:Error calling MCMC_it_wgl\n");
	  return(err);
	}
      if(s_coastal_cod)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"Sampling Weight-parameters: Coastal cod\n");
          #endif
	  err = MCMC_it_wgl(s_D_wgl_CC,s_weight_CC,start_h);
	  if(err)
	    {
	      write_warning("MCMC_model2_it:Error calling MCMC_it_wgl\n");
	      return(err);
	    }
	}
    }  
  else
    {
      s_weight->par->tau_obs = s_wgl_fixed_tau[it_tot];
      s_weight->par->eff[0][0][0][0] = s_wgl_fixed_int[it_tot]; /* Intercept */
      s_weight->par->eff[0][1][0][0] = s_wgl_fixed_slp[it_tot]; /* Slope */
    }


  return(0);
}		/* end of MCMC_model2_it */



/*!
  \author Geir Storvik
  \brief Main part for age model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Simulating missing ages using the ::sample_ages_len_only routine
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating ages if errors in age-readings ::sample_ages_age_error
  - Simulating haul effects using the ::sample_age_haul routine
  - Simulating precision parameters for haul effects
  - Calculating likelihood

  All data are now of the Amigo-type structure.
*/
static int MCMC_it_age(int start_h,int i_force_acc,int i_len_only,
		       int i_it,int i_it_tot,int *o_acc_h)
{
  int err=0;
  int printFile;

  printFile = 0;  
  #ifdef DEBUG_PROG
  printFile = 1;
  #endif

  /* Sample missing ages */
  #ifdef DEBUG_PROG
  printf("Sample missing ages\n");
  #endif
  if(s_D_orig->coastal_cod)
    {
      if(i_it_tot==(s_burn_in+s_num_it_inner*s_num_it_outer-1))
	err = sample_ages_len_only_new_CC(s_D_orig,s_D_CC,s_age,s_D_age,s_length,s_length_CC,
					  s_D_lga,s_D_g_a,s_D_lga_CC,s_D_g_a_CC,printFile);
      else
	err = sample_ages_len_only_new_CC(s_D_orig,s_D_CC,s_age,s_D_age,s_length,s_length_CC,
					  s_D_lga,s_D_g_a,s_D_lga_CC,s_D_g_a_CC,0);
    }
  else if(s_COST)
    {
      if(i_it_tot==(s_burn_in+s_num_it_inner*s_num_it_outer-1))
	err = sample_ages_COST(s_D_orig,s_age,s_D_age,s_length,s_D_lga,s_D_g_a,printFile,i_it_tot);
      else
	err = sample_ages_COST(s_D_orig,s_age,s_D_age,s_length,s_D_lga,s_D_g_a,0,i_it_tot);
    }
  else
    {
      if(i_it_tot==(s_burn_in+s_num_it_inner*s_num_it_outer-1))
	{
	  //store simulated ages for last iteration
	  err = sample_ages_len_only_new(s_D_orig,s_age,s_D_age,s_length,
					 s_D_lga,s_D_g_a,printFile,i_it_tot);
	}
      else
	{
	  err = sample_ages_len_only_new(s_D_orig,s_age,s_D_age,s_length,
					 s_D_lga,s_D_g_a,0,i_it_tot);
	}
    }
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_ages\n");
      return(err);
    }

  /* Sample discarded fish - ages and lengths */
  if(s_length->cens_model) 
    {
      if(s_COST)
	{
          #ifdef DEBUG_PROG
	  printf("\nSampling censoring parameters\n");
          #endif
	  err = sample_cens_par(s_D_lga,s_D_orig,s_D_COST,i_it_tot);
	  if(err)
	    {
	      write_warning("MCMC_it_age:Error calling sample_cens_par\n");
	      return(err);
	    }
          #ifdef DEBUG_PROG
	  printf("\nSampling discarded fish\n");
          #endif
	  if(i_it==(s_num_it_inner*s_num_it_outer-1)) // print simulated discarded in last iteration
	    err = sample_discard_COST(s_age,s_D_age,s_length,s_D_lga,s_D_g_a,s_D_orig,s_D_COST,printFile,i_it_tot);
	  else
	    err = sample_discard_COST(s_age,s_D_age,s_length,s_D_lga,s_D_g_a,s_D_orig,s_D_COST,0,i_it_tot);
	  if(err)
	    {
	      write_warning("MCMC_it_age:Error calling sample_discarded_COST\n");
	      return(err);
	    }
	}
      else /* Sample discarded fish, old version not changed yet */
	{
	  printf("sample_discarded: must change to continuous age\n");
	}	    
    }  
  
  /* Calculate sufficient statistics */
  #ifdef DEBUG_PROG
  printf("Calculate sufficient statistics\n");
  #endif
  err = make_suff_age(s_D_age->glm->ncat,s_age,s_D_age,s_D_lga->haulweight,start_h);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_age:Error calling make_suff_age\n");
      #endif
      write_warning("MCMC_it_age:Error calling make_suff_age\n");
      return(err);
    }


  // Remove first haul effect from sampling
  /* Sample fixed and remaining random effects */
  #ifdef DEBUG_PROG
  printf("Sample fixed and random effects (not haul effect)\n");
  #endif
  s_D_age->glm->xcov[0]->n_cov--;
  err = sample_gauss_eff(s_age->gr_str_f,s_age->par,s_D_age->glm,start_h);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_gauss_eff\n");
      return(err);
    }
  // Add haul effect again
  s_D_age->glm->xcov[0]->n_cov++;


  #ifdef DEBUG_PROG
  int a,h,i,sum,suma;
  double *p, psum;

  p = CALLOC(s_D_age->glm->ncat,double);
  printf("MCMC_it_age:\n");
  psum = G_ZERO;
  for(a=0;a<s_D_age->glm->ncat;a++)
    {
      p[a] = exp(s_age->par->eff[a][0][0][0]);
      psum += p[a];
    }
  for(a=0;a<s_D_age->glm->ncat;a++)
    {
      printf("alpha_const=");
      for(i=0;i<s_D_age->glm->nxcov;i++)
	printf("%lf ",s_age->par->eff[a][i][0][0]);
      printf("p=%lf\n",p[a]/psum);
    }
  printf("sigma_haul=%lf\n",sqrt(G_ONE/s_age->par->tau_obs));

  sum = 0;
  printf("Num_age= ");
  for(a=0;a<s_D_age->glm->ncat;a++)
    {
      suma = 0;
      for(h=0;h<s_D_age->glm->nHaul;h++)
	suma += s_D_age->Ages[h][a];
      printf("%d ",suma);
      sum += suma;
    }
  printf("\n");
  printf("Total aged fish=%d\n",sum);
  FREE(p);
  #endif

  /* Sample complete alpha's */
  #ifdef DEBUG_PROG
  printf("Sample haul effect\n");
  #endif
  err = sample_age_alpha(s_age,s_D_age,start_h,i_force_acc,i_it,o_acc_h);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_age_alpha\n");
      return(err);
    }

  /* Sample precision parameters */
  #ifdef DEBUG_PROG
  printf("Sample precision parameters\n");
  #endif
  err = sample_precision_age_haul(start_h,s_age->par,s_age->alpha,s_D_age->glm);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_precision_haul\n");
      return(err);
    }
  s_D_age->glm->xcov[0]->n_cov--;
  err = sample_precision(start_h,s_age->par,s_D_age->glm);
  s_D_age->glm->xcov[0]->n_cov++;
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_precision\n");
      return(err);
    }

  /* Calculate likelihoods */
  err = calc_lik_age(s_age,s_D_age,&(s_age->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_age\n");
      return(err);
    }
  #ifdef DEBUG_PROG
    printf("Loglik_age=%lf\n",s_age->par->loglik);
    printf("tau_obs=%lf\n",s_age->par->tau_obs);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_age=%lf\n",s_age->par->loglik);
  #endif

  return(0);
}		/* end of MCMC_it_age */



/*!
  \author 
  \brief Main part for lga model inside each MCMC iteration if fixed lga model

  Simulations are performed by the following steps:
  - Calculating sufficient statistics for the lga model
  - Simulating lga parameters (Int, Slp, precision(fish)) using the ::sample_ routine 
  - Calculating likelihood
*/
static int MCMC_it_lga_fixed(int start_h,int i_it)
{
  int err;

  /* Calculate sufficient statistics */
  err = make_suff_lga(s_D_lga,s_D_g_a,start_h);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
      #endif
      write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
      return(err);
    }

  /* "sample" int, slp and tau_obs from a set of realisations */ 
  s_length->par->eff[0][0][0][0] = s_lga_fixed_int[i_it];     
  s_length->par->eff[0][1][0][0] = s_lga_fixed_slp[i_it];     
  s_length->par->tau_obs = s_lga_fixed_tau[i_it]; 

  if(s_D_g_a->g_a_model==1)
    {
      s_D_g_a->g_a_par[0] = s_lga_fixed_g_a_c[i_it];
      s_D_g_a->g_a_par[1] = s_lga_fixed_g_a_theta[i_it];
      s_D_g_a->g_a_par[2] = s_lga_fixed_g_a_gamma[i_it];
    }

  #ifdef DEBUG_PROG
  printf("MCMC_it_lga:\n");
  printf("Int=%lf\n",s_length->par->eff[0][0][0][0]);
  printf("Slope=%lf\n",s_length->par->eff[0][1][0][0]);
  printf("Precision=%lf\n",s_length->par->tau_obs);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Int=%lf\n",s_length->par->eff[0][0][0][0]);
  fprintf(g_caa_log,"Slope=%lf\n",s_length->par->eff[0][1][0][0]);
  fprintf(g_caa_log,"Precision=%lf\n",s_length->par->tau_obs);
  #endif

  err = calc_lik_lin(s_length,s_D_lga,&(s_length->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_lin\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("Loglik_lga=%lf\n",s_length->par->loglik);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_lga=%lf\n",s_length->par->loglik);
  #endif

  return(0);
}		/* end of MCMC_it_lga_fixed */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the age model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_age_fix is called first.
*/
double Qfunc_age_fix(int node, int nnode,char *arg)
{
  return(s_age->gr_str_f->Q[node][nnode]);
}		/* end of Qfunc_age_fix */



/*!
  \author Geir Storvik
  \brief Pick out an element of the the precision matrix for the age model

  For use when sampling using the GMRFLib library.
  For use if blockupdating fixed effects and precisions simultaneously. This
  is currently not implemented.
*/
double Qfunc_age_fix_new(int node, int nnode,char *arg)
{
  return(s_age->gr_str_f->Q_new[node][nnode]);
}		/* end of Qfunc_age_fix_new */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the random effects of the age model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_age_ran is called first. From an earlier version where the
  GMRFLib library was used for simulating haul effects in the age model. Now this is done
  in new routines.
*/
double Qfunc_age_ran(int node, int nnode,char *arg)
{
  return(s_age->par->tau_obs*(node==nnode));
}		/* end of Qfunc_age_ran */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the lga model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_length is called first.
*/
double Qfunc_length(int node, int nnode,char *arg)
{
  return(s_length->gr_str->Q[node][nnode]);
}		/* end of Qfunc_length */


double Qfunc_length_CC(int node, int nnode,char *arg)
{
  return(s_length_CC->gr_str->Q[node][nnode]);
}		/* end of Qfunc_length_CC */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the wgl model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_weight is called first.
*/
/*L:Qfunc_weight*
________________________________________________________________

		Qfunc_weight
________________________________________________________________

Name:		Qfunc_weight
Syntax:		
Description:    Find Q[i,j] for weight
                Requires that make_Q_weight is called first.
Side effects:   
Return value:   Q[i,j]
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
double Qfunc_weight(int node, int nnode,char *arg)
{
  return(s_weight->gr_str->Q[node][nnode]);
}		/* end of Qfunc_weight */
double Qfunc_weight_CC(int node, int nnode,char *arg)
{
  return(s_weight_CC->gr_str->Q[node][nnode]);
}		/* end of Qfunc_weight */

