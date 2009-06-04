/*!
  \file caa_predict.c
  \brief Containing the main routine for prection catch-at-age
  \author Geir Storvik

  The prediction is performed by the ::caa_predict routine.
*/
#include "caa.h"

static int initialize(int i_g_a_model);
static int re_initialize(int i_g_a_model);
static int find_catch_at_age(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     int *i_inc_haul,int i_nMC,int *i_season,
			     Data_totcatch *i_D_totcatch,
			     TC_struct *o_totcatch,double *i_par_haulsize,
			     double *o_mean_l,double *o_mean_w,int iter);
static int calculate_means_age(double **mu_age,int c,Age_struct *i_age,Data_age *i_D_age,
			       Data_totcatch *i_D_totcatch);
static int calculate_means_lga(double *mu_lga,double *o_var_lga,int c,LW_struct *i_length,
			       Data_lin *i_D_lga,int *i_inc_haul,Data_totcatch *i_D_totcatch);
static int calculate_means_wgl(double *mu_wgl,double *o_var_wgl,int c,LW_struct *i_weight,
			       Data_lin *i_D_wgl,int *i_inc_haul,Data_totcatch *i_D_totcatch);


/*(Non-global) External Variables*/

static int                 s_nMCMC;
static int                 s_burnin;
static int                 s_nMC;
static Age_struct         *s_age;
static LW_struct          *s_length;
static LW_struct          *s_weight;
static TC_struct          *s_totcatch;
static double             *s_mean_l;
static double             *s_mean_w;
static Data_age           *s_D_age;
static Data_g_a           *s_D_g_a;
static Data_lin           *s_D_lga;
static Data_lin           *s_D_wgl;
static int                *s_inc_haul;
static Data_totcatch      *s_D_totcatch;
static int      s_N_gauher;  /*!< Number of nodes in Gauss Hermite quadrature */
static double  *s_gauher_x;  /*!< Abscissas in Gauss Hermite quadrature */
static double  *s_gauher_w;  /*!< Weights in Gauss Hermite quadrature */

static int                 s_coastal_cod;
static LW_struct          *s_length_CC;
static LW_struct          *s_weight_CC;
static Data_g_a           *s_D_g_a_CC;
static Data_lin           *s_D_lga_CC;
static Data_lin           *s_D_wgl_CC;

static Data_COST          *s_D_COST;      /*!< Original data for COST project */
static int                 s_COST;        /*!< Indicator =1 if COST project, 0 otherwise */
static TC_struct          *s_totcatch_disc;   
static double             *s_mean_l_disc;
static double             *s_mean_w_disc;
static double             *s_planded;

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif


/*!
  \author Geir Storvik
  \brief Estimates catch-at-age for simulated parameters and random effects.

  Since this is a routine made to be called from R or Splus, all input 
  variables are pointers to vectors. The routine therefore starts to convert 
  input data into approperiate c-structures (as defined in caa.h).

  Thereafter a loop through all simulations is performed with the main
  calculations performed in the ::find_catch_at_age routine.

  The predictions are stored in a int vector, o_mcmc_totcatch with all
  predictions from one iteration in one sequential block. See the
  ::write_it_totcatch routine for the format of this block.

  In addition mean of length-given-age and mean of weight-given-age for each
  simulation is stored in o_mcmc_mean_l and o_mcmc_mean_w. Again the
  ::write_it_totcatch routine described the format.
*/
SEXP caa_predict(SEXP i_mcmc_samp,SEXP i_common_par,
		 SEXP i_data_age,SEXP i_data_lga,SEXP i_data_wgl,SEXP i_data_catch,
		 SEXP i_par_haulsize,SEXP i_dist_cell,
		 SEXP i_N_l_int,SEXP i_l_int,SEXP i_nMC,SEXP i_data_COST)
{
  int    i,it,ind,c;
  int    n_totcatch,n_mean;

  int time_now, time_start;
  time_start = clock();
  
  #ifdef LOG_FILE
  g_caa_log = fopen("caa_logfile_predict.txt","w");
  #endif
   
  SEXP       elmt = R_NilValue;
  int        n_protect=0;
  int        nOutputVar;
 
  /* Variables connect to input data */
  double *mcmc1=NULL, *mcmc2=NULL, *mcmc_COST=NULL;
  int *num_par1, *num_par2, num_par_COST;
  int  age_nBoats, nAges, *a_vec, *season;
  int *n_cov, *ispat, *icell;
  int *age_int_nFac, *age_int_fix, *age_int_c_cov;
  int *age_hsz_nFac, *age_hsz_fix, *age_hsz_c_cov;
  int  lga_nBoats;
  int *lga_int_nFac, *lga_int_fix, *lga_int_c_cov;
  int *lga_slp_nFac, *lga_slp_fix, *lga_slp_c_cov;
  int *lga_hsz_nFac, *lga_hsz_fix, *lga_hsz_c_cov;
  int  lga_g_a_model, lga_g_a_ncat, lga_g_a_nSeason, lga_cens_model;
  int *lga_g_a_a2Age_vec;
  double *lga_g_a_avec;
  double *lga_cens;
  int  wgl_nBoats;
  int *wgl_int_nFac, *wgl_int_fix, *wgl_int_c_cov;
  int *wgl_slp_nFac, *wgl_slp_fix, *wgl_slp_c_cov;
  int *wgl_hsz_nFac, *wgl_hsz_fix, *wgl_hsz_c_cov;
  int *num_adj_area, *adj_area;
  int  tot_nCell, tot_nFactors;
  int *tot_fac_age_int, *tot_fac_age_hsz;
  int *tot_fac_lga_int, *tot_fac_lga_slp, *tot_fac_lga_hsz;
  int *tot_fac_wgl_int, *tot_fac_wgl_slp, *tot_fac_wgl_hsz;
  int *tot_factors;
  double *tot_catch=NULL;
  double *par_haulsize=NULL;
  int N_l_int,  n_MC;
  double *l_int;

  /* Output variables */
  SEXP     mcmc_totcatch, mcmc_mean_l, mcmc_mean_w, err;
  SEXP     mcmc_totcatch_disc, mcmc_mean_l_disc, mcmc_mean_w_disc, mcmc_planded;
  int     *o_err;
  double  *o_mcmc_totcatch, *o_mcmc_mean_l, *o_mcmc_mean_w;
  double  *o_mcmc_totcatch_disc, *o_mcmc_mean_l_disc, *o_mcmc_mean_w_disc, *o_mcmc_planded;
  SEXP     resList, resList_names;
  char    *res_names[4] = {"mcmc_totcatch", "mcmc_mean_l", "mcmc_mean_w", "err"};
  char    *res_names_COST[8] = {"mcmc_totcatch_land", "mcmc_mean_l_land", "mcmc_mean_w_land", "mcmc_totcatch_disc", "mcmc_mean_l_disc", "mcmc_mean_w_disc", "mcmc_planded", "err"};

  /* Connecting input data from R objects to variables */
  /* NOTE: Variables are read-only! */

  s_nMCMC = INTEGER_VALUE(getListElement(i_mcmc_samp,"nMCMC"));
  s_burnin = INTEGER_VALUE(getListElement(i_mcmc_samp,"burnin"));
  mcmc1 = NUMERIC_POINTER(getListElement(i_mcmc_samp,"samples1"));
  mcmc2 = NUMERIC_POINTER(getListElement(i_mcmc_samp,"samples2"));
  num_par1 = INTEGER_POINTER(AS_INTEGER(getListElement(i_mcmc_samp,"numpar1")));
  num_par2 = INTEGER_POINTER(AS_INTEGER(getListElement(i_mcmc_samp,"numpar2")));

  n_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"ncov")));
  ispat = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"ispat")));
  icell = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"icell")));
  num_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"neigh.num")));
  adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"neigh.adj")));
  s_inc_haul = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"inchaul")));

  age_nBoats = INTEGER_VALUE(getListElement(i_data_age,"nBoats"));
  nAges = INTEGER_VALUE(getListElement(i_data_age,"nAges"));
  a_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"avec")));
  season = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"season")));
  age_int_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"int.nfac")));
  age_int_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"int.fix")));
  age_int_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"int.cov")));
  age_hsz_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"hsz.nfac"))); 
  age_hsz_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"hsz.fix")));
  age_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"hsz.cov")));

  lga_nBoats = INTEGER_VALUE(getListElement(i_data_lga,"nBoats"));
  lga_int_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"int.nfac")));
  lga_int_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"int.fix")));
  lga_int_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"int.cov")));
  lga_slp_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"slp.nfac")));
  lga_slp_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"slp.fix")));
  lga_slp_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"slp.cov")));
  lga_hsz_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"hsz.nfac"))); 
  lga_hsz_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"hsz.fix")));
  lga_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"hsz.cov")));
  lga_g_a_model = INTEGER_VALUE(getListElement(i_data_lga,"lgarelInt"));
  lga_g_a_ncat = INTEGER_VALUE(getListElement(i_data_lga,"ga.ncat"));
  lga_g_a_nSeason = INTEGER_VALUE(getListElement(i_data_lga,"ga.nSeason"));
  lga_g_a_a2Age_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_lga,"ga.a2Age.vec")));
  lga_g_a_avec = NUMERIC_POINTER(getListElement(i_data_lga,"ga.avec"));
  lga_cens_model = INTEGER_VALUE(getListElement(i_data_lga,"lga.cens.model"));
  lga_cens = NUMERIC_POINTER(getListElement(i_data_lga,"lga.cens"));
  
  wgl_nBoats = INTEGER_VALUE(getListElement(i_data_wgl,"nBoats"));
  wgl_int_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"int.nfac")));
  wgl_int_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"int.fix")));
  wgl_int_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"int.cov")));
  wgl_slp_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"slp.nfac")));
  wgl_slp_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"slp.fix")));
  wgl_slp_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"slp.cov")));
  wgl_hsz_nFac = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"hsz.nfac"))); 
  wgl_hsz_fix = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"hsz.fix")));
  wgl_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_wgl,"hsz.cov")));
  
  tot_nCell = INTEGER_VALUE(getListElement(i_data_catch,"nCell"));
  tot_nFactors = INTEGER_VALUE(getListElement(i_data_catch,"nfactors"));
  tot_fac_age_int = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.age.int")));
  tot_fac_age_hsz = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.age.hsz")));
  tot_fac_lga_int = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.lga.int")));
  tot_fac_lga_slp = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.lga.slp")));
  tot_fac_lga_hsz = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.lga.hsz")));
  tot_fac_wgl_int = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.wgl.int")));
  tot_fac_wgl_slp = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.wgl.slp")));
  tot_fac_wgl_hsz = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"fac.wgl.hsz")));
  tot_factors = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_catch,"factors")));
  tot_catch =  NUMERIC_POINTER(getListElement(i_data_catch,"catch"));

  par_haulsize = NUMERIC_POINTER(i_par_haulsize);
  N_l_int = INTEGER_VALUE(i_N_l_int);

  s_nMC = INTEGER_VALUE(i_nMC);

  s_coastal_cod = INTEGER_VALUE(getListElement(i_common_par, "coastal.cod"));

  s_COST = INTEGER_VALUE(getListElement(i_data_COST, "COST"));
  if(s_COST) 
    {
      num_par_COST = INTEGER_VALUE(getListElement(i_data_COST, "numpar"));
      mcmc_COST = NUMERIC_POINTER(getListElement(i_data_COST,"samples"));
    }

  /* Allocate space for output */
  n_totcatch = (int) nAges*N_l_int*(s_nMCMC-s_burnin);
  PROTECT(mcmc_totcatch = NEW_NUMERIC(n_totcatch));
  n_protect++;
  o_mcmc_totcatch = NUMERIC_POINTER(mcmc_totcatch);
  n_mean = nAges*(s_nMCMC-s_burnin);
  PROTECT(mcmc_mean_l = NEW_NUMERIC(n_mean));
  n_protect++;
  o_mcmc_mean_l = NUMERIC_POINTER(mcmc_mean_l);
  PROTECT(mcmc_mean_w = NEW_NUMERIC(n_mean));
  n_protect++;
  o_mcmc_mean_w = NUMERIC_POINTER(mcmc_mean_w);
  PROTECT(err = NEW_INTEGER(1));
  n_protect++;
  o_err = INTEGER_POINTER(err);

  if(s_COST)
    {
      PROTECT(mcmc_totcatch_disc = NEW_NUMERIC(n_totcatch));
      n_protect++;
      o_mcmc_totcatch_disc = NUMERIC_POINTER(mcmc_totcatch_disc);
      PROTECT(mcmc_mean_l_disc = NEW_NUMERIC(n_mean));
      n_protect++;
      o_mcmc_mean_l_disc = NUMERIC_POINTER(mcmc_mean_l_disc);
      PROTECT(mcmc_mean_w_disc = NEW_NUMERIC(n_mean));
      n_protect++;
      o_mcmc_mean_w_disc = NUMERIC_POINTER(mcmc_mean_w_disc);
      PROTECT(mcmc_planded = NEW_NUMERIC(tot_nCell*(s_nMCMC-s_burnin)));
      n_protect++;
      o_mcmc_planded = NUMERIC_POINTER(mcmc_planded);
    }

  if(s_COST==0)
    {
      /* Make "names" attribute for output list */ 
      nOutputVar = 4;  
      PROTECT(resList_names = allocVector(STRSXP, nOutputVar));
      n_protect++;
      for(i = 0; i < nOutputVar; i++)
	SET_STRING_ELT(resList_names, i,  mkChar(res_names[i]));
      
      /* Make list with output variables */ 
      PROTECT(resList = allocVector(VECSXP, nOutputVar)); 
      n_protect++;
      SET_VECTOR_ELT(resList, 0, mcmc_totcatch);   
      SET_VECTOR_ELT(resList, 1, mcmc_mean_l);
      SET_VECTOR_ELT(resList, 2, mcmc_mean_w);
      SET_VECTOR_ELT(resList, 3, err);
      setAttrib(resList, R_NamesSymbol, resList_names); //attaching the vector names
    }
  else
    {
      nOutputVar = 8;
      PROTECT(resList_names = allocVector(STRSXP, nOutputVar));
      n_protect++;
      for(i = 0; i < nOutputVar; i++)
	SET_STRING_ELT(resList_names, i,  mkChar(res_names_COST[i]));
      PROTECT(resList = allocVector(VECSXP, nOutputVar)); 
      n_protect++;
      SET_VECTOR_ELT(resList, 0, mcmc_totcatch);   
      SET_VECTOR_ELT(resList, 1, mcmc_mean_l);
      SET_VECTOR_ELT(resList, 2, mcmc_mean_w);
      SET_VECTOR_ELT(resList, 3, mcmc_totcatch_disc);   
      SET_VECTOR_ELT(resList, 4, mcmc_mean_l_disc);
      SET_VECTOR_ELT(resList, 5, mcmc_mean_w_disc);
      SET_VECTOR_ELT(resList, 6, mcmc_planded);
      SET_VECTOR_ELT(resList, 7, err);
      setAttrib(resList, R_NamesSymbol, resList_names); //attaching the vector names

    }

  #ifdef DEBUG_INPUT
  *o_err = write_input_predict_new(i_mcmc_samp,i_common_par,
				   i_data_age,i_data_lga,i_data_wgl,
				   i_data_catch,
				   i_par_haulsize,i_dist_cell,
				   i_N_l_int,i_l_int,i_nMC);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling write_input_predict\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  #endif


  #ifdef LOG_FILE
   fprintf(g_caa_log,"Initializing predict\n");
  #endif

   /* Calculate Gauss hermite abscissas and weights */
  s_N_gauher = 30;
  s_gauher_x = CALLOC(s_N_gauher+1,double);  //Free ok
  s_gauher_w = CALLOC(s_N_gauher+1,double);  //Free ok
  gauher(s_gauher_x,s_gauher_w,s_N_gauher);

  /* Make age data */
  #ifdef DEBUG_PREDICT
  printf("Make age data: makedata_age1\n");
  #endif
  *o_err = makedata_age1(age_nBoats,nAges,a_vec,
			 n_cov[0],age_int_nFac,ispat[0],
			 age_int_fix,age_int_c_cov,
			 num_adj_area,adj_area,
			 n_cov[1],age_hsz_nFac,ispat[1],
			 age_hsz_fix,age_hsz_c_cov,
			 num_adj_area,adj_area,&s_D_age);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling makedata_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Make g_a data */
  #ifdef DEBUG_PREDICT
  printf("Make g_a data: makedata_g_a\n");
  #endif
  *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			lga_g_a_model,s_D_age,&s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      #ifdef DEBUG_PREDICT
      printf("Make extra g_a data: makedata_g_a\n");
      #endif
      *o_err = makedata_g_a(lga_g_a_ncat,lga_g_a_nSeason,lga_g_a_avec,lga_g_a_a2Age_vec,
			    lga_g_a_model,s_D_age,&s_D_g_a_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_g_a\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  //*o_err = make_cell_constr_age(s_D_age,icell);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling make_cell_constr_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Make lga data */
  #ifdef DEBUG_PREDICT
  printf("Make lga data: makedata_lin1\n");
  #endif
  *o_err = makedata_lin1(lga_nBoats,
			 n_cov[2],lga_int_nFac,ispat[2],
			 lga_int_fix,lga_int_c_cov,
			 num_adj_area,adj_area,
			 n_cov[3],lga_slp_nFac,ispat[3],
			 lga_slp_fix,lga_slp_c_cov,
			 num_adj_area,adj_area,
			 n_cov[4],lga_hsz_nFac,ispat[4],lga_hsz_fix,lga_hsz_c_cov,
			 num_adj_area,adj_area,
			 &s_D_lga);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling makedata_lin1 for lga\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      #ifdef DEBUG_PREDICT
      printf("Make extra lga data: makedata_lin1\n");
      #endif
      *o_err = makedata_lin1(lga_nBoats,
			     n_cov[8],lga_int_nFac,ispat[8],
			     lga_int_fix,lga_int_c_cov,
			     num_adj_area,adj_area,
			     n_cov[9],lga_slp_nFac,ispat[9],
			     lga_slp_fix,lga_slp_c_cov,
			     num_adj_area,adj_area,
			     n_cov[10],lga_hsz_nFac,ispat[10],lga_hsz_fix,lga_hsz_c_cov,
			     num_adj_area,adj_area,
			     &s_D_lga_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_lin1 for lga\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  //*o_err = make_cell_constr_lin(s_D_lga,&icell[2]);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling make_cell_constr_lin for lga\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  /* Read weight given length data */
  #ifdef DEBUG_PREDICT
  printf("Make wgl data: makedata_lin1\n");
  #endif
  *o_err = makedata_lin1(wgl_nBoats,
			 n_cov[5],wgl_int_nFac,ispat[5],
			 wgl_int_fix,wgl_int_c_cov,
			 num_adj_area,adj_area,
			 n_cov[6],wgl_slp_nFac,ispat[6],
			 wgl_slp_fix,wgl_slp_c_cov,
			 num_adj_area,adj_area,
			 n_cov[7],wgl_hsz_nFac,ispat[7],
			 wgl_hsz_fix,wgl_hsz_c_cov,
			 num_adj_area,adj_area,
			 &s_D_wgl);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling makedata_lin1 for wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      #ifdef DEBUG_PREDICT
      printf("Make extra wgl data: makedata_lin1\n");
      #endif
      *o_err = makedata_lin1(wgl_nBoats,
			     n_cov[11],wgl_int_nFac,ispat[11],
			     wgl_int_fix,wgl_int_c_cov,
			     num_adj_area,adj_area,
			     n_cov[12],wgl_slp_nFac,ispat[12],
			     wgl_slp_fix,wgl_slp_c_cov,
			     num_adj_area,adj_area,
			     n_cov[13],wgl_hsz_nFac,ispat[13],
			     wgl_hsz_fix,wgl_hsz_c_cov,
			     num_adj_area,adj_area,
			     &s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_lin1 for wgl\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  //*o_err = make_cell_constr_lin(s_D_wgl,&icell[5]);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling make_cell_constr_lin for wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  //Remove haul effects
  if(0)
    {
      s_D_age->glm->xcov[0]->n_cov -= s_inc_haul[0];
      if(n_cov[1]>0)
	s_D_age->glm->xcov[1]->n_cov -= s_inc_haul[1];
      s_D_lga->glm->xcov[0]->n_cov -= s_inc_haul[2];
      s_D_lga->glm->xcov[1]->n_cov -= s_inc_haul[3];
      if(n_cov[4]>0)
	s_D_lga->glm->xcov[2]->n_cov -= s_inc_haul[4];
      s_D_wgl->glm->xcov[0]->n_cov -= s_inc_haul[5];
      s_D_wgl->glm->xcov[1]->n_cov -= s_inc_haul[6];
      if(n_cov[7]>0)
	s_D_wgl->glm->xcov[2]->n_cov -= s_inc_haul[7];
    }

  #ifdef DEBUG_PREDICT
  printf("Initialize model\n");
  #endif
  *o_err = initialize(lga_g_a_model);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  // Cell distributions
  if(s_coastal_cod==0)
    {
      *o_err = makedata_cell_dist(i_dist_cell,icell,s_age,s_D_age,s_length,s_D_lga,s_weight,s_D_wgl);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_cell_dist\n");
      #ifdef LOG_FILE
	  fclose(g_caa_log);
      #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  else
    {
      *o_err = makedata_cell_dist_CC(i_dist_cell,icell,s_age,s_D_age,
				     s_length,s_D_lga,s_weight,s_D_wgl,
				     s_length_CC,s_D_lga_CC,s_weight_CC,s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_cell_dist\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  #ifdef DEBUG_PREDICT
  printf("Make data totcatch\n");
  #endif
  *o_err = makedata_totcatch(s_D_age,s_D_lga,s_D_wgl,
			     n_cov,tot_nCell,tot_nFactors,s_inc_haul,
			     tot_fac_age_int,tot_fac_age_hsz,
			     tot_fac_lga_int,tot_fac_lga_slp,tot_fac_lga_hsz,
			     tot_fac_wgl_int,tot_fac_wgl_slp,tot_fac_wgl_hsz,
			     tot_factors,tot_catch,&s_D_totcatch);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling makedata_totcatch\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
     return(resList);
    }
  s_totcatch = CALLOC(1,TC_struct);     // Free ok
  s_D_totcatch->nlint = N_l_int;
  s_D_totcatch->l_int = NUMERIC_POINTER(i_l_int);
  s_totcatch->catch_at_age = 
     Mmatrix_2d(0,s_D_age->glm->ncat-1,
		0,s_D_totcatch->nlint-1,sizeof(double),1); // Free ok
  s_totcatch->mcmc = o_mcmc_totcatch;
  s_totcatch->mean_l = o_mcmc_mean_l;
  s_totcatch->mean_w = o_mcmc_mean_w;  

  if(s_COST)
    {
      s_totcatch_disc = CALLOC(1,TC_struct);     // Free ok
      s_totcatch_disc->catch_at_age = 
	Mmatrix_2d(0,s_D_age->glm->ncat-1,
		   0,s_D_totcatch->nlint-1,sizeof(double),1); // Free ok
      s_totcatch_disc->mcmc = o_mcmc_totcatch_disc;
      s_totcatch_disc->mean_l = o_mcmc_mean_l_disc;
      s_totcatch_disc->mean_w = o_mcmc_mean_w_disc;  
      #ifdef DEBUG_PREDICT
      printf("Make data COST\n");
      #endif
      *o_err = makedata_COST_predict(s_D_lga->glm->nHaul,&s_D_COST);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling makedata_COST_predict\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      s_planded = CALLOC(tot_nCell,double);
    }

  s_age->par->num_var = num_par1[0];
  s_age->par->mcmc = &mcmc1[0];

  s_length->par->num_var = num_par1[1];
  s_length->par->mcmc = &mcmc1[num_par1[0]*s_nMCMC];

  s_D_g_a->g_a_mcmc = &mcmc1[(num_par1[0]+num_par1[1])*s_nMCMC];

  if(s_coastal_cod)
    {
      s_length_CC->par->num_var = num_par1[3];
      s_length_CC->par->mcmc = &mcmc1[(num_par1[0]+num_par1[1]+num_par1[2])*s_nMCMC];
      s_D_g_a_CC->g_a_mcmc = &mcmc1[(num_par1[0]+num_par1[1]+num_par1[2]+num_par1[3])*s_nMCMC];
    }
  if(s_COST)
    {
      s_D_COST->num_var = num_par_COST;
      s_D_COST->mcmc = &mcmc_COST[0];
    }

  s_length->cens_model = lga_cens_model;
  if(s_length->cens_model)
    {
      s_length->cens_k = lga_cens[0];
      s_length->cens_m = lga_cens[1];
      s_length->cens_r = lga_cens[2];
      s_length->cens_Nlim = lga_cens[3];
    }

  s_weight->par->num_var = num_par2[0];
  s_weight->par->mcmc = &mcmc2[0];

  if(s_coastal_cod)
    {
      s_weight_CC->par->num_var = num_par2[1];
      s_weight_CC->par->mcmc = &mcmc2[num_par2[0]*s_nMCMC];
    }

  time_now = clock();
  //printf("CPU time used %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"CPU time used %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Perform prediction\n");
  #endif

  #ifdef DEBUG_PREDICT
  printf("Start MCMC iterations\n");
  #endif

  for(it=s_burnin;it<s_nMCMC;it++)
    {
      #ifdef DEBUG_PREDICT
      printf("Read mcmc samples\n");
      #endif
      s_D_age->glm->xcov[0]->n_cov--;
      *o_err= read_it(it,s_D_age->glm,s_age->par);
      if(*o_err)
        {
          write_warning("caa_predict:Error calling read_it for age model\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      s_D_age->glm->xcov[0]->n_cov++;
      *o_err = read_it(it,s_D_lga->glm,s_length->par);
      if(*o_err)
        {
          write_warning("caa_predict:Error calling read_it for lga model\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = read_it(it,s_D_lga_CC->glm,s_length_CC->par);
	  if(*o_err)
	    {
	      write_warning("caa_predict:Error calling read_it for lga model\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
      if(s_COST)
	{
	  *o_err = read_it_COST(it,s_D_COST);
	  if(*o_err)
	    {
	      write_warning("caa_predict:Error calling read_it_COST\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
      *o_err = read_it(it,s_D_wgl->glm,s_weight->par);
      if(*o_err)
        {
          write_warning("caa_predict:Error calling read_it for wgl model\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = read_it(it,s_D_wgl_CC->glm,s_weight_CC->par);
	  if(*o_err)
	    {
	      write_warning("caa_predict:Error calling read_it for wgl model\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
      *o_err = read_it_g_a(it,s_D_g_a);
      if(*o_err)
        {
          write_warning("caa_predict:Error calling read_it_g_a\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_coastal_cod)
	{
	  *o_err = read_it_g_a(it,s_D_g_a_CC);
	  if(*o_err)
	    {
	      write_warning("caa_predict:Error calling read_it_g_a\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	}
      #ifdef DEBUG_PREDICT
      printf("find_catch_at_age\n");
      #endif
      time_start = clock();
      *o_err = find_catch_at_age(s_age,s_D_age,s_length,s_D_lga,s_weight,s_D_wgl,s_inc_haul,s_nMC,
				 season,s_D_totcatch,s_totcatch,par_haulsize,s_mean_l,s_mean_w,it);
      time_now = clock();
      //printf("CPU time used in find_catch_at_age in iteration %d: %fs\n", it,
      //	      (double)(time_now-time_start)/CLOCKS_PER_SEC);
      #ifdef LOG_FILE
      time_now = clock();
      fprintf(g_caa_log,"CPU time used in find_catch_at_age in iteration %d: %fs\n", it,
	      (double)(time_now-time_start)/CLOCKS_PER_SEC);
      #endif
      if(*o_err)
        {
          write_warning("caa_predict:Error calling find_catch_at_age\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      #ifdef DEBUG_PREDICT
      printf("write samples\n");
      #endif
      *o_err = write_it_totcatch(it-s_burnin,s_D_totcatch->nCell,
				 s_D_age->glm->ncat,
				 s_D_totcatch->nlint,s_totcatch,
				 s_mean_l,s_mean_w);
      if(*o_err)
        {
          write_warning("caa_predict:Error calling write_totcatch\n");
          #ifdef LOG_FILE
          fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
      if(s_COST)
	{
	  *o_err = write_it_totcatch(it-s_burnin,s_D_totcatch->nCell,
				     s_D_age->glm->ncat,
				     s_D_totcatch->nlint,s_totcatch_disc,
				     s_mean_l_disc,s_mean_w_disc);
	  if(*o_err)
	    {
	      write_warning("caa_predict:Error calling write_totcatch\n");
              #ifdef LOG_FILE
	      fclose(g_caa_log);
              #endif
	      UNPROTECT(n_protect);
	      return(resList);
	    }
	  ind = (it-s_burnin)*s_D_totcatch->nCell;
	  for(c=0; c<s_D_totcatch->nCell; c++)
	    o_mcmc_planded[ind+c] = s_planded[c];
	}
    }
   #ifdef LOG_FILE
   fprintf(g_caa_log,"Cleaning up\n");
   #endif
   #ifdef DEBUG_PREDICT
   printf("Cleaning up\n");
   #endif

  // Clean up
  Fmatrix_2d(&s_totcatch->catch_at_age[0][0],&s_totcatch->catch_at_age[0]);
  FREE(s_totcatch);
  FREE(s_gauher_x);
  FREE(s_gauher_w);

  if(s_COST)
    {
      Fmatrix_2d(&s_totcatch_disc->catch_at_age[0][0],&s_totcatch_disc->catch_at_age[0]);
      FREE(s_totcatch_disc);
      FREE(s_planded);
    }

  if(s_coastal_cod==0)
    {
      *o_err = re_makedata_cell_dist(i_dist_cell,icell,s_age,s_D_age,
				     s_length,s_D_lga,s_weight,s_D_wgl);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling re_makedata_cell_dist\n");
      #ifdef LOG_FILE
	  fclose(g_caa_log);
      #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }
  else
    {
      *o_err = re_makedata_cell_dist_CC(i_dist_cell,icell,s_age,s_D_age,
					s_length,s_D_lga,s_weight,s_D_wgl,
					s_length_CC,s_D_lga_CC,s_weight_CC,s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling re_makedata_cell_dist\n");
      #ifdef LOG_FILE
	  fclose(g_caa_log);
      #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  *o_err = re_initialize(lga_g_a_model);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling re_initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  *o_err = re_makedata_totcatch(s_D_age,s_D_lga,s_D_wgl,
				tot_nCell,tot_nFactors,
				tot_fac_age_int,tot_fac_age_hsz,
				tot_fac_lga_int,tot_fac_lga_slp,tot_fac_lga_hsz,
				tot_fac_wgl_int,tot_fac_wgl_slp,tot_fac_wgl_hsz,
				tot_factors,tot_catch,&s_D_totcatch);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling re_makedata_totcatch\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  *o_err = re_makedata_lin1(wgl_nBoats,
			    n_cov[5],wgl_int_nFac,ispat[5],wgl_int_fix,wgl_int_c_cov,
			    num_adj_area,adj_area,
			    n_cov[6],wgl_slp_nFac,ispat[6],wgl_slp_fix,wgl_slp_c_cov,
			    num_adj_area,adj_area,
			    n_cov[7],wgl_hsz_nFac,ispat[7],
			    wgl_hsz_fix,wgl_hsz_c_cov,
			    num_adj_area,adj_area,
			    &s_D_wgl);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      *o_err = re_makedata_lin1(wgl_nBoats,
				n_cov[5],wgl_int_nFac,ispat[5],wgl_int_fix,wgl_int_c_cov,
				num_adj_area,adj_area,
				n_cov[6],wgl_slp_nFac,ispat[6],wgl_slp_fix,wgl_slp_c_cov,
				num_adj_area,adj_area,
				n_cov[7],wgl_hsz_nFac,ispat[7],
				wgl_hsz_fix,wgl_hsz_c_cov,
				num_adj_area,adj_area,
				&s_D_wgl_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling re_makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  *o_err = re_makedata_lin1(lga_nBoats,
			    n_cov[2],lga_int_nFac,ispat[2],lga_int_fix,lga_int_c_cov,
			    num_adj_area,adj_area,
			    n_cov[3],lga_slp_nFac,ispat[3],lga_slp_fix,lga_slp_c_cov,
			    num_adj_area,adj_area,
			    n_cov[4],lga_hsz_nFac,ispat[4],lga_hsz_fix,lga_hsz_c_cov,
			    num_adj_area,adj_area,
			    &s_D_lga);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }
  if(s_coastal_cod)
    {
      *o_err = re_makedata_lin1(lga_nBoats,
				n_cov[2],lga_int_nFac,ispat[2],lga_int_fix,lga_int_c_cov,
				num_adj_area,adj_area,
				n_cov[3],lga_slp_nFac,ispat[3],lga_slp_fix,lga_slp_c_cov,
				num_adj_area,adj_area,
				n_cov[4],lga_hsz_nFac,ispat[4],lga_hsz_fix,lga_hsz_c_cov,
				num_adj_area,adj_area,
				&s_D_lga_CC);
      if(*o_err)
	{
	  write_warning("caa_predict:Error calling re_makedata_lin1\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  UNPROTECT(n_protect);
	  return(resList);
	}
    }

  *o_err = re_makedata_age1(age_nBoats,nAges,a_vec,
			    n_cov[0],age_int_nFac,ispat[0],
			    age_int_fix,age_int_c_cov,
			    num_adj_area,adj_area,
			    n_cov[1],age_hsz_nFac,ispat[1],
			    age_hsz_fix,age_hsz_c_cov,
			    num_adj_area,adj_area,&s_D_age);
  if(*o_err)
    {
      write_warning("caa_predict:Error calling re_makedata_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      UNPROTECT(n_protect);
      return(resList);
    }

  #ifdef LOG_FILE
  fclose(g_caa_log);
  #endif

  UNPROTECT(n_protect);

  return(resList);    /* End of caa_predict */
}


/*!
  \author Geir Storvik
  \brief Memory allocation and initial values for the sub-models
*/
static int initialize(int i_g_a_model)
{
  int err;

  /* Allocating memory for age structs */
  err = alloc_age(s_D_age,&s_age);
  if(err)
    {
      write_warning("initialize:Error calling alloc_age\n");
      return(err);
    }

  /* Allocating memory for length structs */
  err = alloc_lin(s_D_lga,&s_length);
  if(err)
    {
      write_warning("initialize:Error calling alloc_lin\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = alloc_lin(s_D_lga_CC,&s_length_CC);
      if(err)
	{
	  write_warning("initialize:Error calling alloc_lin\n");
	  return(err);
	}
    }

  /* Allocating memory for weight structs */
  err = alloc_lin(s_D_wgl,&s_weight);
  if(err)
    {
      write_warning("initialize:Error calling init_lin\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = alloc_lin(s_D_wgl_CC,&s_weight_CC);
      if(err)
	{
	  write_warning("initialize:Error calling init_lin\n");
	  return(err);
	}
    }

  if(s_D_g_a->g_a_model>0)
    {
      err = sample_g_a_initialize(s_D_g_a->ncat,s_D_g_a->g_a_model);
      if(err)
	{
	  write_warning("initialize:Error calling sampling_g_a_initialize\n");
	  return(err);
	}
    }

  s_mean_l = CALLOC(s_D_age->glm->ncat,double);
  s_mean_w = CALLOC(s_D_age->glm->ncat,double);
  if(s_COST)
    {
      s_mean_l_disc = CALLOC(s_D_age->glm->ncat,double);
      s_mean_w_disc = CALLOC(s_D_age->glm->ncat,double);
    }

  return(0);
}            /* End of initialize */



/*!
  \brief Reallocate memory allocated in initialize
  \author Geir Storvik
*/
static int re_initialize(int i_g_a_model)
{
  int err;

  if(s_D_g_a->g_a_model>0)
    {
      err = sample_g_a_re_initialize();
      if(err)
	{
	  write_warning("re_initialize:Error calling sampling_g_a_re_initialize\n");
	  return(err);
	}
    }
  /* Allocating memory for age structs */
  err = re_alloc_age(s_D_age,&s_age);
  if(err)
    {
      write_warning("re_initialize:Error calling re_alloc_age\n");
      return(err);
    }

  /* Allocating memory for length structs */
  err = re_alloc_lin(s_D_lga,&s_length);
  if(err)
    {
      write_warning("re_initialize:Error calling re_alloc_lin\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = re_alloc_lin(s_D_lga_CC,&s_length_CC);
      if(err)
	{
	  write_warning("re_initialize:Error calling re_alloc_lin\n");
	  return(err);
	}
    }

  /* Allocating memory for weight structs */
  err = re_alloc_lin(s_D_wgl,&s_weight);
  if(err)
    {
      write_warning("re_initialize:Error calling re_init_lin\n");
      return(err);
    }
  if(s_coastal_cod)
    {
      err = re_alloc_lin(s_D_wgl_CC,&s_weight_CC);
      if(err)
	{
	  write_warning("re_initialize:Error calling re_init_lin\n");
	  return(err);
	}
    }

  FREE(s_mean_l);
  FREE(s_mean_w);
  if(s_COST)
    {
      FREE(s_mean_l);
      FREE(s_mean_w);
    }

  return(0);
}            /* End of re_initialize */



/*!
  \author Geir Storvik
  \brief Estimates catch-at-age for current parameters and random effects.

  The catches are assumed to be sampled from other hauls than those used in
  the fitting. This means that the haul-effects are unknown. For the lga and
  wgl models this effect can be integrated out. For the age model, however,
  the haul-effects are included through Monte Carlo estimation of the
  age-probablities. Using this both expected age-probablities and expected
  weights can be calculated from which catch-at-age can be predicted.

  For the age model the haul effect is assumed always included and is used in the
  Monte Carlo estimate of the age probabilities.
  For the lga and wgl models, haul effects in the intercept are merged together with
  observation effects in the estimation of the conditional means.
*/

static int find_catch_at_age(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     int *i_inc_haul,int i_nMC,int *i_season,
			     Data_totcatch *i_D_totcatch,
			     TC_struct *o_totcatch,double *i_par_haulsize,
			     double *o_mean_l,double *o_mean_w,int iter)
{
  int     a,c,h,i,j,k,l,Hsz,err;
  int     a_start,a_end,season;
  double  mu_lga[3],mu_wgl[3],T,sum_T;
  double  mu_lga_CC[3],mu_wgl_CC[3],mu_lga_hsz_CC,mu_wgl_hsz_CC,var_lga_CC,var_wgl_CC;
  double  mu_hsz,sd_hsz,hsz,mu_lga_hsz,mu_wgl_hsz;
  double  eff,pa_sum,A,B,mean_w,mean_l,E_a,E_a_tot;
  double  var_lga,var_wgl,var,sd;
  double  x,cens_k,cens_m,cens_r;
  Data_cov *xcov;
  double  *E_lga,*E_wga,**mu_age, *pa, *E_pa, **E_pl_a, sum;
  double  sig,arg;
  double  *mu,*int_num_length,*int_num_weight,*E_pa_land,*E_pa_disc,*pland_a,**E_pl_aland,**E_pl_adisc;
  double  *E_lga_land,*E_wga_land,*E_lga_land_tot,*E_wga_land_tot,*E_lga_disc_tot,*E_wga_disc_tot;
  int     *nCell_land_tot,*nCell_disc_tot;
  double  prob,Sconst,Wconst,S,len,pa_sum_land,pa_sum_disc,pland,pland_l,pland_l_old,E_adisc,E_a_tot_disc;
  double  sum_l_aland,sum_l_adisc;
  char    string[150];

  E_lga = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_wga = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  mu_age = Mmatrix_2d(0,i_D_age->glm->nxcov-1,0,i_D_age->glm->ncat-1,sizeof(double),1);     // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  pa = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  E_pa = CALLOC(i_D_age->glm->ncat,double);        // Free ok
  E_pl_a = Mmatrix_2d(0,i_D_age->glm->ncat-1,0,s_D_totcatch->nlint-1,sizeof(double),1);        // Free ok

  E_lga_land = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_wga_land = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_lga_land_tot = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_wga_land_tot = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  nCell_land_tot = CALLOC(i_D_age->glm->ncat,int);      // Free ok
  E_lga_disc_tot = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_wga_disc_tot = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  nCell_disc_tot = CALLOC(i_D_age->glm->ncat,int);      // Free ok
  int_num_length = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  int_num_weight = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  E_pa_land = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  E_pa_disc = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  pland_a = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  E_pl_aland = Mmatrix_2d(0,i_D_age->glm->ncat-1,0,s_D_totcatch->nlint-1,sizeof(double),1);    // Free ok
  E_pl_adisc = Mmatrix_2d(0,i_D_age->glm->ncat-1,0,s_D_totcatch->nlint-1,sizeof(double),1);    // Free ok

  sum_T = G_ZERO;
  
  /* Initialize variables */
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      o_mean_l[a] = G_ZERO;
      o_mean_w[a] = G_ZERO;
      for(l=0;l<i_D_totcatch->nlint;l++)
	o_totcatch->catch_at_age[a][l] = G_ZERO;
    }
  if(s_COST)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  s_mean_l_disc[a] = G_ZERO;
	  s_mean_w_disc[a] = G_ZERO;
	  E_lga_land_tot[a] = G_ZERO;
	  E_wga_land_tot[a] = G_ZERO;
	  nCell_land_tot[a] = 0;
	  E_lga_disc_tot[a] = G_ZERO;
	  E_wga_disc_tot[a] = G_ZERO;
	  nCell_disc_tot[a] = 0;
	  for(l=0;l<i_D_totcatch->nlint;l++)
	    s_totcatch_disc->catch_at_age[a][l] = G_ZERO;
	}
    }

  Hsz = 0;
  if(i_D_age->glm->nxcov>1)
    Hsz = 1;
  if(i_D_lga->glm->nxcov>2)
    Hsz = 1;
  if(i_D_wgl->glm->nxcov>2)
    Hsz = 1;

  if(s_coastal_cod==0)
    {
      err = simulate_cell_effects(i_age,i_D_age,i_length,i_D_lga,i_weight,i_D_wgl);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling simulate_cell_effects\n");
	  return(err);
	}
    }
  else
    {
      err = simulate_cell_effects_CC(i_age,i_D_age,i_length,i_D_lga,i_weight,i_D_wgl,
				     s_length_CC,s_D_lga_CC,s_weight_CC,s_D_wgl_CC);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling simulate_cell_effects\n");
	  return(err);
	}
    }

  for(c=0;c<i_D_totcatch->nCell;c++)
    {
      season = i_season[c];
      if(Hsz)
	{
	  mu_hsz = i_par_haulsize[2*c];
	  sd_hsz = i_par_haulsize[2*c+1];
	}
      else
	{
	  mu_hsz = G_ZERO;
	  sd_hsz = G_ZERO;
	}

      #ifdef DEBUG_PREDICT
      printf("Calculate means for age\n");
      #endif
      err = calculate_means_age(mu_age,c,i_age,i_D_age,i_D_totcatch);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_age\n");
	  return(err);
	}

      #ifdef DEBUG_PREDICT
      printf("Calculate means for lga\n");
      #endif
      /* Calculate means for lga */
      err = calculate_means_lga(mu_lga,&var_lga,c,i_length,i_D_lga,i_inc_haul,i_D_totcatch);
      //printf("mu_lga[0]=%lf,mu_lga[1]=%lf,var_lga=%lf\n",mu_lga[0],mu_lga[1],var_lga);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_lga\n");
	  return(err);
	}
      if(s_coastal_cod)
	{
          #ifdef DEBUG_PREDICT
	  printf("Calculate means for lga - coastal cod\n");
          #endif
	  err = calculate_means_lga(mu_lga_CC,&var_lga_CC,c,s_length_CC,s_D_lga_CC,i_inc_haul,i_D_totcatch);
	  if(err)
	    {
	      write_warning("find_catch_at_age:Error calling calculate_means_lga\n");
	      return(err);
	    }
	}

      #ifdef DEBUG_PREDICT
      printf("Calculate means for wgl\n");
      #endif
      /* Calculate means for wgl */
      err = calculate_means_wgl(mu_wgl,&var_wgl,c,i_weight,i_D_wgl,i_inc_haul,i_D_totcatch);
      //printf("mu_wgl[0]=%lf,mu_wgl[1]=%lf,var_wgl=%lf\n",mu_wgl[0],mu_wgl[1],var_wgl);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_wgl\n");
	  return(err);
	}
      if(s_coastal_cod)
	{
          #ifdef DEBUG_PREDICT
	  printf("Calculate means for wgl- coastal cod\n");
          #endif
	  err = calculate_means_wgl(mu_wgl_CC,&var_wgl_CC,c,s_weight_CC,s_D_wgl_CC,i_inc_haul,i_D_totcatch);
	  if(err)
	    {
	      write_warning("find_catch_at_age:Error calling calculate_means_wgl\n");
	      return(err);
	    }
	}

      if(s_COST)// simulating discarded fish
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    mu[a] = mu_lga[0]+mu_lga[1]*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1];
	  x = s_D_COST->cens->mu[0] + gennor(G_ZERO,1/sqrt(s_D_COST->cens->tau[0]));
	  cens_k = exp(x)/(1+exp(x));
	  cens_m = s_D_COST->cens->mu[1] + gennor(G_ZERO,1/sqrt(s_D_COST->cens->tau[1]));
	}

      #ifdef DEBUG_PREDICT
      printf("Monte Carlo estimation\n");
      #endif
      /* Monte Carlo estimation  */
      sd = G_ONE/sqrt(i_age->par->tau_obs);
      for(a=0;a<i_D_age->glm->ncat;a++)
	E_pa[a] = G_ZERO;
      if(s_COST)
	{
	  pland = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      E_pa_land[a] = G_ZERO;
	      E_pa_disc[a] = G_ZERO;
	    }
	}
      mu_lga_hsz = G_ZERO;
      mu_wgl_hsz = G_ZERO;
      mu_lga_hsz_CC = G_ZERO;
      mu_wgl_hsz_CC = G_ZERO;
      for(h=0;h<i_nMC;h++)
	{
	  hsz = mu_hsz + sd_hsz*gennor(G_ZERO,G_ONE);
          pa_sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
              eff = mu_age[0][a] + gennor(G_ZERO,sd);
	      if(i_D_age->glm->nxcov>1)
		eff += mu_age[1][a]*hsz;
	      pa[a] = exp(eff);
	      pa_sum += pa[a];
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    E_pa[a] += pa[a]/pa_sum;
	  if(s_COST)
	    {
	      /* Calculate p(landed|a) in each age group */
	      /* Draw haul specific r parameter */
	      cens_r = s_D_COST->cens->mu[2] + gennor(G_ZERO,1/sqrt(s_D_COST->cens->tau[2]));
	      //fprintf(stderr,"k=%f,m=%f,r=%f\n",cens_k,cens_m,cens_r);
	      Sconst = cens_k/sqrt(G_PI);
	      sig = sqrt(var_lga);
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{  
		  int_num_length[a] = G_ZERO;
		  int_num_weight[a] = G_ZERO;
		  pland_a[a] = G_ZERO;
		  for(i=1;i<=s_N_gauher;i++)
		    {
		      len = s_gauher_x[i];
		      // using Gaussian cdf as S-function for probability of landing
		      prob = cens_m*(mu[a]+sqrt(2)*sig*len-cens_r);
		      S = pnorm(prob,G_ZERO,G_ONE);
		      int_num_length[a] += s_gauher_w[i] * exp(mu[a]+sqrt(2)*sig*len)*Sconst*S;
		      int_num_weight[a] += s_gauher_w[i] * exp(mu_wgl[1]*(mu[a]+sqrt(2)*sig*len))*Sconst*S;
		      pland_a[a] += s_gauher_w[i] * Sconst*S;
		    }
		  if(pland_a[a]>1.000001)
		    {
		      sprintf(string,"find_catch_at_age:Something is wrong, p(landed|a=%d)=%10.8lf\n",
			      a,pland_a[a]);
		      write_warning(string);
		      return(1);
		    }
		  if(pland_a[a]>0.99999999)
		    pland_a[a]=1.0;
		}
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  E_pa_land[a] += pa[a]/pa_sum*pland_a[a];
		  E_pa_disc[a] += pa[a]/pa_sum*(1-pland_a[a]);
		  pland += pa[a]/pa_sum*pland_a[a];
		}
	    }
	  if(i_D_lga->glm->nxcov>2)
	    mu_lga_hsz += mu_lga[2] * hsz;
	  if(i_D_wgl->glm->nxcov>2)
	    mu_wgl_hsz += mu_wgl[2] * hsz;
	  if(s_coastal_cod)
	    {
	      if(s_D_lga_CC->glm->nxcov>2)
		mu_lga_hsz_CC += mu_lga_CC[2] * hsz;
	      if(s_D_wgl_CC->glm->nxcov>2)
		mu_wgl_hsz_CC += mu_wgl_CC[2] * hsz;
	    }
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	E_pa[a] /= (double) i_nMC;

      if(s_COST)
	{
	  pa_sum_land = G_ZERO;
	  pa_sum_disc = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      pa_sum_land += E_pa_land[a];
	      pa_sum_disc += E_pa_disc[a];
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      E_pa_land[a] /= pa_sum_land;
	      E_pa_disc[a] /= pa_sum_disc;
	    }
	  pland /= (double) i_nMC;
	}
	

      mu_lga_hsz /= (double) i_nMC;
      mu_wgl_hsz /= (double) i_nMC;
      mu_lga[0] += mu_lga_hsz;
      mu_wgl[0] += mu_wgl_hsz;
      if(s_coastal_cod)
	{
	  mu_lga_hsz_CC /= (double) i_nMC;
	  mu_wgl_hsz_CC /= (double) i_nMC;
	  mu_lga_CC[0] += mu_lga_hsz_CC;
	  mu_wgl_CC[0] += mu_wgl_hsz_CC;
	}

      #ifdef DEBUG_PREDICT
      printf("mu_lga=(%lf,%lf),mu_wgl=(%lf,%lf)\n",
	     mu_lga[0],mu_lga[1],mu_wgl[0],mu_wgl[1]);
      xcov = i_D_lga->glm->xcov[0];
      printf("prec_lga_obs=%lf,prec_lga_haul=%lf,var_lga=%lf\n",
	     i_length->par->tau_obs,i_length->par->tau[0][xcov->n_cov-1],var_lga);
      xcov = i_D_wgl->glm->xcov[0];
      printf("prec_wgl_obs=%lf,prec_wgl_haul=%lf,var_wgl=%lf\n",
	     i_weight->par->tau_obs,i_weight->par->tau[0][xcov->n_cov-1],var_wgl);
      if(s_coastal_cod)
	{
	  printf("mu_lga_CC=(%lf,%lf),mu_wgl_CC=(%lf,%lf)\n",
		 mu_lga_CC[0],mu_lga_CC[1],mu_wgl_CC[0],mu_wgl_CC[1]);
	  xcov = s_D_lga_CC->glm->xcov[0];
	  printf("prec_lga_obs=%lf,prec_lga_haul=%lf,var_lga=%lf\n",
		 s_length_CC->par->tau_obs,s_length_CC->par->tau[0][xcov->n_cov-1],var_lga_CC);
	  xcov = s_D_wgl_CC->glm->xcov[0];
	  printf("prec_wgl_obs=%lf,prec_wgl_haul=%lf,var_wgl=%lf\n",
		 s_weight_CC->par->tau_obs,s_weight_CC->par->tau[0][xcov->n_cov-1],var_wgl_CC);
	}
      if(s_COST)
	printf("mu_age E[p(a)] E[p(a|landed)] E_lga E_wga\n");
      else
	printf("mu_age E[p(a)] E_lga E_wga\n");
      #endif

      if(s_COST)// simulating discarded fish
	{
	  mean_w = G_ZERO;
	  mean_l = G_ZERO;  
	  A = mu_wgl[0] + mu_wgl[1]*mu_lga[0];
	  B = mu_wgl[1]*mu_lga[1];
	  var = mu_wgl[1]*mu_wgl[1]*var_lga + var_wgl;
	  sig = sqrt(var_lga);
	  Wconst = exp(mu_wgl[0]+G_HALF*var_wgl);

	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {  
	      E_lga[a] = exp(mu_lga[0]+mu_lga[1]*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var_lga);
	      E_wga[a] = exp(A+B*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var);
	      if(pland_a[a] > 0.0)
		{
		  mean_l += E_pa_land[a] * int_num_length[a] / pland_a[a];
		  mean_w += E_pa_land[a] * Wconst * int_num_weight[a] / pland_a[a];
		  E_lga_land[a] = int_num_length[a] / pland_a[a];
		  E_wga_land[a] = Wconst * int_num_weight[a] / pland_a[a];
		  E_lga_land_tot[a] += E_lga_land[a];
		  E_wga_land_tot[a] += E_wga_land[a];
		  nCell_land_tot[a] += 1;
		}
	      else
		{
		  mean_l += G_ZERO;
		  mean_w += G_ZERO;
		  E_lga_land[a] = G_ZERO;
		  E_wga_land[a] = G_ZERO;
		}
              #ifdef DEBUG_PREDICT
	      printf("%lf %lf %lf %lf %lf\n",mu_age[0][a],E_pa[a],E_pa_land[a],E_lga_land[a],E_wga_land[a]);
              #endif

	      if(s_D_totcatch->nlint>1)
		{
		  E_pl_a[a][0] = pnorm(s_D_totcatch->l_int[1],mu[a],sig);
		  pland_l = cens_k*pnorm(cens_m*(s_D_totcatch->l_int[1]-cens_r),G_ZERO,G_ONE);
		  E_pl_aland[a][0] = E_pl_a[a][0]*pland_l;
		  sum_l_aland = E_pl_aland[a][0];
		  E_pl_adisc[a][0] = E_pl_a[a][0]*(1-pland_l);
		  sum_l_adisc = E_pl_adisc[a][0];
		  sum = E_pl_a[a][0];
		  pland_l_old = pland_l;
		  for(l=1;l<(s_D_totcatch->nlint-1);l++)
		    {
		      E_pl_a[a][l] = pnorm(s_D_totcatch->l_int[l+1],mu[a],sig);
		      E_pl_a[a][l] -= sum;
		      sum += E_pl_a[a][l];
		      pland_l = cens_k*pnorm(cens_m*(s_D_totcatch->l_int[l+1]-cens_r),G_ZERO,G_ONE);
		      E_pl_aland[a][l] = E_pl_a[a][l]*pland_l;
		      sum_l_aland += E_pl_aland[a][l];
		      E_pl_adisc[a][l] = E_pl_a[a][l]*(1-pland_l);
		      sum_l_adisc += E_pl_adisc[a][l];
		      pland_l_old = pland_l;
		    }
		  E_pl_a[a][s_D_totcatch->nlint-1] = G_ONE-sum;
		  E_pl_aland[a][s_D_totcatch->nlint-1] = G_ONE-sum;
		  E_pl_adisc[a][s_D_totcatch->nlint-1] = G_ZERO;
		  if(sum_l_aland == 0)
		    E_pl_aland[a][0] = G_ONE;
		  else
		    {
		      for(l=0;l<s_D_totcatch->nlint;l++)
			E_pl_aland[a][l] /= sum_l_aland;
		    }
		  if(sum_l_adisc == 0)
		    E_pl_adisc[a][0] = G_ONE;
		  else
		    {
		      for(l=0;l<s_D_totcatch->nlint;l++)
			E_pl_adisc[a][l] /= sum_l_adisc;
		    }
		  //if(c==0)
		  if(0)
		    {
		      fprintf(stderr,"a=%d:mu=%f,E_lga=%f,E_wga=%f,mean_l=%f,mean_w=%f\n",
			      a,mu[a],E_lga_land[a],E_wga_land[a],mean_l,mean_w);
		      for(l=0;l<s_D_totcatch->nlint;l++)
			fprintf(stderr,"%5.4f ",E_pl_a[a][l]);
		      fprintf(stderr,"\n");
		      for(l=0;l<s_D_totcatch->nlint;l++)
			fprintf(stderr,"%5.4f ",E_pl_aland[a][l]);
		      fprintf(stderr,"\n");
		      for(l=0;l<s_D_totcatch->nlint;l++)
			fprintf(stderr,"%5.4f ",E_pl_adisc[a][l]);
		      fprintf(stderr,"\n");
		      //for(l=0;l<s_D_totcatch->nlint;l++)
		      //{
		      //  fprintf(stderr,"p(l|a=%d)=%f,p(l|a,land)=%f,p(l|a,disc)=%f\n",
		      //	  a,E_pl_a[a][l],E_pl_aland[a][l],E_pl_adisc[a][l]);
		      //}
		    }
		}
	      else
		E_pl_a[a][0] = G_ONE;
	    }//end for(a=0;a<i_D_age->glm->ncat;a++)

	}//end if(s_length->cens_model)
      else // if not simulating discarded fish 
	{
	  mean_w = G_ZERO;
	  mean_l = G_ZERO;   
	  if(s_coastal_cod)
	    {
	      A = mu_wgl_CC[0] + mu_wgl_CC[1]*mu_lga_CC[0];
	      B = mu_wgl_CC[1]*mu_lga_CC[1];
	      var = mu_wgl_CC[1]*mu_wgl_CC[1]*var_lga_CC + var_wgl_CC;
              #ifdef DEBUG_PREDICT
	      printf("CC:A=%f,B=%f,var=%f\n");
	      #endif
	      a_start = 0;
	      a_end = (int) i_D_age->glm->ncat/2;
	      for(a=a_start;a<a_end;a++)
		{
		  E_lga[a] = exp(mu_lga_CC[0]+mu_lga_CC[1]*s_D_g_a_CC->g_a[s_D_g_a_CC->a2Age_vec[a]+season-1]+G_HALF*var_lga_CC);
		  E_wga[a] = exp(A+B*s_D_g_a_CC->g_a[s_D_g_a_CC->a2Age_vec[a]+season-1]+G_HALF*var);
		  mean_l += E_pa[a] * E_lga[a];
		  mean_w += E_pa[a] * E_wga[a];
                  #ifdef DEBUG_PREDICT
		  printf("%lf %lf %lf %lf\n",mu_age[0][a],E_pa[a],E_lga,E_wga);
                  #endif
		  mu[a] = mu_lga_CC[0]+mu_lga_CC[1]*s_D_g_a_CC->g_a[s_D_g_a_CC->a2Age_vec[a]+season-1];
		  sig = sqrt(var_lga_CC);
		  if(s_D_totcatch->nlint>1)
		    {
		      E_pl_a[a][0] = pnorm(s_D_totcatch->l_int[1],mu[a],sig);
		      sum = E_pl_a[a][0];
		      for(l=1;l<(s_D_totcatch->nlint-1);l++)
			{
			  E_pl_a[a][l] = pnorm(s_D_totcatch->l_int[l+1],mu[a],sig);
			  E_pl_a[a][l] -= sum;
			  sum += E_pl_a[a][l];
			}
		      E_pl_a[a][s_D_totcatch->nlint-1] = G_ONE-sum;
		    }
		  else
		    E_pl_a[a][0] = G_ONE;
		}
	    }//end if(s_coastal_cod)	 
	  A = mu_wgl[0] + mu_wgl[1]*mu_lga[0];
	  B = mu_wgl[1]*mu_lga[1];
	  var = mu_wgl[1]*mu_wgl[1]*var_lga + var_wgl;
	  if(s_coastal_cod)
	    a_start = i_D_age->glm->ncat/2;
	  else
	    a_start = 0;
	  for(a=a_start;a<i_D_age->glm->ncat;a++)
	    {
	      E_lga[a] = exp(mu_lga[0]+mu_lga[1]*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var_lga);
	      E_wga[a] = exp(A+B*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var);
	      mean_l += E_pa[a] * E_lga[a];
	      mean_w += E_pa[a] * E_wga[a];
              #ifdef DEBUG_PREDICT
	      printf("%lf %lf %lf %lf\n",mu_age[0][a],E_p[a],E_lga,E_wga);
              #endif
	      mu[a] = mu_lga[0]+mu_lga[1]*s_D_g_a->g_a[s_D_g_a->a2Age_vec[a]+season-1];
	      sig = sqrt(var_lga);
	      if(s_D_totcatch->nlint>1)
		{
		  E_pl_a[a][0] = pnorm(s_D_totcatch->l_int[1],mu[a],sig);
		  sum = E_pl_a[a][0];
		  for(l=1;l<(s_D_totcatch->nlint-1);l++)
		    {
		      E_pl_a[a][l] = pnorm(s_D_totcatch->l_int[l+1],mu[a],sig);
		      E_pl_a[a][l] -= sum;
		      sum += E_pl_a[a][l];
		    }
		  E_pl_a[a][s_D_totcatch->nlint-1] = G_ONE-sum;
		}
	      else
		E_pl_a[a][0] = G_ONE;
	    }
	}//end else (if not simulating discarded fish)

      T = i_D_totcatch->catch[c]/mean_w;
      if(!(T> -9999.0  && T < 999999999999999999.99))
	{
	  write_warning("find_catch_at_age:Something is wrong 2\n");
	  fprintf(stderr,"mu_age E_p\n");
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      for(i=0;i<i_D_age->glm->nxcov;i++)
		fprintf(stderr,"%lf ",mu_age[i][a]);
	      fprintf(stderr,"%lf\n",E_pa[a]);
	    }
	  fprintf(stderr,"mean_w=%lf,sd=%lf\n",mean_w,sd);
	  fprintf(stderr,"mu_lga = ");
	  for(i=0;i<i_D_lga->glm->nxcov;i++)
	    fprintf(stderr,"%lf ",mu_lga[i]);
	  fprintf(stderr,"mu_wgl = ");
	  for(i=0;i<i_D_wgl->glm->nxcov;i++)
	    fprintf(stderr,"%lf ",mu_wgl[i]);
	  fprintf(stderr,"A=%lf, B=%lf,var=%lf\n",A,B,var);
	  return(1);
	}
      #ifdef DEBUG_PREDICT
      printf("c=%d,catch=%lf,A=%lf,B=%lf,var=%lf,mean_l=%lf,mean_w=%lf,T=%lf\n",
	     c,i_D_totcatch->catch[c],A,B,var,mean_l,mean_w,T);
      #endif
 
      sum_T += T;
      double d1,d2;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  E_a = 0.0;
	  E_adisc = 0.0;
	  sum = G_ZERO;
	  for(l=0;l<i_D_totcatch->nlint;l++)
	    {
	      if(s_COST)
		{
		  d1 = T * E_pa_land[a]*E_pl_aland[a][l];
		  d2 = T*(1-pland)/pland*E_pa_disc[a]*E_pl_adisc[a][l];
		  o_totcatch->catch_at_age[a][l] += T * E_pa_land[a]*E_pl_aland[a][l];
		  s_totcatch_disc->catch_at_age[a][l] += T*(1-pland)/pland*E_pa_disc[a]*E_pl_adisc[a][l];
		  E_a += T * E_pa_land[a]*E_pl_aland[a][l];
		  E_adisc += T*(1-pland)/pland*E_pa_disc[a]*E_pl_adisc[a][l];
		  //if(d1>1E3 || d2 > 1E3)
		  if(d2<0)
		    {
		      printf("iter=%d,c=%d,a=%d,l=%d\n",iter,c,a,l);
		      printf("catch=%f,E_a=%f,T=%f,pland=%f\n",
			     d1,E_a,T,pland);
		      printf("  E_pa_land=%f,E_pl_aland[a][l]=%f\n",
			     E_pa_land[a],E_pl_aland[a][l]);
		      printf("disc=%12.10f,E_adisc=%12.10f,T=%12.10f,pdisc=%12.10f\n",
			     d2,E_adisc,T,(1-pland)/pland);
		      printf(" E_pa_disc=%12.10f,E_pl_adisc[a][l]=%12.10f\n",
			     E_pa_disc[a],E_pl_adisc[a][l]);
		    }
		}
	      else
		{
		  o_totcatch->catch_at_age[a][l] += T * E_pa[a]*E_pl_a[a][l];
		  E_a += T * E_pa[a]*E_pl_a[a][l];
		}
	      sum += o_totcatch->catch_at_age[a][l];
              #ifdef DEBUG_PREDICT
	      printf("T=%lf,E_pa=%lf,E_pl_a=%lf,catch_at_age=%lf\n",
		     T,E_pa[a],E_pl_a[a][l],o_totcatch->catch_at_age[a][l]);
              #endif
	    }
	  if(s_COST)
	    {
	      o_mean_l[a] += E_lga_land[a]*E_a;
	      o_mean_w[a] += E_wga_land[a]*E_a;
	      if(pland_a[a]<1.0)
		{
		  s_mean_l_disc[a] += (E_lga[a]-E_lga_land[a]*pland_a[a])/(1-pland_a[a])*E_adisc;
		  s_mean_w_disc[a] += (E_wga[a]-E_wga_land[a]*pland_a[a])/(1-pland_a[a])*E_adisc;
		  E_lga_disc_tot[a] += (E_lga[a]-E_lga_land[a]*pland_a[a])/(1-pland_a[a]);
		  E_wga_disc_tot[a] += (E_wga[a]-E_wga_land[a]*pland_a[a])/(1-pland_a[a]);
		  nCell_disc_tot[a] += 1;
		}
	    }
	  else
	    {
	      o_mean_l[a] += E_lga[a]*E_a;
	      o_mean_w[a] += E_wga[a]*E_a;
	    }
	}
      if(s_COST)
	s_planded[c] = pland;
    }// end for(c=0;c<i_D_totcatch->nCell;c++)

  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      E_a_tot = 0;
      for(l=0;l<i_D_totcatch->nlint;l++)
	E_a_tot += o_totcatch->catch_at_age[a][l];

      if(E_a_tot<1E-10)
	{
	  if(s_COST)
	    {
	      E_lga_land_tot[a] /= nCell_land_tot[a];
	      E_wga_land_tot[a] /= nCell_land_tot[a];
	      if(E_lga_land_tot[a]>0)
		o_mean_l[a] = E_lga_land[a];
	      else
		o_mean_l[a] /= G_ZERO;
	      if(E_wga_land_tot[a]>0)
		o_mean_w[a] = E_wga_land[a];
	      else
		o_mean_w[a] /= G_ZERO;
	    }
	  else
	    {
	      if(E_lga[a]>0)
		o_mean_l[a] = E_lga[a];
	      else
		o_mean_l[a] /= G_ZERO;
	      if(E_wga[a]>0)
		o_mean_w[a] = E_wga[a];
	      else
		o_mean_w[a] /= G_ZERO;
	    }
	}
      else
	{
	  o_mean_l[a] /= E_a_tot;
	  o_mean_w[a] /= E_a_tot;
	}

      if(s_COST)
	{
	  E_a_tot_disc = 0;
	  for(l=0;l<i_D_totcatch->nlint;l++)
	    E_a_tot_disc += s_totcatch_disc->catch_at_age[a][l];

	  if(E_a_tot_disc<1E-10)
	    {
	      E_lga_disc_tot[a] /= nCell_disc_tot[a];
	      E_wga_disc_tot[a] /= nCell_disc_tot[a];
	      if(E_lga_disc_tot[a]>0)
		s_mean_l_disc[a] = E_lga_disc_tot[a];
	      else
		s_mean_l_disc[a] /= G_ZERO;
	      if(E_wga_disc_tot[a]>0)
		s_mean_w_disc[a] = E_wga_disc_tot[a];
	      else
		s_mean_w_disc[a] /= G_ZERO;
	    }
	  else
	    {
	      s_mean_l_disc[a] /= E_a_tot_disc;
	      s_mean_w_disc[a] /= E_a_tot_disc;
	    }
	}
    }


  #ifdef DEBUG_PREDICT
  sum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      printf("a=%d,p_a=%lf\n",a,E_pa[a]);
      for(l=0;l<i_D_totcatch->nlint;l++)
	{
	  sum += o_totcatch->catch_at_age[a][l];
	  printf("%lf ",o_totcatch->catch_at_age[a][l]);
	}
      printf("\n");
    }
  #endif

  FREE(E_lga);
  FREE(E_wga);
  Fmatrix_2d(&mu_age[0][0],&mu_age[0]);
  FREE(mu);
  FREE(pa);
  FREE(E_pa);
  Fmatrix_2d(&E_pl_a[0][0],&E_pl_a[0]);

  FREE(E_lga_land);
  FREE(E_wga_land);
  FREE(E_lga_land_tot);
  FREE(E_wga_land_tot);
  FREE(nCell_land_tot);
  FREE(E_lga_disc_tot);
  FREE(E_wga_disc_tot);
  FREE(nCell_disc_tot);
  FREE(int_num_length);
  FREE(int_num_weight);
  FREE(E_pa_land);
  FREE(E_pa_disc);
  FREE(pland_a);
  Fmatrix_2d(&E_pl_aland[0][0],&E_pl_aland[0]);
  Fmatrix_2d(&E_pl_adisc[0][0],&E_pl_adisc[0]);

  return(0);
}

static int calculate_means_age(double **mu_age,int c,Age_struct *i_age,Data_age *i_D_age,
			       Data_totcatch *i_D_totcatch)
{
  int a,i,j,j2,k,ncov;
  double eff;
  Data_cov *xcov;
  
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      xcov = i_D_age->glm->xcov[i];
      ncov = xcov->n_cov;
      if(i==0)
	ncov--;  // Excluding haul effects in intercept since this is treated later
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  mu_age[i][a] = G_ZERO;
	  for(j=0;j<ncov;j++)
	    {
	      k = i_D_totcatch->age_xcov[i]->c_cov[c][j];
	      if(j!=xcov->icell && k > -1)
		eff = i_age->par->eff[a][i][j][k];
	      else if(j==xcov->icell)
		eff = i_age->par->cell[i][k];
	      else if(j==xcov->ispat)
		{
                  #ifdef LOG_FILE
		  fprintf(g_caa_log,
			  "calculate_means_age:Missing area c=%d not implemented yet\n",c);
		  fprintf(g_caa_log,"Age model, a=%d,j=%d,k=%d\n",a,j,k);
		  fprintf(g_caa_log,"c=%d",c);
		  for(k=0;k<xcov->n_cov;k++)
		    fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[i]->c_cov[c][k]);
		  fprintf(g_caa_log,"\n");
                  #endif
		  write_warning("calculate_means_age:Something is wrong\n");
		  return(1);
		}
	      else if(i_D_age->glm->xcov[i]->fix[j]==0&&j != xcov->ispat)
		eff = gennor(G_ZERO,G_ONE/sqrt(i_age->par->tau[i][j]));   //cell effect only option
	      else
		{
                  #ifdef LOG_FILE
		  for(j2=0;j2<xcov->n_cov;j2++)
		    fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[i]->c_cov[c][j2]);
		  fprintf(g_caa_log,"\n");
		  fprintf(g_caa_log,"Missing fixed effect in age model (c=%d,j=%d) not allowed\n",c,j);
                  #endif
		  write_warning("calculate_means_age:Something is wrong\n");
		      return(1);
		}
	      mu_age[i][a] += eff;
	    }
	}
    } 
  return(0);
}            /* End of calculate_means_age */

static int calculate_means_lga(double *mu_lga,double *o_var_lga,int c,LW_struct *i_length,Data_lin *i_D_lga,
			       int *i_inc_haul,Data_totcatch *i_D_totcatch)
{
  int i,j,k,ncov;
  double var_lga,eff;
  Data_cov *xcov;
  
  var_lga = G_ONE/i_length->par->tau_obs;
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      mu_lga[i] = G_ZERO;
      ncov = xcov->n_cov;
      if(i==0 && i_inc_haul[2]==1)
	{
	  // Include haul effect in variance and not in effect
	  var_lga += G_ONE/i_length->par->tau[i][xcov->n_cov-1];
	  ncov -= 1;   
	}
      for(j=0;j<ncov;j++)
	{
	  k=i_D_totcatch->lga_xcov[i]->c_cov[c][j];
	  if(j!= xcov->icell && k > -1)
	    eff = i_length->par->eff[0][i][j][k];
	  else if (j==xcov->icell)
	    eff = i_length->par->cell[i][k];
	  else if(i_D_lga->glm->xcov[i]->fix[j]==0&&j != xcov->ispat) //cell and haul effect
	    eff = gennor(G_ZERO,sqrt(G_ONE/i_length->par->tau[i][j]));
	  else if(j==xcov->ispat)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"calculate_means_lga:Missing area c=%d not implemented yet\n",c);
	      fprintf(g_caa_log,"c=%d",c);
	      for(k=0;k<xcov->n_cov;k++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[i]->c_cov[c][k]);
	      fprintf(g_caa_log,"\n");
              #endif
	      write_warning("calculate_means_lga:Missing area effect in age model is not allowed\n");
	      return(1);
	    }
	  else
	    {
	      var_lga += G_ONE/i_length->par->tau[i][j];
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"Missing fixed effect in lga model (c=%d,j=%d) not allowed\n",c,j);
              #endif
	      write_warning("calculate_means_lga:Something is wrong\n");
	      return(1);
	    }
	  mu_lga[i] += eff;
	}
    }
  o_var_lga[0] = var_lga;

  return(0);
}            /* End of calculate_means_lga */

static int calculate_means_wgl(double *mu_wgl,double *o_var_wgl,int c,LW_struct *i_weight,Data_lin *i_D_wgl,
			       int *i_inc_haul,Data_totcatch *i_D_totcatch)
{
  int i,j,k,ncov;
  double var_wgl,eff;
  Data_cov *xcov;

  var_wgl = G_ONE/i_weight->par->tau_obs;
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      ncov = xcov->n_cov;
      mu_wgl[i] = G_ZERO;
      if(i==0 && i_inc_haul[5]==1)
	{
	  // Include haul effect in variance and not in effect
	  var_wgl += G_ONE/i_weight->par->tau[i][xcov->n_cov-1];
	  ncov -= 1;   
	}
      for(j=0;j<ncov;j++)
	{
	  k=i_D_totcatch->wgl_xcov[i]->c_cov[c][j];
	  if(j != xcov->icell && k > -1)
	    eff = i_weight->par->eff[0][i][j][k];
	  else if (j==xcov->icell)
	    eff = i_weight->par->cell[i][k];
	  else if(i_D_wgl->glm->xcov[i]->fix[j]==0&&j != xcov->ispat) //cell and haul effect
	    eff = gennor(G_ZERO,sqrt(G_ONE/i_weight->par->tau[i][j]));
	  else if(j==xcov->ispat)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"calculate_means_wgl:Missing area c=%d not implemented yet\n",c);
	      fprintf(g_caa_log,"c=%d",c);
	      for(k=0;k<xcov->n_cov;k++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[i]->c_cov[c][k]);
	      fprintf(g_caa_log,"\n");
              #endif
	      write_warning("calculate_means_wgl:Something is wrong\n");
	      return(1);
	    }
	  else
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"Missing fixed effect in wgl model (c=%d,j=%d) not allowed\n",c,j);
              #endif
	      write_warning("calculate_means_wgl:Something is wrong\n");
	      return(1);
	    }
	  mu_wgl[i] += eff;
	}
    }
  o_var_wgl[0] = var_wgl;
  return(0);
}            /* End of calculate_means_wgl */
