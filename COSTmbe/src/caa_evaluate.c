/*!
  \file caa_evaluate.c
  \brief Routines for model evaluations
  \author Geir Storvik

*/


#include "caa.h"


static int                 s_nMCMC;
static int                 s_burnin;
static Age_struct         *s_age;
static LW_struct          *s_length;
static LW_struct          *s_weight;
static Data_orig          *s_D_orig;
static Data_age           *s_D_age;
static Data_lin           *s_D_lga;
static Data_lin           *s_D_wgl;
static Data_g_a           *s_D_g_a;
static int                 s_N2=10;       
static double             *s_eps;
static double             *s_d_rep;
static double             *s_g;
static double             *s_entr;

static int initialize(int i_g_a_model);
static int re_initialize(int i_g_a_model);


/*!
  \author Geir Storvik
  \brief Initialization for evaluation routines

  At present only used for opening file for debugging.
*/
int init_evaluate(int i_nHaul,int i_ncat)
{
  int    b,i,j,k1,k2,n;
  double d_i,d_N;
  static int (*cmp)();

  #ifdef DEBUG_EVALUATE_FILE
  s_unit = fopen("caa_evaluate.txt","w");
  #endif

  #ifdef DEBUG_IND_KS
  s_ind_ks = fopen("caa_ind_ks.txt","w");
  #endif


  s_eps = CALLOC(i_nHaul*i_ncat+1,double);

  cmp = compd;
  n = i_nHaul*i_ncat;
  k1 = (int) ((double) n * 0.1);
  k2 = n-k1+1;
  s_d_rep = CALLOC(s_N2,double);
  s_g = CALLOC(s_N2,double);

  for(i=0;i<s_N2;i++)
    {
      for(j=0;j<n;j++)
	s_eps[j] = gennor(G_ZERO,G_ONE);
      qsort(s_eps,n,sizeof(double),cmp);
      s_d_rep[i] = fabs(fabs(s_eps[k2]) - fabs(s_eps[k1]));
    }
  qsort(s_d_rep,s_N2,sizeof(double),cmp);
  d_i = G_ONE;
  d_N = (double) s_N2;
  b = -1;
  while(s_d_rep[b+1]<(d_i*s_d_rep[s_N2-1]/d_N))
    b++;
  s_g[0] = b;
  for(i=0;i<s_N2;i++)
    {
      d_i = (double) (i+1);
      while(s_d_rep[b+1]<(d_i*s_d_rep[s_N2-1]/d_N)&&b<s_N2)
	b++;
      s_g[i] = b;
    }

  s_entr = CALLOC(i_nHaul,double);

  return(0);
}		/* end of init_evaluate */



/*!
  \author Geir Storvik
  \brief Re-initialize thing initialized in initialize_evaluate
*/
int re_init_evaluate(int i_ncat)
{
  #ifdef DEBUG_EVALUATE_FILE
  fclose(s_unit);
  #endif

  #ifdef DEBUG_IND_KS
  fclose(s_ind_ks);
  #endif

  FREE(s_eps);
  FREE(s_g);
  FREE(s_d_rep);
  FREE(s_entr);

   return(0);
}		/* end of re_init_evaluate */


/*!
  \author Geir Storvik
  \brief Calculates Bayes factors using as input simulations
  from a fit.

  Since this is a routine made to be called from R or Splus, all input 
  variables are pointers to vectors. The routine therefore starts to convert 
  input data into approperiate c-structures (as defined in caa.h).

*/
void caa_marg_dens(int *i_Nmcmc,int *i_burnin,double *i_mcmc1,double *i_mcmc2,
		   int *i_num_par1,int *i_num_par2,
		   int *age_nBoats,
		   int *i_totage,double *i_totlength,
		   double *i_totweight,double *i_haulweight,int *i_season,
		   int *i_nFishBoat,int *i_replength,
		   int *i_start_noAge,int *i_num_noAge,
		   int *nAges,int *a_vec,int *i_n_cov,int *i_ispat,
		   int *i_age_int_nFac,int *i_age_int_fix,int *i_age_int_c_cov,
		   int *i_age_hsz_nFac,int *i_age_hsz_fix,int *i_age_hsz_c_cov,
		   int *i_lga_nBoats,
		   int *i_lga_int_nFac,int *i_lga_int_fix,int *i_lga_int_c_cov,
		   int *i_lga_slp_nFac,int *i_lga_slp_fix,int *i_lga_slp_c_cov,
		   int *i_lga_hsz_nFac,int *i_lga_hsz_fix,int *i_lga_hsz_c_cov,
		   int *i_lga_g_a_model,int *i_lga_g_a_ncat,int *i_lga_g_a_nSeason,
		   double *i_lga_g_a_avec,int *i_lga_g_a_a2Age_vec,
		   int *i_wgl_nBoats,
		   int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
		   int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
		   int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
		   int *i_num_adj_area,int *i_adj_area,
		   double *o_loglik,double *o_logprior,double *o_logposterior,
		   int *o_err)
{
  int    i,it;

  #ifdef LOG_FILE
   g_caa_log = fopen("caa_logfile_marg.txt","w");
  #endif

  #ifdef DEBUG_INPUT
     *o_err = write_input_marg_dens(i_Nmcmc,i_burnin,i_mcmc1,i_mcmc2,
                i_num_par1,i_num_par2,age_nBoats,nAges,a_vec,
		i_n_cov,i_ispat,
   	        i_age_int_nFac,i_age_int_fix,i_age_int_c_cov,
   	        i_age_hsz_nFac,i_age_hsz_fix,i_age_hsz_c_cov,
   		i_lga_nBoats,
   		i_lga_int_nFac,i_lga_int_fix,i_lga_int_c_cov,
   		i_lga_slp_nFac,i_lga_slp_fix,i_lga_slp_c_cov,
   		i_lga_hsz_nFac,i_lga_hsz_fix,i_lga_hsz_c_cov,
		i_lga_g_a_model,
   		i_wgl_nBoats,
                i_wgl_int_nFac,i_wgl_int_fix,i_wgl_int_c_cov,
   		i_wgl_slp_nFac,i_wgl_slp_fix,i_wgl_slp_c_cov,
   		i_wgl_hsz_nFac,i_wgl_hsz_fix,i_wgl_hsz_c_cov,
   		i_num_adj_area,i_adj_area);
  #endif

  s_nMCMC = *i_Nmcmc;
  s_burnin = *i_burnin;

  #ifdef LOG_FILE
   fprintf(g_caa_log,"Initializing predict\n");
  #endif


  /* Make original data */
  s_D_orig = CALLOC(1,Data_orig);     // Free ok
  s_D_orig->totage = i_totage;
  s_D_orig->totlength = i_totlength;
  s_D_orig->haulweight = i_haulweight;
  s_D_orig->season = i_season;
  s_D_orig->nFishBoat = i_nFishBoat;
  s_D_orig->replength = i_replength;
  s_D_orig->start_noAge = i_start_noAge;
  s_D_orig->num_noAge = i_num_noAge;
  //  s_D_orig->n_int_len
  //  s_D_orig->int_len_lim

  if(*o_err)
    {
      write_warning("caa_marg_dens:tError calling makedata_orig\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }  

  /* Make age data */
  *o_err = makedata_age1(*age_nBoats,*nAges,a_vec,
			 i_n_cov[0],i_age_int_nFac,i_ispat[0],
			 i_age_int_fix,i_age_int_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[1],i_age_hsz_nFac,i_ispat[1],
			 i_age_hsz_fix,i_age_hsz_c_cov,
			 i_num_adj_area,i_adj_area,
			 &s_D_age);
  if(*o_err)
    {
      write_warning("caa_predic:tError calling makedata_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  /* Make g_a data */
  *o_err = makedata_g_a(*i_lga_g_a_ncat,*i_lga_g_a_nSeason,i_lga_g_a_avec,i_lga_g_a_a2Age_vec,
			*i_lga_g_a_model,s_D_age,&s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_marg_dens:Error calling makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  /* Make lga data */
  *o_err = makedata_lin1(*i_lga_nBoats,
			 i_n_cov[2],i_lga_int_nFac,i_ispat[2],
			 i_lga_int_fix,i_lga_int_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[3],i_lga_slp_nFac,i_ispat[3],
			 i_lga_slp_fix,i_lga_slp_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[4],i_lga_hsz_nFac,i_ispat[4],
			 i_lga_hsz_fix,i_lga_hsz_c_cov,
			 i_num_adj_area,i_adj_area,
			 &s_D_lga);

  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  s_D_lga->haulweight = i_haulweight;

  *o_err = makedata_lga_suff(s_D_lga,s_D_orig,s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling makedata_lga_suff\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  /* Read weight given length data */
  *o_err = makedata_lin1(*i_wgl_nBoats,
			 i_n_cov[5],i_wgl_int_nFac,i_ispat[5],
			 i_wgl_int_fix,i_wgl_int_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[6],i_wgl_slp_nFac,i_ispat[6],
			 i_wgl_slp_fix,i_wgl_slp_c_cov,
			 i_num_adj_area,i_adj_area,
			 i_n_cov[7],i_wgl_hsz_nFac,i_ispat[7],
			 i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			 i_num_adj_area,i_adj_area,
			 &s_D_wgl);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  s_D_wgl->haulweight = i_haulweight;

  *o_err = makedata_wgl_suff(s_D_wgl,i_nFishBoat,i_totlength,i_totweight,i_replength,i_haulweight);
  if(*o_err)
    {
      write_warning("caa_main_model2:Error calling makedata_lin2\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  *o_err = initialize(*i_lga_g_a_model);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
   
  s_age->par->num_var = i_num_par1[0];
  s_age->par->mcmc = &i_mcmc1[0];
  s_length->par->num_var = i_num_par1[1];
  s_length->par->mcmc = &i_mcmc1[i_num_par1[0]*s_nMCMC];
  s_D_g_a->g_a_mcmc = &i_mcmc1[(i_num_par1[0]+i_num_par1[1])*s_nMCMC];
  s_weight->par->num_var = i_num_par2[0];
  s_weight->par->mcmc = &i_mcmc2[0];

  #ifdef LOG_FILE
  fprintf(g_caa_log,"Perform prediction\n");
  #endif

  //Bayes factor for weight-given-length model
  //First pick out last simulation
  it = s_nMCMC-1;
  s_D_age->glm->xcov[0]->n_cov--;
  *o_err= read_it(it,s_D_age->glm,s_age->par);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling read_it\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  s_D_age->glm->xcov[0]->n_cov++;
  *o_err = read_it(it,s_D_lga->glm,s_length->par);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling read_it\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  *o_err = read_it(it,s_D_wgl->glm,s_weight->par);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling read_it\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  *o_err = read_it_g_a(it,s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling read_it_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  for(i=0;i<2;i++)
    {
      o_loglik[i] = G_ZERO;
      o_logprior[i] = G_ZERO;
      o_logposterior[i] = G_ZERO;
    }

  //wgl model

  //Calculate density for data
  *o_err = calc_lik_lin(s_weight,s_D_wgl,&(s_weight->par->loglik));
  o_loglik[1] += s_weight->par->loglik;
 
  //Calculate prior density for parameters

  //Calculate density for fixed and random effects given precisions

  //Calculate density for precisions

  #ifdef LOG_FILE
  fprintf(g_caa_log,"Cleaning up\n");
  #endif

  // Clean up
  *o_err = re_initialize(*i_lga_g_a_model);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling re_initialize\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  *o_err = re_makedata_lga_suff(s_D_lga);
  if(*o_err)
    {
      write_warning("caa_main_model1:Error calling re_makedata_lin2\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  *o_err = re_makedata_lin1(*i_wgl_nBoats,
			    i_n_cov[5],i_wgl_int_nFac,i_ispat[5],
			    i_wgl_int_fix,i_wgl_int_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[6],i_wgl_slp_nFac,i_ispat[6],
			    i_wgl_slp_fix,i_wgl_slp_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[7],i_wgl_hsz_nFac,i_ispat[7],
			    i_wgl_hsz_fix,i_wgl_hsz_c_cov,
			    i_num_adj_area,i_adj_area,
			    &s_D_wgl);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }
  *o_err = re_makedata_lin1(*i_lga_nBoats,
			    i_n_cov[2],i_lga_int_nFac,i_ispat[2],
			    i_lga_int_fix,i_lga_int_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[3],i_lga_slp_nFac,i_ispat[3],
			    i_lga_slp_fix,i_lga_slp_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[4],i_lga_hsz_nFac,i_ispat[4],
			    i_lga_hsz_fix,i_lga_hsz_c_cov,
			    i_num_adj_area,i_adj_area,
			    &s_D_lga);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling re_makedata_lin1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  *o_err = re_makedata_g_a(&s_D_g_a);
  if(*o_err)
    {
      write_warning("caa_marg_dens:Error calling re_makedata_g_a\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }

  *o_err = re_makedata_age1(*age_nBoats,*nAges,a_vec,
			    i_n_cov[0],i_age_int_nFac,i_ispat[0],
			    i_age_int_fix,i_age_int_c_cov,
			    i_num_adj_area,i_adj_area,
			    i_n_cov[1],i_age_hsz_nFac,i_ispat[1],
			    i_age_hsz_fix,i_age_hsz_c_cov,
			    i_num_adj_area,i_adj_area,
			    &s_D_age);
  if(*o_err)
    {
      write_warning("caa_marginal_density:Error calling re_makedata_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      return;
    }


  FREE(s_D_orig);



  #ifdef LOG_FILE
  fclose(g_caa_log);
  #endif

  return;    /* End of caa_marginal_density */
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

  /* Allocating memory for weight structs */
  err = alloc_lin(s_D_wgl,&s_weight);
  if(err)
    {
      write_warning("initialize:Error calling init_lin\n");
      return(err);
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

  /* Allocating memory for weight structs */
  err = re_alloc_lin(s_D_wgl,&s_weight);
  if(err)
    {
      write_warning("re_initialize:Error calling re_init_lin\n");
      return(err);
    }

  return(0);
}            /* End of re_initialize */



/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates log-likelihood for age-model
*/
int calc_lik_age(Age_struct *i_age,Data_age *i_D_age,double *o_loglik)
{
  int  h;
  double sum;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  sum = G_ZERO;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    sum += calc_lik_age_h(h,i_age,i_D_age);
  
  (*o_loglik) = sum;

  return(0);
}		/* end of calc_lik_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the age and
     lga models

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_model1(int i_it,Age_struct *i_age,Data_age *i_D_age,
		    LW_struct *i_lga,Data_lin *i_D_lga,
		    double *o_mean_inv_lik_mod1)
{
  int  err,h;
  double loglik_age,loglik_lga,lik;
  double *res;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  res = CALLOC(i_D_lga->glm->nxcov,double);     // Free ok

  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      loglik_age = calc_lik_age_h(h,i_age,i_D_age);
      loglik_lga = calc_lik_lin_h(h,i_lga,i_D_lga,res);
      lik = exp(loglik_age+loglik_lga);
      err = update_mean(&o_mean_inv_lik_mod1[h],G_ONE/lik,i_it);
      #ifdef DEBUG_EVALUATE_FILE
      fprintf(s_unit,"%d %d %lf\n",i_it,h,lik);
      #endif
    }
  
  FREE(res);

  return(0);
}		/* end of Bayes_CV_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the age model

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_age(int i_it,Age_struct *i_age,Data_age *i_D_age,
                 double *o_mean_inv_lik_age)
{
  int  err,h;
  double lik;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      lik = exp(calc_lik_age_h(h,i_age,i_D_age));
      err = update_mean(&o_mean_inv_lik_age[h],G_ONE/lik,i_it);
      #ifdef DEBUG_EVALUATE_FILE
      fprintf(s_unit,"%d %d %lf\n",i_it,h,lik);
      #endif
    }
  

  return(0);
}		/* end of Bayes_CV_age */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates the log-likelihood for the age model in a specific haul
*/
double calc_lik_age_h(int i_h,Age_struct *i_age,Data_age *i_D_age)
{
  int  a;
  double term1,term2,loglik, n_h;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  n_h = G_ZERO;
  term1 = G_ZERO;
  term2 = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      n_h += (double) i_D_age->Ages_fix[i_h][a];
      term1 += (double) i_D_age->Ages_fix[i_h][a] *i_age->alpha[i_h][a];
      term2 += exp(i_age->alpha[i_h][a]);
    }
  term2 = n_h * log(term2);
  loglik = term1-term2;
  
  return(loglik);
}		/* end of calc_lik_age_h */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates log-likelihood for the lga or wgl model
*/
int calc_lik_lin(LW_struct *i_lin,Data_lin *i_D_lin,double *o_loglik)
{
  int  h;
  double loglik, *res;

  res = CALLOC(i_D_lin->glm->nxcov,double);     // Free ok

  loglik = G_ZERO;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    loglik += calc_lik_lin_h(h,i_lin,i_D_lin,res);
  
  (*o_loglik) = loglik;

  FREE(res);

  return(0);
}		/* end of calc_lik_lin */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Update the mean of the inverse likelihood for the lga or wgl model

  The mean of the inverse likelihood can be used to calculate the 
  Pseudo (cross-validated) Bayesian factor.
*/
int Bayes_CV_lin(int i_it,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *o_mean_inv_lik_lin)
{
  int  err,h;
  double lik, *res;

  res = CALLOC(i_D_lin->glm->nxcov,double);     // Free ok

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_lin er en struct som har data estimater paa forskjellige 
     parametere */

  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      lik = exp(calc_lik_lin_h(h,i_lin,i_D_lin,res));
      err = update_mean(&o_mean_inv_lik_lin[h],G_ONE/lik,i_it);
    }
  

  FREE(res);

  return(0);
}		/* end of Bayes_CV_lin */


/*!
  \author Geir Storvik and Baard Storvik
  \brief Calculates the log-likelihood for the lga or wgl model in a specific haul
*/
double calc_lik_lin_h(int i_h,LW_struct *i_lin,Data_lin *i_D_lin,double *w_res)
{
  int  i,j;
  double term12,sum, N,l_h,loglik;
  Data_cov *xcov;

  /* LW_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_lin er en struct som har data estimater paa forskjellige 
     parametere */
  /* Se write_it_lin for forklaring på i_lin
     Se fish_sim.h for forklaring på Data_lin */
  N = (double) (i_D_lin->glm->suff[i_h][0][0]);

  /* Number of fish in haul times Number of nHaul */
  /* tau_cell,tau_area,tau_haul,tau_fish */
  term12 = -G_HALF * N*(log(G_TWO*G_PI)-log(i_lin->par->tau_obs));
  sum = G_ZERO;
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      xcov = i_D_lin->glm->xcov[i];
      w_res[i] = calc_eff(xcov,i_lin->par->eff[0][i],i_h) - 
	i_D_lin->glm->beta_hat[i_h][0][i];
    }

  l_h = i_D_lin->glm->ssq[i_h];
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    for(j=0;j<i_D_lin->glm->nxcov;j++)
    {
      l_h += i_D_lin->glm->suff[i_h][i][j] * w_res[i] * w_res[j];
    }
  sum += l_h;
  
  loglik = term12-G_HALF*sum*i_lin->par->tau_obs;

  return(loglik);
}		/* end of calc_lik_lin_h */


/*!
  \author Geir Storvik
  \brief Calculates residuals in lga model based on current simulated parameters
*/
int calc_resid_lga(int *i_totage, double *i_totlength,int *i_nFishBoat,
		   int *i_replength,
                   Age_struct *i_age, Data_age *i_D_age,LW_struct *i_lin,Data_lin *i_D_lin,
		   Data_g_a *i_D_g_a,double *o_resid)
{
  int  a,f,h,i,ind;
  double *mu, *beta;
  Data_cov *xcov;

  beta = CALLOC(i_D_lin->glm->nxcov,double);
  mu = CALLOC(i_D_age->glm->ncat,double);

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      for(i=0;i<i_D_lin->glm->nxcov;i++)
	{
	  xcov = i_D_lin->glm->xcov[i];
	  beta[i] = calc_eff(xcov,i_lin->par->eff[0][i],h);
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]];
      for(f=0;f<i_nFishBoat[h];f++)
	{
	  if(i_totage[ind] > -1000 && i_totlength[ind]> -1000.00)
	    o_resid[ind] = i_totlength[ind]-mu[i_totage[ind]-i_D_age->a_vec[0]];
	  else
	    o_resid[ind] = -99999.99;
	  ind++;
	}
    }
  #ifdef DEBUG_EVALUATE
  for(i=0;i<min(ind,100);i++)
    printf("i=%d,a=%d,res_lga=%lf\n",i,(int) i_totage[i],o_resid[i]);
  #endif
  FREE(beta);
  FREE(mu);

  return(0);
}		/* end of calc_resid_lga */



/*!
  \author Geir Storvik
  \brief Calculates residuals in wgl model based on current simulated parameters
*/
int calc_resid_wgl(double *i_totlength, double *i_totweight,int *i_nFishBoat,
                   LW_struct *i_weight,Data_lin *i_D_wgl,
		   double *o_resid)
{
  int  f,h,i,ind;
  double mu, *beta;
  Data_cov *xcov;

  beta = CALLOC(i_D_wgl->glm->nxcov,double);

  ind = 0;
  for(h=0;h<i_D_wgl->glm->nHaul;h++)
    {
      for(i=0;i<i_D_wgl->glm->nxcov;i++)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  beta[i] = calc_eff(xcov,i_weight->par->eff[0][i],h);
	}
      for(f=0;f<i_nFishBoat[h];f++)
	{
	  if(i_totlength[ind]> -1000.00 && i_totweight[ind]> -1000.00)
	    {
	      mu = beta[0]+beta[1]*i_totlength[ind];
	      o_resid[ind] = i_totweight[ind]-mu;
	    }
	  else
	    o_resid[ind] = min(i_totlength[ind],i_totweight[ind]);
	  ind++;
	}
    }
  FREE(beta);

  return(0);
}		/* end of calc_resid_wgl */





/*!
  \author Geir Storvik 
  \brief Calculates the Kolmogorov-Smirnov statistic of haul-randomeffects of age
*/
int calc_KS_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_p)
{
  int    a,h,i,j,n,sum;
  long   ind;
  double alpha,d,d_rep,prob,sqrt_tau;
  //FILE   *unit;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  //unit = fopen("alpha_h.txt","w");

  ind = 1;
  sqrt_tau = sqrt(i_age->par->tau_obs);
  for(h=0;h<i_D_age->glm->nHaul;h++)
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      alpha = calc_eff(i_D_age->glm->xcov[0],i_age->par->eff[a][0],h);
      s_eps[ind] = (i_age->alpha[h][a] - alpha)*sqrt_tau;
      //fprintf(s_unit2,"%lf ",s_eps[ind]);
      //fprintf(unit,"%d %d %lf\n",h,a,s_eps[ind]);
      ind++;
    }
  //fclose(unit);
  n = ind-1;
  //fprintf(s_unit2,"\n");

  ksone(s_eps,n,pnorm0,&d,&prob,&ind);
  #ifdef DEBUG_IND_KS
  fprintf(s_ind_ks,"%ld\n",ind);
  #endif
  
  if(0)
    {
  sum = 0;
  for(j=0;j<s_N2;j++)
    {
      for(i=1;i<=n;i++)
	s_eps[i] = gennor(G_ZERO,G_ONE);
      ksone(s_eps,n,pnorm0,&d_rep,&prob,&ind);
      sum += (d_rep > d);
    }
  *o_p = (double) sum/(double) s_N2;
    }
  //fprintf(s_unit2,"%lf %lf\n",d,d_rep);

  *o_d = d;
  *o_p = prob;

  return(0);
}		/* end of calc_KS_age */


/*!
  \author Geir Storvik 
  \brief Calculates the Robins-etal Discr measure of haul-random effects of age
*/
int calc_D_Robins_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_p)
{
  int    a,h,ind,k1,k2,i,j,n,sum,x,y;
  double alpha,d,d_rep,sqrt_tau;
  static int (*cmp)();

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */


  ind = 0;
  sqrt_tau = sqrt(i_age->par->tau_obs);
  for(h=0;h<i_D_age->glm->nHaul;h++)
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      alpha = calc_eff(i_D_age->glm->xcov[0],i_age->par->eff[a][0],h);
      s_eps[ind] = (i_age->alpha[h][a] - alpha)*sqrt_tau;
      //fprintf(s_unit2,"%lf ",s_eps[ind]);
      ind++;
    }
  n = ind;
  //fprintf(s_unit2,"\n");

  cmp = compd;
  k1 = (int) ((double) n * 0.1);
  k2 = n-k1+1;

  qsort(s_eps,n,sizeof(double),cmp);
  d = fabs(fabs(s_eps[k2]) - fabs(s_eps[k1]));
  *o_d = d;

  x = (int) (s_N2*d/s_d_rep[s_N2-1]+G_ONE);
  y = s_g[x]+1;
  while(s_d_rep[y-1]>d)
    y--;
  *o_p = (double) (y+1)/(double) s_N2;
  if(0)
    {
  sum = 0;
  for(j=0;j<s_N2;j++)
    {
      for(i=0;i<n;i++)
	s_eps[i] = gennor(G_ZERO,G_ONE);
      qsort(s_eps,n,sizeof(double),cmp);
      d_rep = fabs(fabs(s_eps[k2]) - fabs(s_eps[k1]));
      sum += (d_rep > d);
    }
  *o_p = (double) sum/(double) s_N2;
    }
  

  return(0);
}		/* end of calc_D_Robbins_age */



/*!
  \author Geir Storvik 
  \brief Calculates the Kolmogorov-Smirnov statistic of haul-randomeffects of age
*/
int calc_entropy_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_d2,
		     double *o_p)
{
  int    a,h,j,k1,k2;
  double alpha,sum[2],sum1,sum2,d_min,d_max,d_mean,sd;
  double d_min_rep,d_max_rep,d_mean_rep;
  static int (*cmp)();
  //FILE   *unit;

  /* Age_struct er en struct som inneholder parameterne i likelihood*/
  /* Data_age er en struct som har data estimater paa forskjellige 
     parametere */

  cmp = compd;
  k1 = (int) (((double) i_D_age->glm->nHaul)*0.25);
  k2 = (int) (((double) i_D_age->glm->nHaul)*0.75);
  k1 = 0;
  k2 = i_D_age->glm->nHaul-1;

  d_mean = G_ZERO;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      sum1 = G_ZERO;
      sum2 = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  sum1 += exp(i_age->alpha[h][a]);
	  sum2 += exp(i_age->alpha[h][a])*i_age->alpha[h][a];
	}
      s_entr[h] = log(sum1) - sum2/sum1;
      d_mean += s_entr[h];
    }
  qsort(s_entr,i_D_age->glm->nHaul,sizeof(double),cmp);
  d_min = s_entr[k1];
  d_max = s_entr[k2];
  d_mean /= (double) i_D_age->glm->nHaul;
  *o_d = d_min;
  *o_d2 = d_max;
  //*o_d2 = d_mean;

  sum[0]= G_ZERO;
  sum[1]= G_ZERO;
  sd = 1/sqrt(i_age->par->tau_obs);
  //fprintf(stderr,"sd=%lf\n",sd);
  for(j=0;j<s_N2;j++)
    {
      d_max_rep = G_ZERO;
      d_min_rep = log((double) i_D_age->glm->ncat);
      d_mean_rep = G_ZERO;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  sum1 = G_ZERO;
	  sum2 = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      alpha = calc_eff_no_haul(i_D_age->glm->xcov[0],i_age->par->eff[a][0],h);
	      alpha += sd*gennor(G_ZERO,G_ONE);
	      sum1 += exp(alpha);
	      sum2 += exp(alpha)*alpha;
	    }
	  s_entr[h] = log(sum1) - sum2/sum1;
	  //d_min_rep = MIN(d_min_rep,entr);
	  //d_max_rep = MAX(d_max_rep,entr);
          d_mean_rep += s_entr[h];
	}
      qsort(s_entr,i_D_age->glm->nHaul,sizeof(double),cmp);
      d_min_rep = s_entr[k1];
      d_max_rep = s_entr[k2];
      d_mean_rep /= (double) i_D_age->glm->nHaul;
      sum[0] += (d_min_rep > d_min);
      sum[1] += (d_max_rep > d_max);
      //sum[1] += (d_mean_rep > d_mean);
    }
  o_p[0] = (double) sum[0]/(double) s_N2;
  o_p[1] = (double) sum[1]/(double) s_N2;
  //fprintf(stderr,"d= %lf d_rep = %lf\n",d_min,d_min_rep);

  return(0);
}		/* end of calc_entropy_age */

