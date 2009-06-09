/*!
  \file caa_sample_multi.c
  \author Geir Storvik
  \brief Routines for sampling nonlinear part of age model

  This file contains a lot of routines that currently is not used but have
  been tried out in order to obtain better convergence and/or easier implementation
*/
#include "caa.h"


static GMRFLib_hidden_param_tp *s_hidden_par=NULL;
static GMRFLib_optimize_param_tp *s_optimize_param;
static int Loglik_poisson_age(double *logll,double *x,
			      int m,int idx,double *x_vec,char *arg);
static int sample_age_alpha_ages_given(int i_h,Age_struct *i_age,Data_age *i_D_age,
				      int i_acc,int *o_acc);
static double age_haul_calc_posterior(int i_ncat,double *i_alpha,double *i_prob,
				      int *i_Ages,double *i_mu,double i_tau);
static int age_haul_find_mode(int i_ncat,int *i_Ages,double i_N,
			      double *i_mu,double i_tau,
			      double *x_alpha_opt,double *o_prob,double *o_log_opt,
			      double *o_grad,double **o_Hess);
static int age_haul_calc_prob(int i_ncat,double *i_alpha,double *o_prob);


static double  *s_mu;        /*!< Prior means for alpha's */
static double  *s_prob;      /*!< Probabilities for multinomial distribution */
static double  *s_alpha_opt; /*!< Optimal values of alpha's */
static double  *s_alpha_prop;/*!< Proposal values of alpha's */
static double  *s_grad;      /*!< Gradient vector for optimization of haul likelihood */
static double **s_Hess;      /*!< Hesse matrix for optimization of haul likelihood */ 
static double  *s_eps;       /*!< Working vector for sampling alpha's */
static double  *s_eps2;      /*!< Working vector for sampling alpha's */
static double  *s_delta;     /*!< Change in Newton-Raphson for finding mode */
static int      s_N_gauher;  /*!< Number of nodes in Gauss Hermite quadrature */
static double  *s_gauher_x;  /*!< Abscissas in Gauss Hermite quadrature */
static double  *s_gauher_w;  /*!< Weights in Gauss Hermite quadrature */

#define NGIBBS  10        /*! Number of Gibbs iterations for sampling alpha's */

//#define AGE_TEST 1

#ifdef LOG_FILE
extern FILE     *g_caa_log;
extern FILE     *g_caa_alpha;
#endif

/*!
  \author Geir Storvik
  \brief Allocates space for routines used to sample non-linear part of age-model

  This routine must be called before using routines in the caa_sample_multi.c file.
*/
int sample_multi_initialize(int i_ncat)
{
  int    i;
  double sum;

  if(0){
  GMRFLib_default_optimize_param(&s_optimize_param);

  /* sett default verdier */
  if (!s_hidden_par) GMRFLib_default_hidden_par(&s_hidden_par);
  s_hidden_par->cmeanmode   = GMRFLib_COND_MEAN;
  s_hidden_par->gaussapprox = 0;
  /* er denne lik 1, faas approximasjon A1 */
  s_hidden_par->modeoption  = GMRFLib_MODEOPTION_MODE;
  s_hidden_par->nantithetic = 4;
  s_hidden_par->neightype   = GMRFLib_NEIGHTYPE_GRAPH;
  s_hidden_par->norder      = 2;
  s_hidden_par->range       = 6.;
  s_hidden_par->nresolution = 10;
  /* evnt sett lik 10, = anntall regioner i spline */
  s_hidden_par->neighpar    = 1;
  /* 0 betyr ingen korreksjon for integral-leddet.
     1 betyr korreksjon for naboene i grafen, 
     2 naboene og naboenes nabo, etc...   */

  s_hidden_par->nsample     = 10;
  /* hvor mange sample som skal brukes.  */
  if(0)
    {/* approximasjon A2, finnes ved */
      s_hidden_par->neighpar = 0;
      s_hidden_par->nsample  = 0;
    }
  if(1)
    {/* approximasjon A3, billig	*/
      s_hidden_par->neighpar    = 2;
      s_hidden_par->nsample     = 1; /* kun forventningen brukes */
      s_hidden_par->nantithetic = 0; /* og da maa denne vaere null */
    }
  if(0)
    {/* approximasjon A3, bedre*/
      s_hidden_par->nantithetic = 4;  /* default verdi */
      s_hidden_par->neighpar    = 2;
      s_hidden_par->nsample     = 10;  
    }
  }
  s_mu = CALLOC(i_ncat,double);        // Free ok
  if(!s_mu)
    {
      write_warning("sample_multi_initialize:Error allocating s_mu\n");
      return(1);
    }
  s_prob = CALLOC(i_ncat,double);      // Free ok
  if(!s_prob)
    {
      write_warning("sample_multi_initialize:Error allocating s_prob\n");
      return(1);
    }
  s_alpha_opt = CALLOC(i_ncat,double); // Free ok
  if(!s_alpha_opt)
    {
      write_warning("sample_multi_initialize:Error allocating s_alpha_opt\n");
      return(1);
    }
  s_alpha_prop = CALLOC(i_ncat,double);// Free ok
  if(!s_alpha_prop)
    {
      write_warning("sample_multi_initialize:Error allocating s_alpha_prop\n");
      return(1);
    }
  s_grad = CALLOC(i_ncat,double);     // Free ok
  if(!s_grad)
    {
      write_warning("sample_multi_initialize:Error allocating s_grad\n");
      return(1);
    }
  s_Hess = Mmatrix_2d(0,i_ncat-1,0,i_ncat-1,sizeof(double),1); // Free ok
  if(!s_Hess)
    {
      write_warning("sample_multi_initialize:Error allocating s_Hess\n");
      return(1);
    }
  s_eps = CALLOC(i_ncat,double);      // Free ok
  if(!s_eps)
    {
      write_warning("sample_multi_initialize:Error allocating s_eps\n");
      return(1);
    }
  s_eps2 = CALLOC(i_ncat,double);     // Free ok
  if(!s_eps2)
    {
      write_warning("sample_multi_initialize:Error allocating s_eps2\n");
      return(1);
    }
  s_delta = CALLOC(i_ncat,double);    // Free ok
  if(!s_delta)
    {
      write_warning("sample_multi_initialize:Error allocating s_delta\n");
      return(1);
    }

  //Gauss hermite weights
  s_N_gauher = 30;
  s_gauher_x = CALLOC(s_N_gauher+1,double);  //Free ok
  s_gauher_w = CALLOC(s_N_gauher+1,double);  //Free ok
  gauher(s_gauher_x,s_gauher_w,s_N_gauher);

  if(0) 
    {
      sum = G_ZERO;
      for(i=1;i<=s_N_gauher;i++)
	sum += s_gauher_w[i];
      for(i=1;i<=s_N_gauher;i++)
	s_gauher_w[i] /= sum;
    }


  #ifdef DEBUG_MULTI_ALPHA
  g_caa_alpha = fopen("multi_alpha.dat","w");
  #endif

  return(0);
}		/* end of sample_multi_initialize */

/*!
  \author Geir Storvik
  \brief Reallocate space allocated by sample_multi_initialize
*/
int sample_multi_re_initialize()
{
  FREE(s_mu);
  FREE(s_prob);
  FREE(s_alpha_opt);
  FREE(s_alpha_prop);
  FREE(s_grad);
  Fmatrix_2d(&s_Hess[0][0],&s_Hess[0]);
  FREE(s_eps);
  FREE(s_eps2);
  FREE(s_delta);
  FREE(s_gauher_x);
  FREE(s_gauher_w);

  #ifdef DEBUG_MULTI_ALPHA
  fclose(g_caa_alpha);
  #endif
  return(0);
}		/* end of sample_multi_re_initialize */



/*!
  \author Geir Storvik
  \brief Find optimal haul values for h's given
*/
int age_haul_modes(int i_start_h,int i_stop_h,Age_struct *i_age,Data_age *i_D_age)
{
  int     a,err,h;
  double  N_h,tau;
  double  log_opt;

  tau = i_age->par->tau_obs;

  /* Find mode */
  for(h=i_start_h;h<i_stop_h;h++)
    {
      N_h = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  N_h += (double) i_D_age->Ages[h][a];
	  s_mu[a] = G_ZERO;
	  i_age->alpha[h][a] = log((double) i_D_age->Ages[h][a] + 0.1);
	}
      if(N_h > 0)
	{

	  err = age_haul_find_mode(i_D_age->glm->ncat,i_D_age->Ages[h],N_h,s_mu,tau,
				   i_age->alpha[h],s_prob,&log_opt,s_grad,s_Hess);	
	  if(err)
	    {
	      write_warning("age_haul_modes:Error calling age_haul_find_mode\n");
	      return(err);
	    }
	}
    }

          
  return(0);
}		/* end of age_haul_modes */


/*!
  \author Geir Storvik
  \brief Sample ages based on empirical age-given-length distribution

  Assumes lengths are given by categories. For each length category posterior 
  probabilites are calculated using prior and length, and ages are sampled 
  from multinomial distribution.

  Stored in ages per haul.

  This routine is only used for initialization of the parameters to be simulated.

  \todo This is a new version of sample_ages_len_only_init to be used when
  only amigo-type data is available. This is a companion of 
  ::sample_ages_len_only_new, both wich should be the only options when cleaned
  up.
*/
int sample_ages_len_only_init_new(Data_orig *i_D_orig,Age_struct *i_age,
				  Data_age *i_D_age,LW_struct *i_length,
				  Data_lin *i_D_lga,Data_g_a *i_D_g_a,int saveSim)
{
  int            a,a2,f,h,i,s,l_int,ind_a,ind_f,N_int,season,nSeason;
  long          *ages;
  double       **P_al,***P_al_s;
  double         u,sum,sigma,lobs;

  int printHaul=0;
  int printHaulN=442;

  FILE *fp;
  if(saveSim)
    fp = fopen("ages_miss_start.txt","w");

  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  nSeason = i_D_g_a->nSeason;

  /* Initialize */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	}
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
          i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
          i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	}
    }
     

  /* First making a transition matrix P(a|l) for a finite number of lengths */
  N_int = i_D_orig->n_int_len;
  /* Data combined for all seasons, for use if no observations */
  P_al = Mmatrix_2d(0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);
  /* Data for each season */
  P_al_s = Mmatrix_3d(0,nSeason-1,0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok

  for(i=0;i<N_int+1;i++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      P_al[i][a] = G_ZERO;

  for(s=0;s<nSeason;s++)
    for(i=0;i<N_int+1;i++)
      for(a=0;a<i_D_age->glm->ncat;a++)
	P_al_s[s][i][a] = G_ZERO;

  ind_f = 0;  
  for(h=0;h<i_D_age->glm->nHaul;h++)
    for(f=0;f<i_D_orig->nFishBoat[h];f++)
      {
	a = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	lobs = i_D_orig->totlength[ind_f];
	season = i_D_orig->season[h];
	if(a > -1000 && lobs > -1000.0)
	  {
	    l_int = 0;
	    //while(lobs > i_D_orig->int_len_lim[l_int+1])
	    while(lobs > i_D_orig->int_len_lim[l_int])
	      l_int++;
	    P_al[l_int][a] += (double) i_D_orig->replength[ind_f];
	    P_al_s[season-1][l_int][a] += (double) i_D_orig->replength[ind_f];
	  }
	ind_f++;
      }

  // Convert to probabilities
  for(i=0;i<N_int+1;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	sum += P_al[i][a];
      if(sum < 0.0001)
	{
	  if(i==0)
	    P_al[i][0] = G_ONE;
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al[i][a] = P_al[i-1][a];
	    }
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al[i][a] /= sum;
	}
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int+1;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_s[s][i][a];
	  if(sum<0.0001)//use P_al (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] = P_al[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] /= sum;
	    }
	}
    }


  // TESTING
  // Convert to cummulative probabilities
  if(0)
    {
      for(i=0;i<N_int+1;i++)
	{
	  for(a=1;a<i_D_age->glm->ncat;a++)
	    P_al[i][a] += P_al[i][a-1];
	  P_al[i][i_D_age->glm->ncat-1] *= 1.00000002; // to avoid problems in end
	}
    }
  // Start simulation of missing ages
  ind_f = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  a = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[ind_f];
	  if(a < -1000 && lobs > -1000.0)
	    {
	      //Find right interval
	      l_int = 0;
	      //while(lobs > i_D_orig->int_len_lim[l_int+1])
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      //Simulate age
	      //TESTING
	      if(0)
		{
		  u = genunf(G_ZERO,G_ONE);
		  a = 0;
		  while(a<(i_D_age->glm->ncat-1) && u > P_al[l_int][a])
		    a++;
		  i_D_age->Ages[h][a]++;
		  //i_D_lga->Ages[h][a]++; //Must correct index a
		  //i_D_lga->sum_by_cat[h][a] += lobs;
		  //i_D_lga->sqsum_by_cat[h][a] += lobs*lobs;
		}
	      my_genmul(i_D_orig->replength[ind_f],P_al_s[season-1][l_int],i_D_age->glm->ncat,ages);
	      for(a2=0;a2<i_D_age->glm->ncat;a2++)
		{
		  if(ages[a2]>0)
		    {
		      i_D_age->Ages[h][a2] += (int) ages[a2];
		      ind_a = i_D_g_a->a2Age_vec[a2]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a2];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a2];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a2];
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a2,season,
				i_D_g_a->a_vec[ind_a],lobs,(int) ages[a2]);
		      if(printHaul)
			if(h==printHaulN)
			  printf("%d %d %d %f %f %d\n",h,a2,season,
				 i_D_g_a->a_vec[ind_a],lobs,(int) ages[a2]);
		    }
		}
	    }
	  ind_f++;
	}
    }

  if(saveSim)
    fclose(fp);

  // Free allocated memory
  FREE(ages);
  Fmatrix_2d(&P_al[0][0],&P_al[0]);
  Fmatrix_3d(&P_al_s[0][0][0],&P_al_s[0][0],&P_al_s[0]);

  return(0);
}		/* end of sample_ages_len_only_init_new */


/*!
  \author Geir Storvik
  \brief Samples missing ages.

  Here all data are assumed to be of the amigo type, that is long strings of
  age and length. 

  In order to speed up computation, the fish are assumed ordered in hauls. 
  Further, inside each haul, the lenghts are assumed ordered in increasing values. 
  Then the length-range is divided into a finite
  number of intervals in which the age-probabilities are assumed constant for
  all length-values inside an interval.
*/
int sample_ages_len_only_new(Data_orig *i_D_orig,Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
			     int saveSim, int i_it)
{
  int            a,f,h,i,ind_a,ind_f,ind_miss,n,cum_fish,aobs,season;
  int            l_int,l_int_prev;
  double         lobs,r;
  double         sum_p,sigma;
  double         lstart,lend;
  Data_glm      *glm;
  double        *p, *p2, *mu, *beta;
  long          *ages;
  long          *test_ages, *test_ages2;
  

  glm = i_D_lga->glm;

  p = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  p2 = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);      // Free ok

  test_ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  test_ages2 = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  beta = CALLOC(i_D_lga->glm->nxcov,double);      // Free ok      

  sigma = G_ONE/sqrt(i_length->par->tau_obs);

  #ifdef DEBUG_TEST
  FILE *unit;
  unit = fopen("Ages_samples.dat","w");
  #endif
  
  FILE *fp;
  if(saveSim)
    fp = fopen("ages_miss.txt","w");

  int printHaul=0;
  int printHaulN=442;

  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      test_ages[a] = 0;
      test_ages2[a] = 0;
    }
  
  cum_fish = 0;
  ind_miss = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];

      /* Start by initializing sufficient statistics */ 
      if(i_D_age->type_age[0]==1) /* Age errors: All simulated */
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      i_D_age->Ages[h][a] = 0;
	      i_D_age->Ages_disc[h][a] = 0;
	      i_D_age->Ages_land[h][a] = 0;
	    }
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      i_D_lga->Ages[h][a] = 0;
	      i_D_lga->sum_by_cat[h][a] = G_ZERO;
	      i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	    }
	}
      else  /* No age errors: Keep aged fish */
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	      i_D_age->Ages_disc[h][a] = 0;
	      i_D_age->Ages_land[h][a] = 0;
	    }
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
	      i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
	      i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	    }
	}

      /* find prior age-probabilities */
      sum_p = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] = exp(i_age->alpha[h][a]);
	  sum_p += p[a];
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] /= sum_p;
	}
      /* Find intercept and slope for lga model */
      for(i=0;i<i_D_lga->glm->nxcov;i++)
	beta[i] = calc_eff(i_D_lga->glm->xcov[i],i_length->par->eff[0][i],h);
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
	}
      if(i_D_age->type_age[0]==1)/*Ages observed with error*/
	{
	  /* Both age and length observed */
	  for(f=0;f<(i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h]);f++)
	    {
	      ind_f = cum_fish+f;
	      aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	      lobs = i_D_orig->totlength[ind_f];
	      
	      /* Calculate probabilities */
	      sum_p = G_ZERO;
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		{
		  a = i_age->A_neigh[aobs][i];
		  p2[i] = p[a]*i_age->A2A[a][aobs]*dnorm(lobs,mu[a],sigma);
		  sum_p += p2[i];
		}
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		p2[i] /= sum_p;
	      my_genmul(i_D_orig->replength[ind_f],p2,i_age->A_Nneigh[aobs],ages);
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		{
		  a = i_age->A_neigh[aobs][i];
		  i_D_age->Ages[h][a] += (int) ages[i];
		  
		  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		  i_D_lga->Ages[h][ind_a] += (int) ages[i];
		  test_ages2[a] += ages[i];
		  /* Update sufficient statistics */
		  r = (double) ages[i];
		  if(r< 0)
		    {
		      write_warning("Negative ages generated");
		      return(1);
		    }
		  i_D_lga->sum_by_cat[h][ind_a] += r*lobs;
		  i_D_lga->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
		}
	    }
	}

      ind_f = i_D_orig->start_noAge[h];
      l_int_prev = -1;
      l_int = 0;
      /* Loop through all non-aged fish */
      for(f=0;f<i_D_orig->num_noAge[h];f++)
	{
	  if(i_D_orig->replength[ind_f+f]>0)
	    {
	      lobs = i_D_orig->totlength[ind_f+f];
	      lstart = log(round(exp(lobs))-0.5);
	      lend = log(round(exp(lobs))+0.5);
	      //     while(lobs > i_D_orig->int_len_lim[l_int+1])
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      if(l_int!=l_int_prev)
		{
		  //New calculation of age-probabilities needs to be performed
		  /* Calculate probabilities */
		  sum_p = G_ZERO;
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    {
		      //p2[a] = p[a]*dnorm(lobs,mu[a],sigma);
		      p2[a] = p[a]*(pnorm(lend,mu[a],sigma)-pnorm(lstart,mu[a],sigma));
		      sum_p += p2[a];
		    }
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    p2[a] /= sum_p;
		  l_int_prev = l_int;
		}
	      my_genmul(i_D_orig->replength[ind_f+f],p2,i_D_age->glm->ncat,ages);
	      n = 0;
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  if(ages[a]>0)
		    {
		      i_D_age->Ages[h][a] += (int) ages[a];
		      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
		      test_ages[a] += ages[a];
		      n+= ages[a];
		      if(printHaul)
			if(h==printHaulN)
			  printf("f=%d,a=%d,season=%d,l=%f,Nages=%d\n",ind_f+f,a,season,lobs,(int) ages[a]);
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a]);
		    }

              #ifdef DEBUG_TEST
		  fprintf(unit,"%d %d %d %lf\n",h,f,a,lobs);
              #endif
		}
	    }
	}
      ind_miss += i_D_g_a->ncat*(i_D_orig->n_int_len+1);
      cum_fish += i_D_orig->nFishBoat[h];
    }
  if(saveSim)
    fclose(fp);

  #ifdef AGE_TEST
  fprintf(stderr,"\nobs:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(stderr,"%5d ",test_ages2[a]);
  fprintf(stderr,"\nunobs: ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(stderr,"%5d ",test_ages[a]);
  FILE *ageFile;
  ageFile = fopen("../Ages_sampled.dat","a");
  fprintf(ageFile,"\nobs:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(ageFile,"%3d:%5d ",a+1,test_ages2[a]);
  fprintf(ageFile,"\nunobs: ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(ageFile,"%3d:%5d ",a+1,test_ages[a]);
  double testSum;
  fprintf(ageFile,"\nsum:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      testSum = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  //testSum += i_D_age->Ages[h][a];
	  testSum += i_D_lga->Ages_fix[h][a];
	}
      fprintf(ageFile,"%3d:%5.0f ",a+1,testSum);
    }
  fclose(ageFile);
  #endif 

  #ifdef DEBUG_TEST
  fclose(unit);
  #endif


  // Free memory allocated in this routine
  FREE(p);
  FREE(p2);
  FREE(mu);
  FREE(ages);
  FREE(test_ages);
  FREE(test_ages2);
  FREE(beta);

  return(0);
}		/* end of sample_ages_len_only_new */



/*!
  \author Hanne Rognebakke
  \brief Sample ages based on empirical age-given-length distribution

  Assumes lengths are given by categories. For each length category posterior 
  probabilites are calculated using prior and length, and ages are sampled 
  from multinomial distribution.

  Stored in ages per haul.

  This routine is only used for initialization of the parameters to be simulated.
*/
int sample_ages_len_only_init_new_CC(Data_orig *i_D_orig,Data_CC *i_D_CC,
				     Age_struct *i_age,Data_age *i_D_age,
				     Data_lin *i_D_lga,Data_g_a *i_D_g_a,
				     Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,int printSim)
{
  int            a,aobs,ncat_age,f,h,i,s,l_int,ind_f,ind_a,N_int,season;
  int            n_CC,n_Skrei,nSeason;
  long          *ages;
  double       **P_al,***P_al_s;
  double         u,sum,lobs;
  double         p,ptrue,pobs;

  FILE *fp;
  if(printSim)
    fp = fopen("ages_miss_start.txt","w");

  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  ncat_age = (int) i_D_age->glm->ncat/2;
  nSeason = i_D_g_a->nSeason;

  /* Initialize */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
          i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
          i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	}
      for(a=0;a<i_D_g_a_CC->ncat;a++)
	{
	  i_D_lga_CC->Ages[h][a] = i_D_lga_CC->Ages_fix[h][a];
          i_D_lga_CC->sum_by_cat[h][a] = i_D_lga_CC->sum_by_cat_fix[h][a];
          i_D_lga_CC->sqsum_by_cat[h][a] = i_D_lga_CC->sqsum_by_cat_fix[h][a];
	}
    }
  /* First making a transition matrix P(a|l) for a finite number of lengths */
  /* Using unceratain observations as certain for starting values even if classification error */
  N_int = i_D_orig->n_int_len;
  /* Data combined for all seasons, for use if no observations */
  P_al = Mmatrix_2d(0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);
  /* Data for each season */
  P_al_s = Mmatrix_3d(0,nSeason-1,0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok

  for(i=0;i<N_int+1;i++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      P_al[i][a] = G_ZERO;

  for(s=0;s<nSeason;s++)
    for(i=0;i<N_int+1;i++)
      for(a=0;a<i_D_age->glm->ncat;a++)
	P_al_s[s][i][a] = G_ZERO;

  ind_f = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[ind_f];
	  season = i_D_orig->season[h];
	  if(aobs > -1000 && lobs > -1000.0)/*Both age and length observed*/
	    {
	      l_int = 0;
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      if(i_D_orig->tottype[ind_f] == 1) //certain coastal cod
		{
		  P_al[l_int][aobs] += (double) i_D_orig->replength[ind_f];
		  P_al_s[season-1][l_int][aobs] += (double) i_D_orig->replength[ind_f];
		}
	      else if(i_D_orig->tottype[ind_f] == 2) //uncertain coastal cod
		{
		  P_al[l_int][aobs] += (double) i_D_orig->replength[ind_f];
		  P_al_s[season-1][l_int][aobs] += (double) i_D_orig->replength[ind_f];
		}
	      else if(i_D_orig->tottype[ind_f] == 4) //uncertain skrei
		{
		  P_al[l_int][ncat_age+aobs] += (double) i_D_orig->replength[ind_f];
		  P_al_s[season-1][l_int][ncat_age+aobs] += (double) i_D_orig->replength[ind_f];
		}
	      else if(i_D_orig->tottype[ind_f] == 5) //certain skrei
		{
		  P_al[l_int][ncat_age+aobs] += (double) i_D_orig->replength[ind_f];
		  P_al_s[season-1][l_int][ncat_age+aobs] += (double) i_D_orig->replength[ind_f];
		}
	    }
	  ind_f++;
	}
    }
  // Convert to probabilities
  for(i=0;i<N_int+1;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	sum += P_al[i][a];
      if(sum < 0.0001)
	{
	  if(i==0)
	    {
	      P_al[i][0] = 0.5;
	      P_al[i][ncat_age] = 0.5;
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al[i][a] = P_al[i-1][a];
	    }
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al[i][a] /= sum;
	}      
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int+1;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_s[s][i][a];
	  if(sum<0.0001)//use P_al (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] = P_al[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] /= sum;
	    }
	}
    }

  // Start simulation of missing ages, or classification
  // Simulate only ages in year, since season is already known
  ind_f = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[ind_f];
	  season = i_D_orig->season[h];
	  //Find right interval
	  l_int = 0;
	  while(lobs > i_D_orig->int_len_lim[l_int])
	    l_int++;
	  if(aobs > -1000)/*Both age and length observed - add uncertain observations*/
	    {
	      ptrue = P_al[l_int][aobs];
	      n_CC = 0;
	      n_Skrei = 0;
	      if(i_D_orig->tottype[ind_f] == 2) //uncertain coastal cod
		{
		  pobs = i_D_CC->ptype2_CC[aobs];
		  if(pobs==1.0)
		    p = 1.0;
		  else if(ptrue==0.0)
		    p = 0.0;
		  else
		    p = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		  if(p<0 || p>1)
		    write_warning("sample_ages_len_only_init_new_CC: Something is wrong");
		  for(i=0; i<i_D_orig->replength[ind_f]; i++)
		    {
		      if(GMRFLib_uniform()<p)
			n_CC++;
		      else
			n_Skrei++;
		    }
		}
	      else if(i_D_orig->tottype[ind_f] == 4) //uncertain skrei
		{
		  pobs = i_D_CC->ptype4_S[aobs];
		  if(pobs==1.0)
		    p = 1.0;
		  else if(ptrue==0.0)
		    p = 0.0;
		  else
		    p = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		  if(p<0 || p>1)
		    write_warning("sample_ages_len_only_init_new_CC: Something is wrong");
		  for(i=0; i<i_D_orig->replength[ind_f]; i++)
		    {
		      if(GMRFLib_uniform()<p)
			n_Skrei++;
		      else
			n_CC++;
		    }
		}
	      //else (tottype ==3 or tottype == -99999)
	      if(n_CC > 0 || n_Skrei > 0)
		{
		  // Add coastal cod
		  i_D_age->Ages[h][aobs] += n_CC;
		  ind_a = i_D_g_a_CC->a2Age_vec[aobs]+(season-1);
		  i_D_lga_CC->Ages[h][ind_a] += n_CC;
		  i_D_lga_CC->sum_by_cat[h][ind_a] += lobs * n_CC;
		  i_D_lga_CC->sqsum_by_cat[h][ind_a] += lobs*lobs * n_CC;
		  // Add skrei
		  i_D_age->Ages[h][ncat_age+aobs] += n_Skrei;
		  ind_a = i_D_g_a->a2Age_vec[aobs]+(season-1);
		  i_D_lga->Ages[h][ind_a] += n_Skrei;
		  i_D_lga->sum_by_cat[h][ind_a] += lobs * n_Skrei;
		  i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * n_Skrei;
		  if(printSim)
		    {
		      if(n_CC>0)
			fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				h+1,i_D_g_a_CC->a_vec[ind_a],lobs,n_CC,aobs,season,i_D_orig->tottype[ind_f],0);
		      if(n_Skrei>0)
			fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				h+1,i_D_g_a->a_vec[ind_a],lobs,n_Skrei,ncat_age+aobs,season,i_D_orig->tottype[ind_f],0);
		    }
		}
	    }
	  if(aobs < -1000)
	    {
	      //Simulate age
	      //my_genmul(i_D_orig->replength[ind_f],P_al[l_int],i_D_age->glm->ncat,ages);
	      my_genmul(i_D_orig->replength[ind_f],P_al_s[season-1][l_int],i_D_age->glm->ncat,ages);
	      
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  if(ages[a]>0)
		    {
		      i_D_age->Ages[h][a] += (int) ages[a];
		      if(a<ncat_age) // Coastal cod
			{
			  ind_a = i_D_g_a_CC->a2Age_vec[a]+(season-1);
			  i_D_lga_CC->Ages[h][ind_a] += (int) ages[a];
			  i_D_lga_CC->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			  i_D_lga_CC->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			  if(printSim)
			    fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				    h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,1,1);
			}
		      else // Skrei
			{
			  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
			  i_D_lga->Ages[h][ind_a] += (int) ages[a];
			  i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			  i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			  if(printSim)
			    fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				    h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,5,1);
			}
		    }
		}
	    }
	  ind_f++;
	}
    }
      
  if(printSim)
    fclose(fp);

  
  // Free allocated memory
  FREE(ages);
  Fmatrix_2d(&P_al[0][0],&P_al[0]);
  Fmatrix_3d(&P_al_s[0][0][0],&P_al_s[0][0],&P_al_s[0]);

  return(0);
}		/* end of sample_ages_len_only_init_new_CC */


/*!
  \author Geir Storvik
  \brief Samples missing ages.

  Here all data are assumed to be of the amigo type, that is int strings of
  age and length. 

  In order to speed up computation, the fish are assumed ordered in hauls. 
  Further, inside each haul, the lenghts are assumed ordered in increasing values. 
  Then the length-range is divided into a finite
  number of intervals in which the age-probabilities are assumed constant for
  all length-values inside an interval.
*/
int sample_ages_len_only_new_CC(Data_orig *i_D_orig,Data_CC *i_D_CC,
				Age_struct *i_age,Data_age *i_D_age,
				LW_struct *i_length,LW_struct *i_length_CC,
				Data_lin *i_D_lga,Data_g_a *i_D_g_a,
				Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,
				int printSim)
{
  int            a,ncat_age,a2,f,h,i,ind_f,ind_a,cum_fish,aobs,season,type;
  int            l_int,l_int_prev,n_CC,n_Skrei;
  double         lobs,r;
  double         sum_p;
  double        *p,*p2,*mu,*sigma,*beta;
  long          *ages;
  double         u,prob,ptrue,pobs;
  
  FILE *fp;
  if(printSim)
    fp = fopen("ages_miss.txt","w");

 
  ncat_age = (int) i_D_age->glm->ncat/2;

  p = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  p2 = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  sigma = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  beta = CALLOC(i_D_lga->glm->nxcov,double);      // Free ok      

  for(a=0;a<ncat_age;a++)
    sigma[a] = G_ONE/sqrt(i_length_CC->par->tau_obs);
  for(a=0;a<ncat_age;a++)
    sigma[ncat_age+a] = G_ONE/sqrt(i_length->par->tau_obs);


  cum_fish = 0;
  ind_f = 0;

  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];
      /* Start by initializing sufficient statistics */ 
      /* All are simulated due to classification error */
      for(a=0;a<i_D_age->glm->ncat;a++)
	i_D_age->Ages[h][a] = 0;
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = 0;
	  i_D_lga->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	}
      for(a=0;a<i_D_g_a_CC->ncat;a++)
	{
	  i_D_lga_CC->Ages[h][a] = 0;
	  i_D_lga_CC->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga_CC->sqsum_by_cat[h][a] = G_ZERO;
	}
      /* find prior age-probabilities */
      sum_p = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] = exp(i_age->alpha[h][a]);
	  sum_p += p[a];
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] /= sum_p;
	}
      /* Find intercept and slope for lga model */
      //Coastal cod
      for(i=0;i<i_D_lga_CC->glm->nxcov;i++)
	beta[i] = calc_eff(i_D_lga_CC->glm->xcov[i],i_length_CC->par->eff[0][i],h);	  
      for(a=0;a<ncat_age;a++)
	{
	  mu[a] = beta[0] + beta[1]*i_D_g_a_CC->g_a[i_D_g_a_CC->a2Age_vec[a]+season-1];	
	}
      //Skrei
      for(i=0;i<i_D_lga->glm->nxcov;i++)
	beta[i] = calc_eff(i_D_lga->glm->xcov[i],i_length->par->eff[0][i],h);
      for(a=0;a<ncat_age;a++)
	{
	  mu[ncat_age+a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
	}

      l_int_prev = -1;
      l_int = 0;
      /* Loop through all fish */
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[ind_f];
	  type = i_D_orig->tottype[ind_f];
	  if(aobs > -1000 && lobs > -1000.0)/* Both age and length observed */
	    {
	      if(i_D_age->type_age[0]==0)/* No age error, simulate type only */
		{
		  n_CC = 0;
		  n_Skrei = 0;
		  if(type == 1) //coastal cod
		    {
		      pobs = i_D_CC->ptype1_CC[aobs];
		      ptrue = p[aobs];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      if(prob<0 || prob>1)
			write_warning("sample_ages_len_only_init_new_CC: Something is wrong");
		      for(i=0; i<i_D_orig->replength[ind_f]; i++)
			{
			  if(GMRFLib_uniform()<prob)
			    n_CC++;
			  else
			    n_Skrei++;
			}
		    }
		  else if(type == 2) //uncertain coastal cod
		    {
		      pobs = i_D_CC->ptype2_CC[aobs];
		      ptrue = p[aobs];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      for(i=0; i<i_D_orig->replength[ind_f]; i++)
			{
			  if(GMRFLib_uniform()<prob)
			    n_CC++;
			  else
			    n_Skrei++;
			}
		    }
		  else if(type == 4) //uncertain skrei
		    {
		      pobs = i_D_CC->ptype4_S[aobs];
		      ptrue = p[ncat_age+aobs];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      for(i=0; i<i_D_orig->replength[ind_f]; i++)
			{
			  if(GMRFLib_uniform()<prob)
			    n_Skrei++;
			  else
			    n_CC++;
			}
		    }
		  else if(type == 5) //skrei
		    {
		      pobs = i_D_CC->ptype5_S[aobs];
		      ptrue = p[ncat_age+aobs];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      for(i=0; i<i_D_orig->replength[ind_f]; i++)
			{
			  if(GMRFLib_uniform()<prob)
			    n_Skrei++;
			  else
			    n_CC++;
			}
		    }
		  // Add coastal cod
		  i_D_age->Ages[h][aobs] += n_CC;
		  ind_a = i_D_g_a_CC->a2Age_vec[aobs]+(season-1);
		  i_D_lga_CC->Ages[h][ind_a] += n_CC;
		  i_D_lga_CC->sum_by_cat[h][ind_a] += lobs * n_CC;
		  i_D_lga_CC->sqsum_by_cat[h][ind_a] += lobs*lobs * n_CC;
		  // Add skrei
		  i_D_age->Ages[h][ncat_age+aobs] += n_Skrei;
		  ind_a = i_D_g_a->a2Age_vec[aobs]+(season-1);
		  i_D_lga->Ages[h][ind_a] += n_Skrei;
		  i_D_lga->sum_by_cat[h][ind_a] += lobs * n_Skrei;
		  i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * n_Skrei;
		  if(printSim)
		    {
		      if(n_CC>0)
			fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				h+1,i_D_g_a_CC->a_vec[ind_a],lobs,n_CC,aobs,season,i_D_orig->tottype[ind_f],0);
		      if(n_Skrei>0)
			fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				h+1,i_D_g_a->a_vec[ind_a],lobs,n_Skrei,ncat_age+aobs,season,i_D_orig->tottype[ind_f],0);
		    }
		}
	      if(i_D_age->type_age[0]==1)/*Ages observed with error*/
		{
		  /* Both age and length observed */
		  /* Calculate probabilities */
		  sum_p = G_ZERO;
		  for(i=0;i<i_D_age->glm->ncat;i++)
		    p2[i] = G_ZERO;
		  for(i=0;i<i_age->A_Nneigh[aobs];i++) //calculate age probabilities for coastal cod
		    {
		      a = i_age->A_neigh[aobs][i];
		      if(type == 1) //coastal cod
			pobs = i_D_CC->ptype1_CC[a];
		      else if(type == 2) //uncertain coastal cod
			pobs = i_D_CC->ptype2_CC[a];
		      else if(type == 4) //uncertain skrei
			pobs = i_D_CC->ptype4_CC[a];
		      else if(type == 5) //skrei
			pobs = i_D_CC->ptype5_CC[a];
		      ptrue = p[a];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      p2[a] = p[a]*i_age->A2A[a][aobs]*dnorm(lobs,mu[a],sigma[a]);
		      p2[a] *= prob;
		      sum_p += p2[a];
		    }
		  for(i=0;i<i_age->A_Nneigh[aobs];i++) //calculate age probabilities for skrei
		    {
		      a = i_age->A_neigh[aobs][i];
		      if(type == 1) //coastal cod
			pobs = i_D_CC->ptype1_S[a];
		      else if(type == 2) //uncertain coastal cod
			pobs = i_D_CC->ptype2_S[a];
		      else if(type == 4) //uncertain skrei
			pobs = i_D_CC->ptype4_S[a];
		      else if(type == 5) //skrei
			pobs = i_D_CC->ptype5_S[a];
		      ptrue = p[ncat_age+a];
		      if(pobs==1)
			prob = 1;
		      else if(ptrue==0)
			prob = 0;
		      else
			prob = pobs*ptrue/(pobs*ptrue+(1-pobs)*(1-ptrue));
		      p2[ncat_age+a] = p[ncat_age+a]*i_age->A2A[a][aobs]*dnorm(lobs,mu[ncat_age+a],sigma[ncat_age+a]);
		      p2[ncat_age+a] *= prob;
		      sum_p += p2[ncat_age+a];
		    }
		  for(i=0;i<i_D_age->glm->ncat;i++)
		    p2[i] /= sum_p;
		  my_genmul(i_D_orig->replength[ind_f],p2,i_D_age->glm->ncat,ages);
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    {
		      if(ages[a] > 0)
			{
			  i_D_age->Ages[h][a] += (int) ages[a];
			  if(a<ncat_age) // Coastal cod
			    {		      
			      ind_a = i_D_g_a_CC->a2Age_vec[a]+(season-1);
			      i_D_lga_CC->Ages[h][ind_a] += (int) ages[a];
			      i_D_lga_CC->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			      i_D_lga_CC->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			      if(printSim)
				fprintf(fp,"%d %f %f %d %d %d %d %d\n",
					h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,1,1);
			    }
			  else // Skrei
			    {
			      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
			      i_D_lga->Ages[h][ind_a] += (int) ages[a];
			      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			      if(printSim)
				fprintf(fp,"%d %f %f %d %d %d %d %d\n",
					h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,5,1);
			    }
			}
		    }
		}
	    }
	  
	  if(aobs < -1000 && lobs > -1000.0)/* Age missing*/
	    {
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      if(l_int!=l_int_prev)
		{
		  //New calculation of age-probablities needs to be performed
		  /* Calculate probabilities */
		  sum_p = G_ZERO;
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    {
		      p2[a] = p[a]*dnorm(lobs,mu[a],sigma[a]);
		      sum_p += p2[a];
		    }
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    p2[a] /= sum_p;		  
		  l_int_prev = l_int;
		}
	      my_genmul(i_D_orig->replength[ind_f],p2,i_D_age->glm->ncat,ages);
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  if(ages[a]>0)
		    {
		      i_D_age->Ages[h][a] += (int) ages[a];
		      if(a<ncat_age) // Coastal cod
			{
			  ind_a = i_D_g_a_CC->a2Age_vec[a]+(season-1);
			  i_D_lga_CC->Ages[h][ind_a] += (int) ages[a];
			  i_D_lga_CC->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			  i_D_lga_CC->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			  if(printSim)
			    fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				    h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,1,1);
			}
		      else // Skrei
			{
			  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
			  i_D_lga->Ages[h][ind_a] += (int) ages[a];
			  i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
			  i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
			  if(printSim)
			    fprintf(fp,"%d %f %f %d %d %d %d %d\n",
				    h+1,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a],a,season,5,1);
			}
		    }
		}
	    }
	  ind_f++;
	}
      cum_fish += i_D_orig->nFishBoat[h];
    }
  if(printSim)
    fclose(fp);

  // Free memory allocated in this routine
  FREE(p);
  FREE(p2);
  FREE(mu);
  FREE(sigma);
  FREE(ages);
  FREE(beta);


  return(0);
}		/* end of sample_ages_len_only_new_CC */




/*!
  \author Geir Storvik
  \brief Calculate sufficient statistics for age model

  Sufficient statistics are for each haul
  - the number of fish per age
  - the sum of $g(a)$
  - the square sum of \f$g(a)\f$
  - Least squares intercept, slope and sum of squares
*/
int make_suff_age(int i_ncat,Age_struct *i_age,Data_age *i_D_age,double *i_haulweight,int i_start_h)
{
  int      a,h;

  if(i_D_age->glm->nxcov==1)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_ncat;a++)
	    {
	      i_D_age->glm->beta_hat[h][a][0] = i_age->alpha[h][a];
	      i_D_age->glm->suff[h][0][0] = G_ONE;
	    }
	}
    }
  else if(i_D_age->glm->nxcov==2)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_ncat;a++)
	    {
	      i_D_age->glm->beta_hat[h][a][0] = i_age->alpha[h][a];
              //In order to make update of b to work with old equations, 
              // the following is put to zero
	      i_D_age->glm->beta_hat[h][a][1] = G_ZERO;  
	    }
	  i_D_age->glm->suff[h][0][0] = G_ONE;
	  i_D_age->glm->suff[h][0][1] = i_haulweight[h];
	  i_D_age->glm->suff[h][1][0] = i_D_age->glm->suff[h][0][1];
	  i_D_age->glm->suff[h][1][1] = i_haulweight[h]*i_haulweight[h];
	}
    }
  else
    {
      write_warning("age->nxcov different from 1 or 2 is not implemented\n");
      return(1);
    }



  return(0);
}               /* end of make_suff_age */



/*!
  \author Geir Storvik
  \brief Calculate sufficient statistics for length given age model

  The sufficient statistics are conditional on the \f$g(a)\f$ function.

  Sufficient statistics are for each haul
  - the number of fish per age
  - the sum of $g(a)$
  - the square sum of \f$g(a)\f$
  - Least squares intercept, slope and sum of squares
*/
int make_suff_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,int i_start_h)
{
  int      a,h,maxAges;
  double   g_a,N_h,N_h_a,sum_g,sum_g2,sum_l,sum_l2,sum_gl,beta0,beta1,ssq;
  //  double  *sum_length, *sum_N_a;

  int err;
  for(h=i_start_h;h<i_D_lga->glm->nHaul;h++)
    {
      sum_g = G_ZERO;
      sum_g2 = G_ZERO;
      sum_gl = G_ZERO;
      sum_l = G_ZERO;
      sum_l2 = G_ZERO;
      maxAges = 0;
      N_h = G_ZERO;
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  g_a = i_D_g_a->g_a[a];
          N_h_a = (double) i_D_lga->Ages[h][a];
          sum_g +=  N_h_a * g_a;
          maxAges = max(maxAges,i_D_lga->Ages[h][a]);
          sum_g2 +=  N_h_a * g_a * g_a;
          sum_gl += g_a * i_D_lga->sum_by_cat[h][a];
	  sum_l += i_D_lga->sum_by_cat[h][a];
	  sum_l2 += i_D_lga->sqsum_by_cat[h][a];
          N_h += N_h_a;
	  //          sum_N_a[a] += N_h_a;
	  //          sum_length[a] += i_D_lga->sum_by_cat[h][a];
	}
      if(N_h < 0.0001)
	{
	  /* No data */
          beta0 = G_ZERO;
          beta1 = G_ZERO;
          ssq = G_ZERO;
	}
      else if(fabs(maxAges-N_h)<0.001)
	{ /* Only one age group sampled */
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g2-sum_g)<0.0001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g*sum_g-N_h*sum_g2)<0.0000000000001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else
	{
	  beta1 = (sum_l*sum_g-N_h*sum_gl)/(sum_g*sum_g-N_h*sum_g2);
	  beta0 = (sum_l-beta1*sum_g)/N_h;
	  ssq = sum_l2+beta0*beta0*N_h+beta1*beta1*sum_g2-
	    G_TWO*(beta0*sum_l+beta1*sum_gl-beta0*beta1*sum_g)+0.00000001;
	}
      if(ssq < G_ZERO || !(beta0 > -9999999.99 && beta0 < 999999999.99) ||
                         !(beta0 > -9999999.99 && beta0 < 999999999.99))
	{
	  printf("h=%d,beta0=%lf,beta1=%lf,ssq=%lf\n",h,beta0,beta1,ssq);
          for(a=0;a<i_D_g_a->ncat;a++)
	    printf("%d %f %f %f\n",
                    i_D_lga->Ages[h][a],g_a,
		    i_D_lga->sum_by_cat[h][a],i_D_lga->sqsum_by_cat[h][a]);
	  write_warning("make_suff_lga:Something is wrong\n");
	  return(1);
	}
      i_D_lga->glm->beta_hat[h][0][1] = beta1;
      i_D_lga->glm->beta_hat[h][0][0] = beta0;
      i_D_lga->glm->ssq[h] = ssq;
      i_D_lga->glm->suff[h][0][0] = N_h;
      i_D_lga->glm->suff[h][0][1] = sum_g;
      i_D_lga->glm->suff[h][1][0] = sum_g;
      i_D_lga->glm->suff[h][1][1] = sum_g2;

      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = N_h * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff[h][0][1]*i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][2] = N_h * i_D_lga->haulweight[h]*i_D_lga->haulweight[h];
	  
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff[h][0][2];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff[h][1][2];
	  i_D_lga->glm->beta_hat[h][0][2] = 0.0;
	}
    }


  return(0);
}               /* end of make_suff_lga */




/*!
  \author Geir Storvik
  \brief  Make graph for haul effects in age model, to be used in GMRFLib

  This routine is currently not used because other routines specialized for
  our problem is used instead
*/
int make_graph_age_ran(Age_struct *i_age,Data_age *i_D_age)
{
  int           i,err;
  Graph2_str   *gr;

  gr = i_age->gr_str_r;
  /* Find number of nodes */
  /* Allocate space */
  gr->acc = CALLOC(i_D_age->glm->nHaul,int);  // Free ok
  for(i=0;i<i_D_age->glm->nHaul;i++)
    gr->acc[i] = 0;
  gr->N = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  gr->mean = CALLOC(i_D_age->glm->ncat,double);   // Free ok
  gr->mu_old = CALLOC(i_D_age->glm->ncat,double); // Free ok
  gr->mu_new = CALLOC(i_D_age->glm->ncat,double); // Free ok
  gr->x_new = CALLOC(i_D_age->glm->ncat,double);  // Free ok
  gr->d = CALLOC(i_D_age->glm->ncat,double);      // Free ok

  err = GMRFLib_create_graph(&gr->graph);         // Free ok
  if(err)
    {
      write_warning("make_graph_age_ran:Error calling GMRFLib_create_graph\n");
      return(err);
    }
  gr->graph->n = i_D_age->glm->ncat;              
  gr->graph->nnbs = CALLOC(gr->graph->n,int);     // Free ok
  gr->graph->nbs = CALLOC(gr->graph->n, int *);   // Free ok

  /* Build neighbor structure to be independent */
  for(i=0;i<gr->graph->n;i++)
      gr->graph->nnbs[i] = 0;

  err = GMRFLib_prepare_graph(gr->graph);
  if(err)
    {
      write_warning("make_graph_age_ran:Error calling GMRFLib_prepare_graph\n");
      return(err);
    }

  for(i=0;i<gr->graph->n;i++)
      gr->d[i] = G_ONE;
  
  return(0);
}		/* end of make_graph_age_ran */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated in make_graph_age_ran

  Not in use now for the same reasons as make_graph_age_ran
*/
int re_make_graph_age_ran(Age_struct *i_age,Data_age *i_D_age)
{
  int           err;
  Graph2_str   *gr;

  gr = i_age->gr_str_r;
  /* Allocate space */
  FREE(gr->acc);
  FREE(gr->N);
  FREE(gr->mean);
  FREE(gr->mu_old);
  FREE(gr->mu_new);
  FREE(gr->x_new);
  FREE(gr->d);

  err = GMRFLib_free_graph(gr->graph); 
  if(err)
    {
      write_warning("re_make_graph_age_ran:Error calling GMRFLib_create_graph\n");
      return(err);
    }

  return(0);
}		/* end of re_make_graph_age_ran */


/*!
  \author Geir Storvik
  \brief  Samples random effects in age model using the GMRFLib library

  Include a haul effect as random effects
  Starts sampling for h=i_start_h

  This routine is currently not used because other routines specialized for
  our problem is used instead
*/
int sample_age_ran(Age_struct *i_age,Data_age *i_D_age,
		   int i_start_h)
{
  int  a,h,i,acc;
  double  sum,n,mu_prop,fac=1.5,sum_old,sum_new;
  double lacc;
  char *arg_loglik_old[2];
  char *arg_loglik_new[2];
  Graph2_str *gr;
  static int it=0;
 
  if(it==0)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	{
	  sum = G_ZERO;
	  n = G_ONE;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      sum += exp(i_age->alpha[h][a]);
	      n += (double) i_D_age->Ages[h][a];
	    }
	  i_age->mu[h] = gengam(sum,n);
	}
    }

  /*s_optimize_param->fp = fopen("opt.txt","w");*/
  gr = i_age->gr_str_r;

  arg_loglik_old[0] = (char *) gr->N;
  arg_loglik_old[1] = (char *) gr->mu_old;
  arg_loglik_new[0] = (char *) gr->N;
  arg_loglik_new[1] = (char *) gr->mu_new;
  acc = 0;
  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
    {
      fprintf(stderr,"%d",h);
      /* First sample i_age->mu using scale proposal */
      mu_prop = scale_proposal(i_age->mu[h],fac,NULL);
      /* Set parameters */
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
          gr->N[a] = (double) i_D_age->Ages[h][a];
	  gr->mu_old[a] = i_age->mu[h];
	  gr->mu_new[a] = mu_prop;
	  gr->mean[a] = G_ZERO;
	  for(i=0;i<i_D_age->glm->nxcov;i++)
	    gr->mean[a] += calc_eff(i_D_age->glm->xcov[i],i_age->par->eff[a][i],h)*i_D_age->glm->suff[h][0][i];
	}
      /* Then sample alpha */
      GMRFLib_blockupdate2(&lacc,gr->x_new,i_age->alpha[h],
			   NULL,NULL,  /* b=0 */
			   NULL,NULL,  /* c=0 */
			   gr->mean,gr->mean,  /* Means=0 */
			   gr->d,gr->d,  
			   Loglik_poisson_age,(char *)arg_loglik_new,
			   Loglik_poisson_age,(char *)arg_loglik_old,
			   NULL,gr->graph,
			   gr->Qfunc,NULL,  /* Not upd hyperparameters */
			   gr->Qfunc,NULL,
			   NULL,NULL,NULL,NULL, /* No Qfunc2_old2_new */
			   s_optimize_param,
			   s_hidden_par);
      sum_old = G_ZERO;
      sum_new = G_ZERO;
      n = G_ONE;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  sum_old += exp(i_age->alpha[h][a]);
	  sum_new += exp(gr->x_new[a]);
	  n += (double) i_D_age->Ages[h][a];
	}
      /*
      l_old = n * log(i_age->mu[h]) - i_age->mu[h] * sum_old;
      l_new = n * log(mu_prop) - mu_prop * sum_new;
      lacc += l_new - l_old;
      */
      if(it<2 || (*GMRFLib_uniform)()<exp(MIN(0,lacc)))
	{
	  /* Accepting new samples */
          i_age->mu[h] = mu_prop;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_age->alpha[h][a] = gr->x_new[a];
	  gr->acc[h]++;
	  /*          fprintf(stderr,"Accepted,lacc=%lf\n",lacc);*/
          acc++;
	}
      /*
      else
          fprintf(stderr,"Not accepted,lacc=%lf\n",lacc);
      if(it>1)
        fprintf(stderr,"h=%d,acc=%lf\n",h,(double) gr->acc[h]/(double) (it+1));
      */
      
    }
  fprintf(stderr,"it=%d,acc=%lf\n",it,(double) acc/
	  (double) (i_D_age->glm->nHaul*(it+1)));

  it++;

  return(0);
}		/* end of sample_age_ran */


/*!
  \author Geir Storvik
  \brief Samples random effects in age model.

  Include a haul effect as random effect. Starts sampling for h=i_start_h

  The main simulation is performed per haul through the 
  sample_age_alpha_ages_given routine.
*/
int sample_age_alpha(Age_struct *i_age,Data_age *i_D_age,
			int i_start_h,int i_acc,int i_it,int *acc_h)
{
  int     err,h,acc,acc_tot,acc_min,acc_max;


  acc_min = 1000000;
  acc_max = 0;
  acc_tot = G_ZERO; 
  if(i_it==0)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	acc_h[h] = 0;
    }
  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
    {
      // Alternative using Gibbs sampling on each alpha's separately
      // Used for testing the sample_age_alpha_ages_given routine
      //err = sample_age_alpha_ages_given_gibbs(h,i_age,i_D_age,i_acc,&acc);

      err = sample_age_alpha_ages_given(h,i_age,i_D_age,i_acc,&acc);
      if(err)
	{
	  write_warning("sample_age_alpha:Error calling sample_age_alpha_ages_given\n");
	  return(err);
	}

      acc_tot += acc;
      acc_h[h] += acc;
      acc_min = MIN(acc_min,acc_h[h]);
      acc_max = MAX(acc_max,acc_h[h]);
    }
  
  return(0);
}		/* end of sample_age_alpha */


/*!
  \author Geir Storvik
  \brief Routine for sampling alpha's inside a haul.

  This routine samples \f$\alpha_{h,a},a=1,...,A\f$ conditioned on 
  all fixed and random effects (except the haul effect), other hyper-parameters
  and observed ages. Note that simulating the \f$\alpha_{h,a}\f$'s is
  equivalent to sampling the haul effects.

  The simulation is performed through Metropolis-Hastings independence sampling
  steps by construting a proposal distribution through optimization of the
  conditional distribution. A Gaussian proposal with mean in the mode and
  covariance given from the hessian is used.

  Since the optimization requires some computer time,
  NGIBBS iterations of the Metropolis-Hastings independence sampling is performed.

  A problem with the independence sampling can be that the density for the
  backwards proposal is very small, making acceptance very unlikely. A possibility
  in that case is to scale the covariance matrix such that both jump forwards
  and backwards becomes more likely. This is implemented through the sc_old
  and sc_new variables. Currently these are set to 1, though.

  080401: Modified the routine to sample only the i_ncat-1 first alpha's under the
  constraint that all sum to zero. The age_haul_find_mode routine is changed accordingly.
*/
static int sample_age_alpha_ages_given(int i_h,Age_struct *i_age,
				      Data_age *i_D_age,int i_force_acc,int *o_acc)
{
  int     a,a2,it,err,noacc,i;
  double  N_h,sum,u,tau,ssq_old,ssq_new;
  double  log_q_old,log_q_new,log_d_old,log_d_new,log_d_opt,log_old,log_new;
  double  sc_old,sc_new;

  tau = i_age->par->tau_obs;
  /* Initializing */
  N_h = G_ZERO;
  sum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      N_h += (double) i_D_age->Ages[i_h][a];
      s_mu[a] = 0;
      for(i=0;i<i_D_age->glm->nxcov;i++)
	s_mu[a] += calc_eff(i_D_age->glm->xcov[i],i_age->par->eff[a][i],i_h)*i_D_age->glm->suff[i_h][0][i];
      s_mu[a] -= i_age->par->eff[a][0][i_D_age->glm->xcov[0]->n_cov-1][i_h]; //Not haul effect
    }

  /* Find mode */
  err = age_haul_find_mode(i_D_age->glm->ncat,i_D_age->Ages[i_h],N_h,s_mu,tau,
			   s_alpha_opt,s_prob,&log_d_opt,s_grad,s_Hess);
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      if(!(s_alpha_opt[a] > -99999999999.99 && s_alpha_opt[a] < 999999999.99))
	{
	  printf("h=%d,tau=%lf\n",i_h,tau);
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    printf("a=%d,Ages=%d, mu=%lf,s_alpha_opt=%lf\n",
		   a2,i_D_age->Ages[i_h][a2],s_mu[a2],s_alpha_opt[a2]);
	  write_warning("sample_age_alpha_ages_given:something is wrong\n");
	  return(1);
	}
    }

  /* Calculate ssq for old sample */ 
  for(a=0;a<(i_D_age->glm->ncat-1);a++)
    s_eps2[a] = i_age->alpha[i_h][a]-s_alpha_opt[a];
  ssq_old = cholssq0(s_Hess,i_D_age->glm->ncat-1,s_eps2);

  /* Calculate proposal log-likelihood for old sample */
  err = age_haul_calc_prob(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob);
  log_d_old = age_haul_calc_posterior(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob,
				    i_D_age->Ages[i_h],s_mu,tau);

  // Make a proper scaling so that jump backwards become more likely
  if(log_d_opt > log_d_old)
    sc_old = ssq_old/(G_TWO*(log_d_opt-log_d_old));
  else
    sc_old = G_ONE;
  // Turning of the scaling...
  sc_old = G_ONE;
  log_q_old = -G_HALF*ssq_old/sc_old;
  log_old = log_d_old - log_q_old;

  /* Sample */
  
  noacc = 1;
  it = 0;
  while(it < NGIBBS)
    {
      /* Draw from proposal */
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	s_eps[a] = gennor(G_ZERO,G_ONE/sqrt(sc_old));
      chollTl0(s_Hess,i_D_age->glm->ncat-1,s_eps,s_eps2);  //chollTl?
      #ifdef DEBUG_HAUL_OPT
      int a2;
      FILE  *unit;
      unit = fopen("Hchol.dat","w");
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	{
	  for(a2=0;a2<(i_D_age->glm->ncat-1);a2++)
	    fprintf(unit,"%lf ",s_Hess[a][a2]);
	  fprintf(unit,"\n");
	}
      fclose(unit);
      unit = fopen("eps_eps2.dat","w");
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	fprintf(unit,"%lf %lf\n",s_eps[a],s_eps2[a]);
      fclose(unit);
      #endif       
      
      sum = G_ZERO;
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	{
	  s_alpha_prop[a] = s_alpha_opt[a] + s_eps2[a];
          if(!(s_alpha_prop[a] > -99999999999.99 && s_alpha_prop[a] < 999999999.99))
	    {
              printf("h=%d\n",i_h);
              for(a2=0;a2<i_D_age->glm->ncat;a2++)
		printf("Ages[%d]=%d, s_alpha_opt[%d]=%lf\n",
		       a2,i_D_age->Ages[i_h][a2],a2,s_alpha_opt[a]);
	      write_warning("sample_age_alpha_ages_given:something is wrong\n");
	      return(1);
	    }
	  sum += s_alpha_prop[a];
	}
      s_alpha_prop[i_D_age->glm->ncat-1] = -sum;
      // Calculate age-probabilities for new sample
      err = age_haul_calc_prob(i_D_age->glm->ncat,s_alpha_prop,s_prob);
      log_d_new = age_haul_calc_posterior(i_D_age->glm->ncat,s_alpha_prop,s_prob,
					i_D_age->Ages[i_h],s_mu,tau);

      /* Calculate ssq for new sampl */
      ssq_new = G_ZERO;
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	ssq_new += s_eps[a]*s_eps[a];
      // Scaling for forward jump
      sc_new = ssq_new/(G_TWO*(log_d_opt-log_d_new));
      // Scaling turned off....
      sc_new = G_ONE;
      /* Calculate proposal likelihoods */
      log_q_new = -G_HALF*ssq_new/sc_new;
      log_new = log_d_new - log_q_new;

      // M_H acceptance
      u = genunf(G_ZERO,G_ONE);
      if(i_force_acc || u < exp(log_new-log_old))
	{
	  /* Accept sample */
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_age->alpha[i_h][a] = s_alpha_prop[a];
	  log_old = log_new;
	  noacc = 0;
	}
      it++;
    }
  err = age_haul_calc_prob(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob);
  *o_acc = 1-noacc;          
  return(0);
}		/* end of sample_age_alpha_ages_given */



/*!
  \author Geir Storvik
  \brief Finds optimal values of alpha for a given haul.

  Optimization is performed through Newton-Raphson steps with initial
  values based on observed ages and hyper-parameters (i.e not on old value
  which is important in order to make it possible to use in an independence sampler)

  080401: Optimizing on first A-1 alpha's, using that the last one is minus the sum 
  of the   other. Note that the covariance in the prior for the first A-1 alpha's 
  conditional on that the sum of all are zero is I+11^T.
  
  Returns both mode, gradient, hessian, the age-probabilities and 
  conditional log-posterior density
*/
static int age_haul_find_mode(int i_ncat,int *i_Ages,double i_N,
			      double *i_mu,double i_tau,
			      double *x_alpha_opt,double *o_prob,double *o_log_opt,
			      double *o_grad,double **o_Hess)
{
  int      it,more,err,a,a2;
  double   log_opt_old,log_opt,scale,mean_alpha,sum;

  // Initial values
  mean_alpha = G_ZERO;
  sum = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      x_alpha_opt[a] = G_HALF*(i_tau*i_mu[a]+log(0.01 + (double) i_Ages[a]));
      mean_alpha += x_alpha_opt[a];
      sum += i_mu[a];
    }
  //if(fabs(sum)>0.000001)
  //  write_warning("sum of mu not equal to zero\n");
  mean_alpha /= (double) i_ncat;
  for(a=0;a<i_ncat;a++)
    x_alpha_opt[a] -= mean_alpha;
  err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);
  log_opt_old = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);

  it=0;
  more = 1;
  // Run maximum 4 iterations with NR on first i_ncat-1 alpha's
  while(more && it < 4)
     {
       for(a=0;a<(i_ncat-1);a++)
	 {
	   o_grad[a] = -i_tau*(x_alpha_opt[a]-i_mu[a]-x_alpha_opt[i_ncat-1]+i_mu[i_ncat-1]) + 
	     (double) (i_Ages[a]-i_Ages[i_ncat-1])-i_N*(o_prob[a]-o_prob[i_ncat-1]);
	   o_Hess[a][a] = G_TWO*i_tau + i_N*o_prob[a]*(G_ONE-o_prob[a]) +
	                  i_N*o_prob[i_ncat-1]*(G_ONE-o_prob[i_ncat-1])+
                          G_TWO*i_N*o_prob[a]*o_prob[i_ncat-1];
	   for(a2=a+1;a2<(i_ncat-1);a2++)
	     {
	       o_Hess[a][a2] = i_tau + i_N*o_prob[i_ncat-1]-
                               i_N*(o_prob[a]-o_prob[i_ncat-1])*
                                   (o_prob[a2]-o_prob[i_ncat-1]);
	       o_Hess[a2][a] = o_Hess[a][a2];
	     }
	 }
#ifdef DEBUG_HAUL_OPT
       FILE *unit;
       unit = fopen("Hess.dat","w");
       for(a=0;a<(i_ncat-1);a++)
	 {
	   for(a2=0;a2<(i_ncat-1);a2++)
	     fprintf(unit,"%lf ",o_Hess[a][a2]);
	   fprintf(unit,"\n");
	 }
       fclose(unit);
#endif
       err = choldc0(o_Hess,i_ncat-1);
       if(err)
	 {
	   write_warning("sample_age_find_mode:Error calling choldc0\n");
	   return(err);
	 }
       cholsl0_old(o_Hess,i_ncat-1,o_grad,s_delta);
#ifdef DEBUG_HAUL_OPT
       unit = fopen("grad_detaold_delta.dat","w");
       for(a=0;a<(i_ncat-1);a++)
	 fprintf(unit,"%lf %lf\n",o_grad[a],s_delta[a]);
       fclose(unit);
#endif       

       sum = G_ZERO;
       for(a=0;a<(i_ncat-1);a++)
	 {
	   x_alpha_opt[a] = x_alpha_opt[a]+s_delta[a];
	   sum += x_alpha_opt[a];
	 }
       x_alpha_opt[i_ncat-1] = -sum;

       err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);

       log_opt = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);
       scale = G_ONE;
       // If new value do not give a better value, use a smaller Newton-Raphson step
       while(log_opt < (log_opt_old-0.00000001))
	 {
	   scale /= G_TWO;
	   sum = G_ZERO;
	   for(a=0;a<(i_ncat-1);a++)
	     {
	       x_alpha_opt[a] = x_alpha_opt[a]-scale*s_delta[a];
	       sum += x_alpha_opt[a];
	     }
	   x_alpha_opt[i_ncat-1] = -sum;

	   err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);

	   log_opt = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);
	 }
       more = fabs(log_opt_old-log_opt) > 0.001;
       log_opt_old = log_opt;
       it++;
    }

  *o_log_opt = log_opt;

  return(0);
}		/* end of age_haul_find_mode */


/*!
  \author Geir Storvik
  \brief Calculates posterior for alpha given observed ages and prior expectations
*/
static double age_haul_calc_posterior(int i_ncat,double *i_alpha,double *i_prob,
				      int *i_Ages,double *i_mu,double i_tau)
{
  int       a;
  double    log_pr,res,sum,log_d;

  /* Prior */
  log_pr = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      res = i_alpha[a] - i_mu[a];
      log_pr -= res*res;
    }
  log_pr *= G_HALF*i_tau;
  /* Likelihood */
  sum = G_ZERO;
  log_d = G_ZERO;
  for(a=0;a<i_ncat;a++)
    log_d += log(i_prob[a]) * (double) i_Ages[a];

  return(log_pr+log_d);
}		/* end of age_haul_calc_posterior */


/*!
  \author Geir Storvik
  \brief Calculates probabilities for age-groups given alpha's
*/
static int age_haul_calc_prob(int i_ncat,double *i_alpha,double *o_prob)
{
  int    a;
  double sum;
  
  sum = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      o_prob[a] = exp(i_alpha[a]);
      sum += o_prob[a];
    }
  for(a=0;a<i_ncat;a++)
      o_prob[a] /= sum;
  
  return(0);
}


/*!
  \author Geir Storvik
  \brief Calculates loglikelihood of poisson distribution for age-data

  This routine is used in sample_age_ran where the multinomial distribution 
  is extended by a intensity variable in order to make the age-data
  poisson distributed.
*/
static int Loglik_poisson_age(double *logll,double *x,
		       int m,int idx,double *x_vec,char *arg)
{
  int i;
  double *N, *mu;
  char **args;

  args = (char **)arg;
  N = (double *)args[0];
  mu = (double *)args[1];

  for(i=0;i<m;i++) logll[i] = N[idx]*x[i]-mu[idx]*exp(x[i]);

  return GMRFLib_SUCCESS;
}		/* end of Loglik_poisson_age */



/*!
  \author Geir Storvik
  \brief Samples precision parameters for haul effect in age model
  080401:  Modified in order to take into account that sum_a alpha_{h,a}=0
*/
int sample_precision_age_haul(int i_start_h,Eff_str *i_par,
                              double **i_alpha,Data_glm *i_glm)
{
  int    n,a,h,i;
  double mu,ssq;
  Data_cov *xcov;

  n = (i_glm->nHaul-i_start_h)*(i_glm->ncat-1);
  ssq = 0;
  for(h=i_start_h;h<i_glm->nHaul;h++)
    {
      xcov = i_glm->xcov[0];
      for(a=0;a<i_glm->ncat;a++)
	{
	  mu = 0;
	  for(i=0;i<i_glm->nxcov;i++)
	    mu += calc_eff(i_glm->xcov[i],i_par->eff[a][i],h)*i_glm->suff[h][0][i];
          mu -= i_par->eff[a][0][i_glm->xcov[0]->n_cov-1][h]; //Not haul effect
	  ssq += pow(i_alpha[h][a] - mu,2);
	  //xcov->n_cov++;  
	}
      if(!(ssq >= 0 && ssq < 9999999999999.99))
	{
	  fprintf(stderr,"ssq=%lf\n",ssq);
	  write_warning("sample_precision_age_haul:Something is wrong\n");
	  return(1);
	}
    }
  if(!(ssq > 0 && ssq < 9999999999999.99))
    {
      fprintf(stderr,"ssq=%lg\n",ssq);
      write_warning("sample_precision_age_haul:Something is wrong\n");
      return(1);
    }
  i_par->tau_obs = gengam(i_par->prior_prec_obs[0]+G_HALF * ssq,
			  i_par->prior_prec_obs[1]+G_HALF * (double) n);

  return(0);
}		/* end of sample_precision_age_haul */



/*!
  \author Hanne Rognebakke
  \brief Initial sampling of lengths for discarded fish

  To be used in COST project
*/
int sample_discard_init2_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,
			      Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_orig *i_D_orig,
			      Data_COST *i_D_COST,int printDisc)
{
  int     a,i,h,f,ind,ind_a,ind_l;
  int f_start,f_end,N_ml_land;
  int     N_int,ncat,season,Ndisc;
  double  age,ldisc,Int,Slp,log_a0,log_adiff,cens_k,cens_m,cens_r;
  double  sig_fish,ldisc_mid,ldisc_end,ldisc_start,p_disc,sum_a,sum_p,pnorm_old,tmp;
  int    *Ages;
  double *p_l,*mu;
  long   *lengths;
  FILE   *fp;

  if(printDisc)
    fp = fopen("discard_init.txt","w");

  /* Initialize */
  N_int = i_D_COST->mland->N_int_disc;
  ncat = i_D_age->glm->ncat;

  Ages = CALLOC(ncat,int); // Free ok
  p_l = CALLOC(N_int,double); // Free ok
  mu = CALLOC(ncat,double); // Free ok
  lengths = CALLOC(N_int,long); // Free ok

  ind = 0;
  for(h=0;h<i_D_COST->mland->n_trip;h++)
    for(f=0;f<N_int;f++)
      {
	i_D_COST->mland->lfreq_disc[ind] = 0;
	ind++;
      }
  for(a=0;a<ncat;a++)
    Ages[a] = 0;
  for(h=i_D_COST->obs->n_trip;h<i_D_age->glm->nHaul;h++)
    for(a=0;a<ncat;a++)
      Ages[a] += i_D_age->Ages[h][a];
  //Increase probability of being in smallest age groups due to no observed discards
  for(a=0;a<ncat;a++)
    if(Ages[a]==0 && Ages[a+1]>0)
      {
	Ages[a+1] *= 2;
	Ages[a] = Ages[a+1];
	break;
      }
  sum_a = 0;
  for(a=0;a<ncat;a++)
    sum_a += Ages[a];
  
  Int = i_length->par->eff[0][0][0][0];
  Slp = i_length->par->eff[0][1][0][0];
  log_a0 = log(i_D_g_a->a_vec[0]);
  log_adiff = log(i_D_g_a->a_vec[i_D_g_a->ncat-1])-log(i_D_g_a->a_vec[0]);
  cens_k = i_D_COST->cens->k;
  cens_m = i_D_COST->cens->m;
  sig_fish = sqrt(G_ONE/i_length->par->tau_obs);
  sig_fish = 0.5;

  /* Simulate discards*/
  ind_l = 0;
  for(h=i_D_COST->obs->n_trip;h<i_D_age->glm->nHaul;h++)
    {      
      season = i_D_orig->season[h];
      cens_r = i_D_COST->cens->r[h];
 
      //Start simulating number of discards equal to number of landed
      Ndisc = i_D_orig->n_landed[h];
      
      sum_p = G_ZERO;
      pnorm_old = G_ZERO;
      p_l[0] = G_ZERO;
      for(a=0;a<ncat;a++)
	mu[a] = Int + Slp*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
      for(i=1;i<N_int;i++)
	{
	  if(i_D_COST->mland->l_disc[ind_l+i]>1.6)
	    {
	  ldisc_mid = i_D_COST->mland->l_disc[ind_l+i];
	  ldisc_end = i_D_COST->mland->int_len_lim[i];
	  ldisc_start = i_D_COST->mland->int_len_lim[i-1];
	  p_disc = cens_function(ldisc_start,ldisc_mid,ldisc_end,cens_k,cens_m,cens_r);
	  p_l[i] = G_ZERO;
	  for(a=0;a<ncat;a++)
	    {
	      tmp = pnorm(ldisc_end,mu[a],sig_fish);
	      tmp = dnorm(ldisc_mid,mu[a],sig_fish);
	      p_l[i] += tmp*Ages[a]/sum_a;
	      pnorm_old = tmp;
	    }
	  p_l[i] *= p_disc;
	  sum_p += p_l[i];
	    }
	  else
	    p_l[i] = G_ZERO;
	}
      for(i=1;i<N_int;i++)
	p_l[i] /= sum_p;

      my_genmul(Ndisc,p_l,N_int,lengths);
      for(i=0;i<N_int;i++)
	{
	  if(lengths[i]>0)
	    {
	      ldisc = i_D_COST->mland->l_disc[ind_l+i];
	      age = exp(log_a0 + log_adiff*(ldisc-Int)/Slp);
	      if(age < i_D_age->a_vec[0])
		age = i_D_age->a_vec[0];
	      if(age > i_D_age->a_vec[i_D_age->glm->ncat-1])
		age = i_D_age->a_vec[i_D_age->glm->ncat-1];		
	      a = (int) age;
	      i_D_COST->mland->lfreq_disc[ind_l+i] += (int) lengths[i];
	      i_D_age->Ages[h][a] += (int) lengths[i];
	      i_D_age->Ages_disc[h][a] += (int) lengths[i];
	      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
	      i_D_lga->Ages[h][ind_a] += (int) lengths[i];
	      i_D_lga->sum_by_cat[h][ind_a] += ldisc * (double) lengths[i];
	      i_D_lga->sqsum_by_cat[h][ind_a] += ldisc*ldisc * (double) lengths[i];
	      if(printDisc)
		fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],ldisc,(int) lengths[i]);
	    }
	}
      ind_l += i_D_COST->mland->N_int_disc;
    } // end for(h=i_D_COST->obs->n_trip;h<i_D_age->glm->nHaul;h++)
  
  FREE(Ages);
  FREE(p_l);
  FREE(mu);
  FREE(lengths);
  
  return(0);
}		/* end of sample_discard_init2_COST */


/*!
  \author Hanne Rognebakke
  \brief Initial sampling of lengths for discarded fish

  To be used in COST project
*/
int sample_discard_init_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,
			     Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_orig *i_D_orig,
			     Data_COST *i_D_COST,int printDisc)
{
  int    a,i,h,s,f,f_start,f_end,l_int,l_int_start,ind,ind_l,ind_f,nMiss,miss;
  int    N_int,Ndisc,ncat,nSeason,season,ind_a;
  int    N_obs,N_obs_disc,N_obs_land,N_ml_land;
  double sum,cens_m,cens_r,lobs,fac,ldisc;
  double *P_al;
  int    **N_al;
  int    *N_l;
  long   *al_sim;
  FILE   *fp;
  if(printDisc)
    fp = fopen("discard_init.txt","w");

  /* Initialize */
  N_int = i_D_orig->n_int_len;
  ncat = i_D_age->glm->ncat;
  nSeason = i_D_g_a->nSeason;
  ind = 0;
  for(h=0;h<i_D_COST->mland->n_trip;h++)
    for(f=0;f<i_D_COST->mland->N_int_disc;f++)
      {
	i_D_COST->mland->lfreq_disc[ind] = 0;
	ind++;
      }

  /* Probability matrix P(l) for observer data */
  N_l = CALLOC(N_int,int); //Free ok
  for(i=0;i<N_int;i++)
    N_l[i] = 0;

  /* Frequency matrix P(a,l) for discards using observer data */
  P_al = CALLOC(ncat*N_int,double); //Free ok
  N_al = Mmatrix_2d(0,nSeason,0,ncat*N_int,sizeof(int),1); //Free ok
  al_sim = CALLOC(ncat*N_int,long); //Free ok
  for(s=0;s<nSeason;s++)
    for(i=0;i<ncat*N_int;i++)
      N_al[s][i] = 0;
  
  N_obs_disc = 0;
  N_obs_land = 0;
  for(h=0;h<i_D_COST->obs->n_trip;h++)
    {
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      season = i_D_orig->season[h];
      for(f=f_start;f<f_end;f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[f];
	  if(a > -1000 && lobs > -1000.0)
	    {
	      l_int = 0;
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      N_al[season-1][a*N_int+l_int] += i_D_orig->discard[f];
	    }
	  if(lobs > -1000.0)
	    {
	      N_l[l_int] += i_D_orig->discard[f];
	      N_l[l_int] += i_D_orig->landed[f];
	    }
	  N_obs_disc += i_D_orig->discard[f];
	  N_obs_land += i_D_orig->landed[f];
	}
    }
  N_obs = N_obs_disc+N_obs_land;

  //Copy for missing seasons
  miss = 1;
  nMiss = 0;
  for(s=0;s<nSeason;s++)
    {
      sum = G_ZERO;
      for(i=0;i<N_int*ncat;i++)
	sum += (double) N_al[s][i];
      if(sum<0.0001)
	{
	  if(s==0)
	    {
	      nMiss++;
	      miss = 1;
	    }
	  else if(s>0 & miss==1)
	    nMiss++;
	  else
	    for(i=0;i<N_int*ncat;i++)
	      N_al[s][i] = N_al[s-1][i];
	}
      else
	miss = 0;
    }
  for(s=0;s<nMiss;s++)
    for(i=0;i<N_int*ncat;i++)
      N_al[s][i] = N_al[nMiss][i];

  /* Simulate discards*/
  ind_l = 0;
  for(h=i_D_COST->obs->n_trip;h<i_D_age->glm->nHaul;h++)
    {      
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      N_ml_land = 0;
      for(f=f_start;f<f_end;f++)
	N_ml_land += i_D_orig->landed[f];

      season = i_D_COST->mland->season[h-i_D_COST->obs->n_trip];
      cens_m = i_D_COST->cens->m;
      cens_r = i_D_COST->cens->r[h];
      sum = G_ZERO;
      l_int = 0;
      while(cens_r*(1-1/sqrt(cens_m)) > i_D_orig->int_len_lim[l_int])
	{
	  sum += (double) N_l[l_int];
	  l_int++;
	}	  
      l_int_start = l_int;
      while(cens_r > i_D_orig->int_len_lim[l_int])
	{
	  sum += (double) N_l[l_int];
	  l_int++;
	}
      if(cens_r < i_D_COST->cens->mu[2])
	Ndisc = (int)(N_ml_land*(double)sum/(N_obs-sum));
      else
	Ndisc = (int)(N_ml_land);
      if(Ndisc==0 && N_ml_land>0)
	Ndisc = 1;
      
      // Convert to probabilities
      ind = 0;
      for(a=0;a<ncat;a++)
	for(i=0;i<N_int;i++)
	  {
	    if(i>l_int_start && i<l_int)
	      P_al[ind] = (double) N_al[season-1][ind];
	    else
	      P_al[ind] = G_ZERO;
	    ind++;
	  }
      sum = G_ZERO;
      for(i=0;i<ncat*N_int;i++)
	sum += P_al[i];
      if(sum<0.0001)
	{
	  write_warning("sample_discard_init_COST:Something is wrong\n");
	}
      else
	{
	  for(i=0;i<ncat*N_int;i++)
	    P_al[i] /= sum;
	}
      
      my_genmul(Ndisc,P_al,N_int*ncat,al_sim);
      
      ind=0;
      for(a=0;a<ncat;a++)
	{
	  for(i=0;i<N_int;i++)
	    {
	      if(al_sim[ind]>0)
		{
		  ldisc = i_D_orig->int_len[i];
		  l_int = 0;
		  while(ldisc > i_D_COST->mland->l_disc[l_int])
		    l_int++;
		  i_D_COST->mland->lfreq_disc[ind_l+l_int] += (int) al_sim[ind];
		  i_D_age->Ages[h][a] += (int) al_sim[ind];
		  i_D_age->Ages_disc[h][a] += (int) al_sim[ind];
		  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		  i_D_lga->Ages[h][ind_a] += (int) al_sim[ind];
		  i_D_lga->sum_by_cat[h][ind_a] += ldisc * (double) al_sim[ind];
		  i_D_lga->sqsum_by_cat[h][ind_a] += ldisc*ldisc * (double) al_sim[ind];
		  if(printDisc)
		    fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],
			    i_D_COST->mland->l_disc[ind_l+l_int],i_D_COST->mland->lfreq_disc[ind_l+l_int]);
		}
	      ind++;
	    }
	}
      ind_l += i_D_COST->mland->N_int_disc;
    }
   
  if(printDisc)
    fclose(fp);

  //printf("Simulate discard finished\n");
  Fmatrix_2d(&N_al[0][0],&N_al[0]);
  FREE(P_al);
  FREE(al_sim);
  FREE(N_l);

  return(0);
}		/* end of sample_discard_init_COST */



/*!
  \author Hanne Rognebakke
  \brief Samples lengths of discarded fish

  To be used in COST project
*/
int sample_discard_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,Data_lin *i_D_lga,
			Data_g_a *i_D_g_a,Data_orig *i_D_orig,Data_COST *i_D_COST,int printDisc, int i_it)
{
  int     i,a,f,h,N,nacc,sim_Ndisc,nit,season,ind_a,ind_l,lobs,ncat,N_int,err;
  long   *Ndisc;
  double  sig_fish,ldisc,len,prob,S,Sconst,p,z,sum_p,p_land,p_disc,arg;
  double  cens_k,cens_m,cens_r,c,d;
  double *mu,*beta,*p_a,*p_land_a,*p_a_disc,*p_l;
  long   *lengths;
  char    string[150];
  double  tmp,pnorm_old,ldisc_mid,p_disc_mid,ldisc_end,ldisc_start;
  int     Ndisc_tot,Nobs_tot,Nobs_land,f_start,f_end;

  double lambda;

  FILE   *fp,*fp2;
  if(printDisc)
    fp = fopen("discard.txt","w");
  

  ncat = i_D_age->glm->ncat;
  N_int = i_D_COST->mland->N_int_disc;

  p_a = CALLOC(ncat,double);      // Free ok
  p_land_a = CALLOC(ncat,double);      // Free ok
  p_a_disc = CALLOC(ncat,double);      // Free ok
  Ndisc = CALLOC(ncat,long);  // Free ok
  beta = CALLOC(i_D_lga->glm->nxcov,double); // Free ok
  mu = CALLOC(ncat,double); // Free ok
  p_l = CALLOC(N_int,double); // Free ok
  lengths = CALLOC(N_int,long); // Free ok

  Ndisc_tot = 0;
  Nobs_tot = 0;

  sig_fish = sqrt(G_ONE/i_length->par->tau_obs);


  /* Initialize */
  ind_l = 0;
  for(h=0;h<i_D_COST->mland->n_trip;h++)
    for(f=0;f<N_int;f++)
      {
	i_D_COST->mland->lfreq_disc[ind_l] = 0;
	ind_l++;
      }

  err = sample_lambda_prior_COST(i_D_COST);

  ind_l = 0;
  for(h=i_D_COST->obs->n_trip;h<i_D_age->glm->nHaul;h++)
    {
      /* find prior age-probabilities */
      sum_p = G_ZERO;
      for(a=0;a<ncat;a++)
	{
	  p_a[a] = exp(i_age->alpha[h][a]);
	  sum_p += p_a[a];
	}
      for(a=0;a<ncat;a++)
	{
	  p_a[a] /= sum_p;
	}
      
      cens_k = i_D_COST->cens->k;
      cens_m = i_D_COST->cens->m;
      cens_r = i_D_COST->cens->r[h];
      Sconst = cens_k/sqrt(G_PI);
      season = i_D_orig->season[h];

      for(i=0;i<i_D_lga->glm->nxcov;i++) // calculate b0 and b1
	beta[i] = calc_eff(i_D_lga->glm->xcov[i],i_length->par->eff[0][i],h);
      
      p_land = G_ZERO;
      for(a=0;a<ncat;a++)
	{
	  mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
	  p_land_a[a] = G_ZERO;
	  for(i=1;i<=s_N_gauher;i++)
	    {
	      len = s_gauher_x[i];
	      // using Gaussian cdf as S-function for probability of landing
	      arg = cens_m*(mu[a]+sqrt(2)*sig_fish*len-cens_r);
	      if(arg<-100)
		S = G_ZERO;
	      else if(arg>100)
		S = G_ONE;
	      else
		S = pnorm(arg,G_ZERO,G_ONE);
	      p_land_a[a] += s_gauher_w[i] * Sconst*S;
	    }
	  if(p_land_a[a]>1.000001)
	    {
	      sprintf(string,
		      "sample_discarded_COST:Something is wrong, p(landed|a=%d)=%10.8lf\n",
		      a,p_land_a[a]);
	      write_warning(string);
	      return(1);
	    }
	  p_land += p_a[a]*p_land_a[a];

	}
      if(p_land>1.000001)
	{
	  sprintf(string,"sample_discarded_COST:Something is wrong, p(landed)=%10.8lf\n",
		  p_land);
	  write_warning(string);
	  return(1);
	}

      Nobs_land = 0;
      for(a=0;a<ncat;a++)
	Nobs_land += i_D_age->Ages_land[h][a];
      Nobs_tot += Nobs_land;

      c = i_D_COST->mland->c;
      d = i_D_COST->mland->d;
      i_D_COST->mland->lambda[h-i_D_COST->obs->n_trip] = gengam(d+p_land,c+Nobs_land);
      i_D_COST->mland->lambda[h-i_D_COST->obs->n_trip] = 2*Nobs_land;

      i_D_orig->n_discard[h] = 0;
      for(a=0;a<ncat;a++)
	{
	  lambda = i_D_COST->mland->lambda[h-i_D_COST->obs->n_trip]*p_a[a]*(1-p_land_a[a]);
	  if(lambda < 0.0000001)
	    {
	      //printf("lambda=%f,p_land_a=%f\n",lambda,p_land_a[a]);
	      lambda = 0.0;
	      Ndisc[a] = 0;
	    }
	  else
	    Ndisc[a] = ignpoi(lambda);  
    
	  //if(Ndisc[a] > i_length->cens_Nlim*i_D_age->Ages_land[h][a] && i_D_age->Ages_land[h][a]>0)
	  if(0)
	    // if number of discarded to large, then adjust to limit
	    {
	      printf("h=%d, a=%d: Nland=%d, Ndisc=%d, ",h,a,i_D_age->Ages_land[h][a],Ndisc[a]);
	      Ndisc[a] = i_length->cens_Nlim*i_D_age->Ages_land[h][a];
	      printf("Ndisc_new=%d, Nlim=%f\n p(a)=%f, p(land|a)=%f, lambda=%f, prod=%f,r=%f\n",
		     Ndisc[a],i_length->cens_Nlim,p_a[a],p_land_a[a],i_D_COST->mland->lambda[h-i_D_COST->obs->n_trip],lambda,cens_r);
	    }

	  if(Ndisc[a]>0)
	    {
	      // use all the discards for estimating the age parameters
	      i_D_age->Ages[h][a] += Ndisc[a];
	      i_D_age->Ages_disc[h][a] += Ndisc[a];
	      i_D_orig->n_discard[h] += Ndisc[a];
	      Ndisc_tot += Ndisc[a];
	      // use only discards in age group where landed fish for estimating lga parameters
	      if(i_D_age->Ages_land[h][a]>0)
		{
		  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		  i_D_lga->Ages[h][ind_a] += Ndisc[a];
		}
	      // sample Ndisc lengths
	      sum_p = G_ZERO;
	      pnorm_old = G_ZERO;
	      p_l[0] = G_ZERO;
	      for(i=1;i<N_int;i++)
		{
		  ldisc_mid = i_D_COST->mland->l_disc[ind_l+i];
		  ldisc_end = i_D_COST->mland->int_len_lim[i];
		  ldisc_start = i_D_COST->mland->int_len_lim[i-1];
		  p_disc = cens_function(ldisc_start,ldisc_mid,ldisc_end,cens_k,cens_m,cens_r);
		  tmp = pnorm(ldisc_end,mu[a],sig_fish);
		  p_l[i] = (tmp-pnorm_old)*p_disc/(1-p_land_a[a]);
		  pnorm_old = tmp;
		  sum_p += p_l[i];
		}	      
	      for(i=0;i<N_int;i++)
		p_l[i] /= sum_p;
	      my_genmul(Ndisc[a],p_l,N_int,lengths);	  
	      for(i=0;i<N_int;i++)
		{
		  if(lengths[i]>0)
		    {
		      ldisc = i_D_COST->mland->l_disc[ind_l+i];
		      i_D_COST->mland->lfreq_disc[ind_l+i] += (int) lengths[i];
		      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		      
		      // use only discards in age group where landed fish for estimating lga parameters
		      if(i_D_age->Ages_land[h][a]>0)
			{
			  i_D_lga->sum_by_cat[h][ind_a] += ldisc*(double)lengths[i];
			  i_D_lga->sqsum_by_cat[h][ind_a] += ldisc*ldisc*(double)lengths[i];
			}
		      if(printDisc)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],ldisc,(int)lengths[i]);
		    }
		}
	    }
	}
      ind_l += N_int;
    }

  if(printDisc)
    fclose(fp); 

  FREE(Ndisc);
  FREE(beta);
  FREE(mu)
  FREE(p_a);
  FREE(p_land_a);
  FREE(p_a_disc);
  FREE(lengths);
  FREE(p_l);

  return(0);
}               /* end of sample_discard_COST */

double cens_function(double lstart,double lmid,double lend,double k,double m,double r)
{
  double arg,lint,p_disc;
  int i;
  int nint = 1;

  lint = (lend-lstart)/nint;

  if(1)
    {
      p_disc = G_ZERO;
      for(i=0;i<=nint;i++)
	{        
	  arg = m*(lstart+i*lint-r);
	  if(arg < -100)
	    p_disc += 1;
	  else if(arg > 100)
	    p_disc += 1-k;
	  else
	    p_disc += 1-k*pnorm(arg,G_ZERO,G_ONE);
	}
      p_disc = p_disc/(nint+1);
    }

  if(0)//using only midpoint
    {
      arg = m*(lmid-r);
      if(arg<-100)
	p_disc += 1;
      else if(arg>100)
	p_disc += 1-k;
      else
	p_disc += 1-k*pnorm(arg,G_ZERO,G_ONE);
    }

  return(p_disc);
}


/*!
  \author Hanne Rognebakke
  \brief Samples censoring parameters k, m, r
 */
int sample_cens_par(Data_lin *i_D_lga, Data_orig *i_D_orig, Data_COST *i_D_COST, int i_it)
{
  int i,h,p,t,ind_l,f,f_start,f_end;
  double fac,b,d,u,tau,x,arg;
  double cens_k,cens_m,cens_r,k_old,k_new,cens_func_old,cens_func_new;
  double log_prior_diff,log_disc_new,log_disc_old,log_land_new,log_land_old;
  double par_new,par_old;
  double sum;
  double lstart,lmid,lend;

  int N_disc,N_land;

  int printFile = 0;
  FILE *fp;

  if(printFile)
    {
      fp = fopen("cens_test.txt","w");
      fprintf(fp,"\n new it\n");
    }

  fac = 1.01;

  /* Sample hyperparameters, expected value for r */
  sum = G_ZERO;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      sum += i_D_COST->cens->r[h];
    }
  b = i_D_COST->cens->mu_prior_mean*i_D_COST->cens->mu_prior_prec+i_D_COST->cens->tau[2]*sum;
  tau = i_D_COST->cens->mu_prior_prec+(double)i_D_lga->glm->nHaul*i_D_COST->cens->tau[2];
  u = GMRFLib_stdnormal();
  i_D_COST->cens->mu[2] = b/tau + u/sqrt(tau);

  /* Sample hyperparameters, precision for r */
  sum = G_ZERO;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      sum += (i_D_COST->cens->r[h]-i_D_COST->cens->mu[2])*
	(i_D_COST->cens->r[h]-i_D_COST->cens->mu[2]);
    }
  i_D_COST->cens->tau[2] = gengam(i_D_COST->cens->b_prior+G_HALF*sum,
				  i_D_COST->cens->a_prior+G_HALF*(double)i_D_lga->glm->nHaul);

  if(printFile)
    {
      fprintf(fp,"Hyperparameters, expected value for r\n");
      fprintf(fp,"r: mu=%f,b=%f,tau=%f\n",i_D_COST->cens->mu[2],b,tau);
      fprintf(fp,"Hyperparameters, precision for r\n");
      fprintf(fp,"r: tau=%f,a=%f,b=%f\n\n",i_D_COST->cens->tau[2],
	      i_D_COST->cens->a_prior+G_HALF*(double)i_D_lga->glm->nHaul,i_D_COST->cens->b_prior+G_HALF*sum);
    }

  /* Sample parameters in censoring function */

  /* Sample k, same for all hauls */
  k_old = i_D_COST->cens->k;
  if(k_old = 1.0)
    par_old = 100;
  else
    par_old = log(k_old)-log(1-k_old);
  cens_m = i_D_COST->cens->m;
  par_new = scale_proposal(par_old,fac,NULL);
  k_new = exp(par_new)/(G_ONE+exp(par_new));

  /* Find probabilities for old and current sample */
  log_disc_new = G_ZERO;
  log_disc_old = G_ZERO;
  ind_l = 0;
  N_disc = 0;
  N_land = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      cens_r = i_D_COST->cens->r[h];
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      for(f=f_start;f<f_end;f++)
	{  
	  N_land += i_D_orig->landed[f];
	  if(i_D_orig->discard[f]>0)
	    {
	      lstart = log(round(exp(i_D_orig->totlength[f]))-0.5);
	      lmid = i_D_orig->totlength[f];
	      lend = log(round(exp(i_D_orig->totlength[f]))+0.5);
	      cens_func_old = cens_function(lstart,lmid,lend,k_old,cens_m,cens_r);
	      cens_func_new = cens_function(lstart,lmid,lend,k_new,cens_m,cens_r);
	      if(cens_func_old > G_ZERO)
		log_disc_old += i_D_orig->discard[f]*log(cens_func_old);
	      if(cens_func_new > G_ZERO)
		log_disc_new += i_D_orig->discard[f]*log(cens_func_new);
	      N_disc += i_D_orig->discard[f];
	    }
	}
      // add simulated discards for market landing data
      if(h>=i_D_COST->obs->n_trip)
	{
	  for(i=0;i<i_D_COST->mland->N_int_disc;i++)
	    {
	      if(i_D_COST->mland->lfreq_disc[ind_l+i]>0)
		{
		  lmid = i_D_COST->mland->l_disc[ind_l+i];
		  lend = i_D_COST->mland->int_len_lim[i];
		  lstart = i_D_COST->mland->int_len_lim[i-1];
		  cens_func_old = cens_function(lstart,lmid,lend,k_old,cens_m,cens_r);
		  cens_func_new = cens_function(lstart,lmid,lend,k_new,cens_m,cens_r);
		  if(cens_func_old > G_ZERO)
		    log_disc_old += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_old);
		  if(cens_func_new > G_ZERO)
		    log_disc_new += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_new);
		  N_disc += i_D_COST->mland->lfreq_disc[ind_l+i];
		}
	    }
	  ind_l += i_D_COST->mland->N_int_disc;
	}
    }
  /* log accept prior */
  log_prior_diff = -G_HALF*i_D_COST->cens->tau[0]*(par_new*(par_new-2*i_D_COST->cens->mu[0])-
						   par_old*(par_old-2*i_D_COST->cens->mu[0]));
  d = log_prior_diff;

  /* log accept data */
  d += N_land*(log(k_new)-log(k_old))+log_disc_new-log_disc_old;

  u = genunf(G_ZERO,G_ONE);
  if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
    i_D_COST->cens->k = k_new;

  if(printFile)
    {
      fprintf(fp,"k_old=%f,k_new=%f,x_old=%f,x_new=%f\n",k_old,k_new,par_old,par_new);
      fprintf(fp,"log_lik_old_discard=%f,log_lik_new_discard=%f,N_discard=%d\n",
	      log_disc_old,log_disc_new,N_disc);
      fprintf(fp,"log_lik_old_landed=%f,log_lik_new_landed=%f,N_landed=%d\n",
	      log_land_old,log_land_new,N_land);
      fprintf(fp,"log_lik_prior=%f\n",log_prior_diff);
      fprintf(fp,"log_accept_prob=%f,log(rand_numb)=%f\n\n",d,log(u));
    }
 

  /* Sample m, same for all hauls */
  cens_k = i_D_COST->cens->k;
  par_old = i_D_COST->cens->m;
  par_new = scale_proposal(par_old,fac,NULL);

  /* Find probabilities for old and current sample */
  log_disc_new = G_ZERO;
  log_disc_old = G_ZERO;
  log_land_new = G_ZERO;
  log_land_old = G_ZERO;
  ind_l = 0;
  N_disc = 0;
  N_land = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      cens_r = i_D_COST->cens->r[h];
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      for(f=f_start;f<f_end;f++)
	if(i_D_orig->landed[f]>0)
	  {
	    lstart = log(round(exp(i_D_orig->totlength[f]))-0.5);
	    lmid = i_D_orig->totlength[f];
	    lend = log(round(exp(i_D_orig->totlength[f]))+0.5);
	    
	    cens_func_old = cens_function(lstart,lmid,lend,cens_k,par_old,cens_r);
	    cens_func_new = cens_function(lstart,lmid,lend,cens_k,par_new,cens_r);

	    if(1-cens_func_old > G_ZERO)
	      log_land_old += i_D_orig->landed[f]*log(1-cens_func_old);
	    if(1-cens_func_new > G_ZERO)
	      log_land_new += i_D_orig->landed[f]*log(1-cens_func_new);
	    N_land += i_D_orig->landed[f];
	  }
      for(f=f_start;f<f_end;f++)
	{  
	  if(i_D_orig->discard[f]>0)
	    {
	      lstart = log(round(exp(i_D_orig->totlength[f]))-0.5);
	      lmid = i_D_orig->totlength[f];
	      lend = log(round(exp(i_D_orig->totlength[f]))+0.5);
	    
	      cens_func_old = cens_function(lstart,lmid,lend,cens_k,par_old,cens_r);
	      cens_func_new = cens_function(lstart,lmid,lend,cens_k,par_new,cens_r);

	      if(cens_func_old > G_ZERO)
		log_disc_old += i_D_orig->discard[f]*log(cens_func_old);
	      if(cens_func_new > G_ZERO)
		log_disc_new += i_D_orig->discard[f]*log(cens_func_new);
	      N_disc += i_D_orig->discard[f];
	    }
	}
      // add simulated discards for market landing data
      if(h>=i_D_COST->obs->n_trip)
	{
	  for(i=1;i<i_D_COST->mland->N_int_disc;i++)
	    {
	      if(i_D_COST->mland->lfreq_disc[ind_l+i]>0)
		{
		  lmid = i_D_COST->mland->l_disc[ind_l+i];
		  lend = i_D_COST->mland->int_len_lim[i];
		  lstart = i_D_COST->mland->int_len_lim[i-1];
		  cens_func_old = cens_function(lstart,lmid,lend,cens_k,par_old,cens_r);
		  cens_func_new = cens_function(lstart,lmid,lend,cens_k,par_new,cens_r);
		  if(cens_func_old > G_ZERO)
		    log_disc_old += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_old);
		  if(cens_func_new > G_ZERO)
		    log_disc_new += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_new);
		  N_disc += i_D_COST->mland->lfreq_disc[ind_l+i];
		}
	    }
	  ind_l += i_D_COST->mland->N_int_disc;
	}
    }
  /* log accept prior */
  log_prior_diff = -G_HALF*i_D_COST->cens->tau[1]*(par_new*(par_new-2*i_D_COST->cens->mu[1])-
						   par_old*(par_old-2*i_D_COST->cens->mu[1]));
  d = log_prior_diff;

  /* log accept data */
  d += log_land_new-log_land_old+log_disc_new-log_disc_old;

  u = genunf(G_ZERO,G_ONE);
  if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
    i_D_COST->cens->m = par_new;

  if(printFile)
    {
      fprintf(fp,"m_old=%f,m_new=%f\n",par_old,par_new);
      fprintf(fp,"log_lik_old_discard=%f,log_lik_new_discard=%f,N_discard=%d\n",
	      log_disc_old,log_disc_new,N_disc);
      fprintf(fp,"log_lik_old_landed=%f,log_lik_new_landed=%f,N_landed=%d\n",
	      log_land_old,log_land_new,N_land);
      fprintf(fp,"log_lik_prior=%f\n",log_prior_diff);
      fprintf(fp,"log_accept_prob=%f,log(rand_numb)=%f\n\n",d,log(u));
    }


  /* Sample r parameter */
  ind_l = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      cens_k = i_D_COST->cens->k;
      cens_m = i_D_COST->cens->m;

      par_old = i_D_COST->cens->r[h];
      par_new = scale_proposal(par_old,fac,NULL);

      if(printFile)
	{
	  fprintf(fp,"h=%d\n",h);
	  fprintf(fp,"k=%f,m=%f,r_old=%f,r_new=%f\n",
		  i_D_COST->cens->k,i_D_COST->cens->m,par_old,par_new);
	}
	  
      /* Find probabilities for old and current sample */
      log_disc_new = G_ZERO;
      log_disc_old = G_ZERO;
      log_land_new = G_ZERO;
      log_land_old = G_ZERO;
      N_disc = 0;
      N_land = 0;

      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      

      // landed data
      for(f=f_start;f<f_end;f++)
	if(i_D_orig->landed[f]>0)
	  {
	    lstart = log(round(exp(i_D_orig->totlength[f]))-0.5);
	    lmid = i_D_orig->totlength[f];
	    lend = log(round(exp(i_D_orig->totlength[f]))+0.5);
	    
	    cens_func_old = cens_function(lstart,lmid,lend,cens_k,cens_m,par_old);
	    cens_func_new = cens_function(lstart,lmid,lend,cens_k,cens_m,par_new);

	    if(1-cens_func_old > G_ZERO)
	      log_land_old += i_D_orig->landed[f]*log(1-cens_func_old);
	    if(1-cens_func_new > G_ZERO)
	      log_land_new += i_D_orig->landed[f]*log(1-cens_func_new);
	    N_land += i_D_orig->landed[f];
	  }
      // discard data
      for(f=f_start;f<f_end;f++)
	{  
	  if(i_D_orig->discard[f]>0)
	    {
	      lstart = log(round(exp(i_D_orig->totlength[f]))-0.5);
	      lmid = i_D_orig->totlength[f];
	      lend = log(round(exp(i_D_orig->totlength[f]))+0.5);
	      
	      cens_func_old = cens_function(lstart,lmid,lend,cens_k,cens_m,par_old);
	      cens_func_new = cens_function(lstart,lmid,lend,cens_k,cens_m,par_new);
	      if(cens_func_old > G_ZERO)
		log_disc_old += i_D_orig->discard[f]*log(cens_func_old);
	      if(cens_func_new > G_ZERO)
		log_disc_new += i_D_orig->discard[f]*log(cens_func_new);
	      N_disc += i_D_orig->discard[f];
	    }
	}
      // add simulated discards for market landing data
      if(h>=i_D_COST->obs->n_trip)
	{
	  for(i=1;i<i_D_COST->mland->N_int_disc;i++)
	    {
	      if(i_D_COST->mland->lfreq_disc[ind_l+i]>0)
		{
		  lmid = i_D_COST->mland->l_disc[ind_l+i];
		  lend = i_D_COST->mland->int_len_lim[i];
		  lstart = i_D_COST->mland->int_len_lim[i-1];
		  cens_func_old = cens_function(lstart,lmid,lend,cens_k,cens_m,par_old);
		  cens_func_new = cens_function(lstart,lmid,lend,cens_k,cens_m,par_new);

		  if(cens_func_old > G_ZERO)
		    log_disc_old += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_old);
		  if(cens_func_new > G_ZERO)
		    log_disc_new += i_D_COST->mland->lfreq_disc[ind_l+i]*log(cens_func_new);
		  N_disc += i_D_COST->mland->lfreq_disc[ind_l+i];
		}
	    }
	  ind_l += i_D_COST->mland->N_int_disc;
	}
      /* log accept prior */
      log_prior_diff = -G_HALF*i_D_COST->cens->tau[2]*(par_new*(par_new-2*i_D_COST->cens->mu[2])-
						       par_old*(par_old-2*i_D_COST->cens->mu[2]));
      d = log_prior_diff;

      /* log accept data */
      d += log_land_new-log_land_old+log_disc_new-log_disc_old;

      u = genunf(G_ZERO,G_ONE);
      if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
	i_D_COST->cens->r[h] = par_new;

      if(printFile)
	{
	  fprintf(fp,"r_old=%f,r_new=%f\n",par_old,par_new);
	  fprintf(fp,"log_lik_old_discard=%f,log_lik_new_discard=%f,N_discard=%d\n",
		  log_disc_old,log_disc_new,N_disc);
	  fprintf(fp,"log_lik_old_landed=%f,log_lik_new_landed=%f,N_landed=%d\n",
		  log_land_old,log_land_new,N_land);
	  fprintf(fp,"log_lik_prior=%f,mu=%f,tau=%f\n",log_prior_diff,
		  i_D_COST->cens->mu[2],i_D_COST->cens->tau[2]);
	  fprintf(fp,"log_accept_prob=%f,log(rand_numb)=%f\n",d,log(u));
	}
    }
  
  if(printFile)
    fclose(fp);
  
  return(0);
}               /* end of sample_cens_par */


int sample_ages_init2_COST(Data_orig *i_D_orig,Age_struct *i_age,
			   Data_age *i_D_age,LW_struct *i_length,
			   Data_lin *i_D_lga,Data_g_a *i_D_g_a,
			   Data_COST *i_D_COST,int saveSim)
{
  int a,a2,h,season,f,f_start,f_end,ind_a;
  double age,lobs,lmin,lmax,Int,Slp,log_a0,log_adiff;

  FILE *fp;
  if(saveSim)
    fp = fopen("ages_miss_start.txt","w");

  /* Initialize */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	  i_D_age->Ages_disc[h][a] = 0;
	  i_D_age->Ages_land[h][a] = 0;
	}
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
          i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
          i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	}
    }

  Int = i_length->par->eff[0][0][0][0];
  Slp = i_length->par->eff[0][1][0][0];
  log_a0 = log(i_D_g_a->a_vec[0]);
  log_adiff = log(i_D_g_a->a_vec[i_D_g_a->ncat-1])-log(i_D_g_a->a_vec[0]);

  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];
      f_start = i_D_orig->start_noAge[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      for(f=f_start;f<f_end;f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[f];
	  if(a < -1000 && lobs > -1000.0)
	    {
	      age = exp(log_a0 + log_adiff*(lobs-Int)/Slp);
	      if(age < i_D_age->a_vec[0])
		age = i_D_age->a_vec[0];
	      if(age > i_D_age->a_vec[i_D_age->glm->ncat-1])
		age = i_D_age->a_vec[i_D_age->glm->ncat-1];		
	      a2 = (int) age;
	      i_D_age->Ages[h][a2] += i_D_orig->replength[f];
	      i_D_age->Ages_disc[h][a2] += i_D_orig->discard[f];
	      i_D_age->Ages_land[h][a2] += i_D_orig->landed[f];
	      ind_a = i_D_g_a->a2Age_vec[a2]+(season-1);
	      i_D_lga->Ages[h][ind_a] += i_D_orig->replength[f];
	      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) i_D_orig->replength[f];
	      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) i_D_orig->replength[f];
	      if(saveSim)		
		fprintf(fp,"%d %d %d %f %f %d %d\n",h,a2,season,i_D_g_a->a_vec[ind_a],
			lobs,i_D_orig->discard[f],i_D_orig->landed[f]);
	    }
	}
    }

  if(saveSim)
    fclose(fp);


  return(0);
}		/* end of sample_ages_init2_COST */


/*!
  \brief Sample ages based on empirical age-given-length distribution

  This routine is only used for initialization of the parameters to be simulated.
*/
int sample_ages_init_COST(Data_orig *i_D_orig,Age_struct *i_age,
			  Data_age *i_D_age,LW_struct *i_length,
			  Data_lin *i_D_lga,Data_g_a *i_D_g_a,
			  Data_COST *i_D_COST,int saveSim)
{
  int            a,a2,f,f_start,f_end,h,i,s,l_int,ind_a,N_int,season,nSeason;
  long          *ages;
  double       **P_al_disc,***P_al_disc_s,**P_al_land,***P_al_land_s,**P_al,***P_al_s;
  double         u,sum,sigma,lobs,fac;
  int Ndisc=0,Nland=0;
  
  FILE *fp;
  if(saveSim)
    fp = fopen("ages_miss_start.txt","w");
  
  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  nSeason = i_D_g_a->nSeason;
  
  /* Initialize */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	  i_D_age->Ages_disc[h][a] = 0;
	  i_D_age->Ages_land[h][a] = 0;
	}
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
          i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
          i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	}
    }
  
  /* First making a transition matrix P(a|l) for a finite number of lengths */
  N_int = i_D_orig->n_int_len;
  P_al_disc = Mmatrix_2d(0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);
  P_al_disc_s = Mmatrix_3d(0,nSeason-1,0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok
  P_al_land = Mmatrix_2d(0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);
  P_al_land_s = Mmatrix_3d(0,nSeason-1,0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok
  P_al = Mmatrix_2d(0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);
  P_al_s = Mmatrix_3d(0,nSeason-1,0,N_int-1,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok

  for(i=0;i<N_int;i++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      P_al[i][a] = G_ZERO;

  for(s=0;s<nSeason;s++)
    for(i=0;i<N_int;i++)
      for(a=0;a<i_D_age->glm->ncat;a++)
	P_al_s[s][i][a] = G_ZERO;

  for(i=0;i<N_int;i++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      {
	P_al_disc[i][a] = G_ZERO;
	P_al_land[i][a] = G_ZERO;
      }
  for(s=0;s<nSeason;s++)
    for(i=0;i<N_int;i++)
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  P_al_disc_s[s][i][a] = G_ZERO;
	  P_al_land_s[s][i][a] = G_ZERO;
	}

  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      for(f=f_start;f<f_end;f++)
	{
	  Ndisc += i_D_orig->discard[f];
	  Nland += i_D_orig->landed[f];
	}
    }
  fac = (double) Nland/Ndisc;
  //printf("Ndisc=%d,Nland=%d,fac=%f\n",Ndisc,Nland,fac);

  // P(a|l) for discards - using observer data:
  for(h=0;h<i_D_COST->obs->n_trip;h++)
    {
      season = i_D_orig->season[h];
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h];
      for(f=f_start;f<f_end;f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[f];
	  if(a > -1000 && lobs > -1000.0)
	    {
	      i_D_age->Ages_disc[h][a] += i_D_orig->discard[f];
	      i_D_age->Ages_land[h][a] += i_D_orig->landed[f];
	      l_int = 0;
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      P_al_disc[l_int][a] += (double) i_D_orig->discard[f];
	      P_al_disc_s[season-1][l_int][a] += (double) i_D_orig->discard[f];
	      P_al[l_int][a] += (double) fac* i_D_orig->replength[f];
	      P_al_s[season-1][l_int][a] += (double) fac*i_D_orig->replength[f];
	    }
	}
    }
  // P(a|l) for landed - using market landing data:
  for(h=i_D_COST->obs->n_trip;h<i_D_lga->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];
      f_start = i_D_orig->start_Age[h];
      f_end = i_D_orig->start_noAge[h];
      for(f=f_start;f<f_end;f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[f];
	  if(a > -1000 && lobs > -1000.0)
	    {
	      i_D_age->Ages_disc[h][a] += i_D_orig->discard[f];
	      i_D_age->Ages_land[h][a] += i_D_orig->landed[f];
	      l_int = 0;
	      //while(lobs > i_D_orig->int_len_lim[l_int+1])
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      P_al_land[l_int][a] += i_D_orig->landed[f];
	      P_al_land_s[season-1][l_int][a] += (double) i_D_orig->landed[f];
	      P_al[l_int][a] += (double) i_D_orig->replength[f];
	      P_al_s[season-1][l_int][a] += (double) i_D_orig->replength[f];
	    }
	}
    }
  // Convert to probabilities
  for(i=0;i<N_int;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  sum += P_al_disc[i][a];
	}
      if(sum < 0.0001)
	{
	  if(i==0)
	    P_al_disc[i][0] = G_ONE;
	  else
	    for(a=0;a<i_D_age->glm->ncat;a++)
	      P_al_disc[i][a] = P_al_disc[i-1][a];
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al_disc[i][a] /= sum;
	}
    } 
  for(i=0;i<N_int;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  sum += P_al_land[i][a];
	}
      if(sum < 0.0001)
	{
	  if(i==0)
	    P_al_land[i][0] = G_ONE;
	  else
	    for(a=0;a<i_D_age->glm->ncat;a++)
	      P_al_land[i][a] = P_al_land[i-1][a];
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al_land[i][a] /= sum;
	}
    }
  for(i=0;i<N_int;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	sum += P_al[i][a];
      if(sum < 0.0001)
	{
	  if(i==0)
	    P_al[i][0] = G_ONE;
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al[i][a] = P_al[i-1][a];
	    }
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al[i][a] /= sum;
	}
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_s[s][i][a];
	  if(sum<0.0001)//use P_al (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] = P_al[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] /= sum;
	    }
	}
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_disc_s[s][i][a];
	  if(sum<0.0001)//use P_al_disc (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_disc_s[s][i][a] = P_al_disc[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_disc_s[s][i][a] /= sum;
	    }
	}
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_land_s[s][i][a];
	  if(sum<0.0001)//use P_al_land (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_land_s[s][i][a] = P_al_land[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_land_s[s][i][a] /= sum;
	    }
	}
    }
  // Start simulation of missing ages
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];
      f_start = i_D_orig->start_noAge[h];
      f_end = i_D_orig->start_noAge[h]+i_D_orig->num_noAge[h];
      for(f=f_start;f<f_end;f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  lobs = i_D_orig->totlength[f];
	  if(a < -1000 && lobs > -1000.0)
	    {
	      //Find right interval
	      l_int = 0;
	      //while(lobs > i_D_orig->int_len_lim[l_int+1])
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      //Simulate age
	      //Use P_al_disc if discarded and P_al_land if landed
	      my_genmul(i_D_orig->discard[f],P_al_disc_s[season-1][l_int],
			i_D_age->glm->ncat,ages);
	      for(a2=0;a2<i_D_age->glm->ncat;a2++)
		{
		  if(ages[a2]>0)
		    {
		      i_D_age->Ages[h][a2] += (int) ages[a2];
		      i_D_age->Ages_disc[h][a2] += (int) ages[a2];
		      ind_a = i_D_g_a->a2Age_vec[a2]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a2];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a2];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a2];
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f %d 0\n",h,a2,season,i_D_g_a->a_vec[ind_a],
				lobs,(int) ages[a2]);
		    }
		}
	      my_genmul(i_D_orig->landed[f],P_al_land_s[season-1][l_int],
			i_D_age->glm->ncat,ages);
	      for(a2=0;a2<i_D_age->glm->ncat;a2++)
		{
		  if(ages[a2]>0)
		    {
		      i_D_age->Ages[h][a2] += (int) ages[a2];
		      i_D_age->Ages_land[h][a2] += (int) ages[a2];
		      ind_a = i_D_g_a->a2Age_vec[a2]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a2];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a2];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a2];
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f 0 %d\n",h,a2,season,i_D_g_a->a_vec[ind_a],
				lobs,(int) ages[a2]);
		    }
		}
	    }
	}
    }
 
  if(saveSim)
    fclose(fp);

  // Free allocated memory
  FREE(ages);
  Fmatrix_2d(&P_al_disc[0][0],&P_al_disc[0]);
  Fmatrix_3d(&P_al_disc_s[0][0][0],&P_al_disc_s[0][0],&P_al_disc_s[0]);
  Fmatrix_2d(&P_al_land[0][0],&P_al_land[0]);
  Fmatrix_3d(&P_al_land_s[0][0][0],&P_al_land_s[0][0],&P_al_land_s[0]);
  Fmatrix_2d(&P_al[0][0],&P_al[0]);
  Fmatrix_3d(&P_al_s[0][0][0],&P_al_s[0][0],&P_al_s[0]);

  return(0);
}		/* end of sample_ages_init_COST */

/*!
  \author Geir Storvik
  \brief Samples missing ages.

  Here all data are assumed to be of the amigo type, that is long strings of
  age and length. 

  In order to speed up computation, the fish are assumed ordered in hauls. 
  Further, inside each haul, the lenghts are assumed ordered in increasing values. 
  Then the length-range is divided into a finite
  number of intervals in which the age-probabilities are assumed constant for
  all length-values inside an interval.
*/
int sample_ages_COST(Data_orig *i_D_orig,Age_struct *i_age,Data_age *i_D_age,
		     LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		     int saveSim, int i_it)
{
  int            a,f,h,i,ind_a,ind_f,ind_miss,n,cum_fish,aobs,season;
  int            l_int,l_int_prev;
  double         lobs,r;
  double         sum_p,sigma;
  double         lstart,lend;
  Data_glm      *glm;
  double        *p, *p2, *mu, *beta;
  long          *ages;
  long          *test_ages, *test_ages2;
  

  glm = i_D_lga->glm;

  p = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  p2 = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);      // Free ok

  test_ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  test_ages2 = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  beta = CALLOC(i_D_lga->glm->nxcov,double);      // Free ok      

  sigma = G_ONE/sqrt(i_length->par->tau_obs);

  #ifdef DEBUG_TEST
  FILE *unit;
  unit = fopen("Ages_samples.dat","w");
  #endif
  
  FILE *fp;
  if(saveSim)
    fp = fopen("ages_miss.txt","w");

  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      test_ages[a] = 0;
      test_ages2[a] = 0;
    }
  
  cum_fish = 0;
  ind_miss = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      season = i_D_orig->season[h];

      /* Start by initializing sufficient statistics */ 
      if(i_D_age->type_age[0]==1) /* Age errors: All simulated */
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      i_D_age->Ages[h][a] = 0;
	      i_D_age->Ages_disc[h][a] = 0;
	      i_D_age->Ages_land[h][a] = 0;
	    }
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      i_D_lga->Ages[h][a] = 0;
	      i_D_lga->sum_by_cat[h][a] = G_ZERO;
	      i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	    }
	}
      else  /* No age errors: Keep aged fish */
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	      i_D_age->Ages_disc[h][a] = 0;
	      i_D_age->Ages_land[h][a] = 0;
	    }
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
	      i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
	      i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	    }
	}

      /* find prior age-probabilities */
      sum_p = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] = exp(i_age->alpha[h][a]);
	  sum_p += p[a];
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  p[a] /= sum_p;
	}
      /* Find intercept and slope for lga model */
      for(i=0;i<i_D_lga->glm->nxcov;i++)
	beta[i] = calc_eff(i_D_lga->glm->xcov[i],i_length->par->eff[0][i],h);
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
	}
      if(i_D_age->type_age[0]==1)/*Ages observed with error*/
	{
	  /* Both age and length observed */
	  for(f=0;f<(i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h]);f++)
	    {
	      ind_f = cum_fish+f;
	      aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	      lobs = i_D_orig->totlength[ind_f];
	      
	      /* Calculate probabilities */
	      sum_p = G_ZERO;
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		{
		  a = i_age->A_neigh[aobs][i];
		  p2[i] = p[a]*i_age->A2A[a][aobs]*dnorm(lobs,mu[a],sigma);
		  sum_p += p2[i];
		}
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		p2[i] /= sum_p;
	      my_genmul(i_D_orig->replength[ind_f],p2,i_age->A_Nneigh[aobs],ages);
	      for(i=0;i<i_age->A_Nneigh[aobs];i++)
		{
		  a = i_age->A_neigh[aobs][i];
		  i_D_age->Ages[h][a] += (int) ages[i];
		  
		  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		  i_D_lga->Ages[h][ind_a] += (int) ages[i];
		  test_ages2[a] += ages[i];
		  /* Update sufficient statistics */
		  r = (double) ages[i];
		  if(r< 0)
		    {
		      write_warning("Negative ages generated");
		      return(1);
		    }
		  i_D_lga->sum_by_cat[h][ind_a] += r*lobs;
		  i_D_lga->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
		}
	    }
	}
      /* Loop through all aged fish */
      for(f=i_D_orig->start_Age[h];f< i_D_orig->start_noAge[h];f++)
	{
	  a = i_D_orig->totage[f]-i_D_age->a_vec[0];
	  i_D_age->Ages_disc[h][a] += i_D_orig->discard[f];
	  i_D_age->Ages_land[h][a] += i_D_orig->landed[f];
	}

      ind_f = i_D_orig->start_noAge[h];
      l_int_prev = -1;
      l_int = 0;
      /* Loop through all non-aged fish */
      for(f=0;f<i_D_orig->num_noAge[h];f++)
	{
	  if(i_D_orig->discard[ind_f+f]>0)
	    {
	      lobs = i_D_orig->totlength[ind_f+f];
	      lstart = log(round(exp(lobs))-0.5);
	      lend = log(round(exp(lobs))+0.5);
	      //     while(lobs > i_D_orig->int_len_lim[l_int+1])
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      if(l_int!=l_int_prev)
		{
		  //New calculation of age-probabilities needs to be performed
		  /* Calculate probabilities */
		  sum_p = G_ZERO;
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    {
		      //p2[a] = p[a]*dnorm(lobs,mu[a],sigma);
		      p2[a] = p[a]*(pnorm(lend,mu[a],sigma)-pnorm(lstart,mu[a],sigma));
		      sum_p += p2[a];
		    }
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    p2[a] /= sum_p;
		  l_int_prev = l_int;
		}
	      my_genmul(i_D_orig->discard[ind_f+f],p2,i_D_age->glm->ncat,ages);
	      n = 0;
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  if(ages[a]>0)
		    {
		      i_D_age->Ages[h][a] += (int) ages[a];
		      i_D_age->Ages_disc[h][a] += (int) ages[a];
		      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
		      test_ages[a] += ages[a];
		      n+= ages[a];
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],lobs,(int) ages[a]);
		    }

                  #ifdef DEBUG_TEST
		  fprintf(unit,"%d %d %d %lf\n",h,f,a,lobs);
                  #endif
		}
	    }
	  if(i_D_orig->landed[ind_f+f]>0)
	    {
	      if(i_D_orig->discard[ind_f+f]==0)
		{
		  lobs = i_D_orig->totlength[ind_f+f];
		  lstart = log(round(exp(lobs))-0.5);
		  lend = log(round(exp(lobs))+0.5);
		  //     while(lobs > i_D_orig->int_len_lim[l_int+1])
		  while(lobs > i_D_orig->int_len_lim[l_int])
		    l_int++;
		  if(l_int!=l_int_prev)
		    {
		      //New calculation of age-probabilities needs to be performed
		      /* Calculate probabilities */
		      sum_p = G_ZERO;
		      for(a=0;a<i_D_age->glm->ncat;a++)
			{
			  //p2[a] = p[a]*dnorm(lobs,mu[a],sigma);
			  p2[a] = p[a]*(pnorm(lend,mu[a],sigma)-pnorm(lstart,mu[a],sigma));
			  sum_p += p2[a];
			}
		      for(a=0;a<i_D_age->glm->ncat;a++)
			p2[a] /= sum_p;
		      l_int_prev = l_int;
		    }
		}
	      my_genmul(i_D_orig->landed[ind_f+f],p2,i_D_age->glm->ncat,ages);
	      n = 0;
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  if(ages[a]>0)
		    {
		      i_D_age->Ages[h][a] += (int) ages[a];
		      i_D_age->Ages_land[h][a] += (int) ages[a];
		      ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs * (double) ages[a];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs * (double) ages[a];
		      test_ages[a] += ages[a];
		      n+= ages[a];
		      if(saveSim)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,
				i_D_g_a->a_vec[ind_a],lobs,(int) ages[a]);
		    }

                  #ifdef DEBUG_TEST
		  fprintf(unit,"%d %d %d %lf\n",h,f,a,lobs);
                  #endif
		}
	    }
	}
      ind_miss += i_D_g_a->ncat*(i_D_orig->n_int_len+1);
      cum_fish += i_D_orig->nFishBoat[h];
    }
  if(saveSim)
    fclose(fp);

  #ifdef AGE_TEST
  fprintf(stderr,"\nobs:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(stderr,"%5d ",test_ages2[a]);
  fprintf(stderr,"\nunobs: ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(stderr,"%5d ",test_ages[a]);
  FILE *ageFile;
  ageFile = fopen("../Ages_sampled.dat","a");
  fprintf(ageFile,"\nobs:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(ageFile,"%3d:%5d ",a+1,test_ages2[a]);
  fprintf(ageFile,"\nunobs: ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    fprintf(ageFile,"%3d:%5d ",a+1,test_ages[a]);
  double testSum;
  fprintf(ageFile,"\nsum:   ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      testSum = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  //testSum += i_D_age->Ages[h][a];
	  testSum += i_D_lga->Ages_fix[h][a];
	}
      fprintf(ageFile,"%3d:%5.0f ",a+1,testSum);
    }
  fclose(ageFile);
  #endif 

  #ifdef DEBUG_TEST
  fclose(unit);
  #endif


  // Free memory allocated in this routine
  FREE(p);
  FREE(p2);
  FREE(mu);
  FREE(ages);
  FREE(test_ages);
  FREE(test_ages2);
  FREE(beta);

  return(0);
}		/* end of sample_ages_COST */

