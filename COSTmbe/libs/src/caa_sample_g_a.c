/*!
  \file caa_sample_g_a.c
  \author Geir Storvik
  \brief Routines for sampling nonlinear part of lga model

*/
#include "caa.h"


static double calc_S_unorm(double *i_par);
static double calc_S_beta(double *i_par);
static double calc_S_unorm0(double *i_par);
static double calc_S_beta0(double *i_par);
static int find_opt_g_a(double *o_par_opt,double *o_g);
static int calc_g_a_S_R(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int calc_g_a_S_R_unorm(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int prior_g_a_S_R(int i_npar,double *i_par);
static int calc_g_a_polyn3(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int calc_g_a_polyn3_unorm(int i_ncat,double *i_a_vec,
                                 double *i_par,double *o_g);
static int (*s_calc_g_a)(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int (*s_calc_g_a_unorm)(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int (*s_prior_g_a_par)(int i_npar,double *i_par);

static int       s_ncat,s_npar;
static double   *s_a_vec;
static double   *s_parmin, *s_parmax;
static double   *s_par, *s_par_opt, *s_par_new, *s_p, **s_xi;
static double  **s_suff;
static double   *s_g_a;

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif


/*!
  \author Geir Storvik
  \brief Allocates space for routines used to sample non-linear part of lga-model

  This routine must be called before using routines in the caa_sample_g_a.c file.
*/
int sample_g_a_initialize(int i_ncat,int i_g_a_model)
{
  if(i_g_a_model==1)
    {
      s_calc_g_a = calc_g_a_S_R;
      s_calc_g_a_unorm = calc_g_a_S_R_unorm;
      s_prior_g_a_par = prior_g_a_S_R;
      s_npar = 3;
    }
  else if(i_g_a_model==2)
    {
      s_calc_g_a = calc_g_a_polyn3;
      s_calc_g_a_unorm = calc_g_a_polyn3_unorm;
      s_npar = 2;
    }
  s_par = CALLOC(s_npar,double);    // Free ok
  if(!s_par)
    {
      write_warning("sample_g_a_initialize:Error allocating s_par\n");
      return(1);
    }
  s_par_new = CALLOC(s_npar,double);// Free ok
  if(!s_par_new)
    {
      write_warning("sample_g_a_initialize:Error allocating s_par_new\n");
      return(1);
    }
  s_par_opt = CALLOC(s_npar,double);// Free ok
  if(!s_par_opt)
    {
      write_warning("sample_g_a_initialize:Error allocating s_par_opt\n");
      return(1);
    }
  s_p = CALLOC(s_npar+1,double);     // Free ok
  if(!s_p)
    {
      write_warning("sample_g_a_initialize:Error allocating s_p\n");
      return(1);
    }
  s_xi = Mmatrix_2d(0,s_npar,0,s_npar,sizeof(double),1);  // Free ok
  if(!s_xi)
    {
      write_warning("sample_g_a_initialize:Error allocating s_xi\n");
      return(1);
    }
  s_g_a = CALLOC(i_ncat,double);    // Free ok
  if(!s_g_a)
    {
      write_warning("sample_g_a_initialize:Error allocating s_g_a\n");
      return(1);
    }

  s_parmin = CALLOC(s_npar,double);
  s_parmax = CALLOC(s_npar,double);
  s_parmin[0] = 0.000000001;
  s_parmax[0] = 7.0;
  s_parmin[1] = 0.00000001;;
  s_parmax[1] = 1000000000.00;
  s_parmin[2] = 0.00000001;;
  s_parmax[2] = 1000000000.00;

  return(0);
}		/* end of sample_g_a_initialize */




/*!
  \author Geir Storvik
  \brief Reallocate space allocated by sample_g_a_initialize
*/
int sample_g_a_re_initialize()
{
  FREE(s_par);
  FREE(s_par_new);
  FREE(s_par_opt);
  FREE(s_p);
  Fmatrix_2d(&s_xi[0][0],&s_xi[0]);
  FREE(s_g_a);
  FREE(s_parmin);
  FREE(s_parmax);

  return(0);
}		/* end of sample_g_a_re_initialize */


/*!
  \author Geir Storvik
  \brief Calculates sufficient statistics for sampling the function \f$g(a)\f$

  The function \f$g(a)\f$ describes the non-linear link between age and length
  So far the only possible choises are log-linear and the Schnute-Richards
  models. 

Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int suff_g_a_init(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		  Data_g_a *i_D_g_a,int i_start_h,double **o_suff)
{
  int         a,h;
  double      N_h_a;
  Data_glm   *glm;

  
  glm = i_D_lga->glm;
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      o_suff[0][a] = G_ZERO;
      o_suff[1][a] = G_ZERO;
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
        {
          N_h_a = (double) i_D_lga->Ages[h][a];
	  o_suff[0][a] += i_D_lga->sum_by_cat[h][a];
	  o_suff[1][a] += N_h_a;
	}
    }

  return(0);
}               /* end of suff_g_a_init */



/*L:suff_g_a*
________________________________________________________________

		suff_g_a
________________________________________________________________

Name:		suff_g_a
Syntax:		
Description:    Calculates sufficient statistics for sampling the function g(a) 
                describing the non-linear
                link between age and length
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int suff_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
	     Data_g_a *i_D_g_a,int i_start_h,double **o_suff)
{
  int         a,h,i;
  double      N_h_a, res;
  double     *beta;
  Data_glm   *glm;

  glm = i_D_lga->glm;
  beta = CALLOC(glm->nxcov,double); //Free ok

  if(glm->nxcov==2)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	    {
	      N_h_a = (double) i_D_lga->Ages[h][a];
	      for(i=0;i<glm->nxcov;i++)
		beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	      res = i_D_lga->sum_by_cat[h][a] - beta[0]*N_h_a;
	      o_suff[0][a] += beta[1]*res;
	      o_suff[1][a] += beta[1]*beta[1]*N_h_a;
	    }
	}
    }
  else if(glm->nxcov==3)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	    {
	      N_h_a = (double) i_D_lga->Ages[h][a];
	      for(i=0;i<glm->nxcov;i++)
		beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	      res = i_D_lga->sum_by_cat[h][a] - beta[0]*N_h_a - beta[2]*i_D_lga->haulweight[h];
	      o_suff[0][a] += beta[1]*res;
	      o_suff[1][a] += beta[1]*beta[1]*N_h_a;
	    }
	}
    }

  FREE(beta);

  return(0);
}               /* end of suff_g_a */



/*L:fit_g_a_init*
________________________________________________________________

		fit_g_a_init
________________________________________________________________

Name:		fit_g_a_init
Syntax:		
Description:    Samples the function g(a) describing the non-linear
                link between age and length, assuming g(a)=log(1-exp(-K(a-a_0)))
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int fit_g_a_init(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		 Data_g_a *i_D_g_a,double **i_suff)
{
  int       a,err,i,j,k,p,iter;
  double   *par_opt, **xi;
  double    fret,ftol,S,minS,S0;

  s_suff = i_suff;
  s_ncat = i_D_g_a->ncat;
  s_a_vec = i_D_g_a->a_vec;

  p = i_D_g_a->g_a_npar;
  par_opt = CALLOC(p+1,double);
  xi = Mmatrix_2d(0,p,0,p,sizeof(double),1);

  for(j=0;j<i_D_g_a->g_a_npar;j++)
    par_opt[1+j] = log(max(0.00001,i_D_g_a->g_a_par[j]));


  for(j=1;j<=p;j++)
  for(i=1;i<=p;i++)
    xi[j][i] = G_ZERO;
  for(j=1;j<=p;j++)
    xi[j][j] = 1.0;

  ftol = 0.000001;

  //par_opt[3] = 1.0;

  // First try log model
  par_opt[1] = 0.0000000000000001;
  for(i=2;i<=s_npar;i++)
    par_opt[i] = G_ONE;
  S0 = calc_S_beta0(par_opt);
  i_D_g_a->g_a_par[0] = G_ZERO;
  for(j=1;j<s_npar;j++)
    i_D_g_a->g_a_par[j] = exp(par_opt[1+j]);

  minS = S0;
  for(k=0;k<10;k++)
    {  
      for(j=0;j<i_D_g_a->g_a_npar;j++)
        par_opt[1+j] = log(max(0.00001,i_D_g_a->g_a_par[j]));
      par_opt[1] = (double) k/4.0 + 0.00000001;
      err = powell(par_opt,xi,p,ftol,&iter,&fret,calc_S_beta);
      S = calc_S_beta(par_opt);
      if(S<minS)
	{
	  for(j=0;j<s_npar;j++)
	    i_D_g_a->g_a_par[j] = exp(par_opt[1+j]);
	}
    }


  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);

  fprintf(g_caa_log, "fit_g_a_init:\n");
  fprintf(g_caa_log, "amin = %d  amax = %d  c = %lf  theta = %lf  gamma = %lf \n", 
	  i_D_g_a->a_vec[0], i_D_g_a->a_vec[i_D_g_a->ncat-1], i_D_g_a->g_a_par[0], 
	  i_D_g_a->g_a_par[1], i_D_g_a->g_a_par[2]);

  for(a=0;a<i_D_g_a->ncat;a++)
    i_D_g_a->g_a[a] = s_g_a[a];

  #ifdef DEBUG_G_A
  printf("g_a=");
  for(a=0;a<i_D_g_a->ncat;a++)
    printf("%lf, ",s_g_a[a]);
  printf("\n");
 for(j=0;j<s_npar;j++)
   par_opt[1+j] = log(i_D_g_a->g_a_par[j]);
  S = calc_S_beta(par_opt);
  printf("c=%lg,theta=%lg,gamma=%lg, S=%lg\n",
         i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],
         i_D_g_a->g_a_par[2],S);
  #endif


  FREE(par_opt);
  Fmatrix_2d(&xi[0][0],&xi[0]);

  return(0);
}               /* end of fit_g_a_init */




/*L:sample_g_a*
________________________________________________________________

		sample_g_a
________________________________________________________________

Name:		sample_g_a
Syntax:		
Description:    Samples the function g(a) describing the non-linear
                link between age and length, assuming g(a)=log(1-exp(-K(a-a_0)))

                The algorithm is a systematic scan Metropolis-Hastings algorithm
	        where for c and theta a scale-proposal

                Because ths sampling takes very little time, num=1000 iteration are
                performe in order to get reasonable estimates from the 
		conditional posterior. This could probably be made more
		efficient. 
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int sample_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
               Data_g_a *i_D_g_a,double **i_suff,int i_it)
{
  int       p,h;
  int       num=1000;
  int       a,i,i_new,j,k, *acc_par,err,N_c;
  double    d=0,fac,sig,max_c,**Q, mu[2],beta[2],N_h_a,eps[2],eps2[2],bvec[2];
  double    S, S_new,u,tau;
  double    logprior,logprior_new;
  double    *N_a, **N_a_beta, *sum_by_cat, *sum_beta1_by_cat, **suff;
  double    aprior[3],bprior[3];
  Data_glm  *glm;

  Q = Mmatrix_2d(0,1,0,1,sizeof(double),1);

  s_suff = i_suff;
  s_ncat = i_D_g_a->ncat;
  s_a_vec = i_D_g_a->a_vec;
  glm = i_D_lga->glm;
  tau = i_length->par->tau_obs;
  aprior[0] = 0.01;
  aprior[1] = 0.01;
  aprior[2] = 0.01;
  bprior[0] = 0.01;
  bprior[1] = 0.01;
  bprior[2] = 0.01;


  if(0)
    {  
      for(j=0;j<i_D_g_a->g_a_npar;j++)
	s_par_opt[j] = max(0.00001,i_D_g_a->g_a_par[j]);
      err = find_opt_g_a(s_par_opt,&S);
      for(j=0;j<i_D_g_a->g_a_npar;j++)
	i_D_g_a->g_a_par[j] = s_par_opt[j];
      err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
      for(a=0;a<i_D_g_a->ncat;a++)
	i_D_g_a->g_a[a] = s_g_a[a];
      #ifdef DEBUG_G_A
      printf("c=%lg,theta=%lg,gamma=%lg,S=%lg\n",
	     s_par_opt[0],s_par_opt[1],
	     s_par_opt[2],S);
      printf("g_a=\n");
      for(a=0;a<i_D_g_a->ncat;a++)
	printf("%lf %lf\n",s_g_a[a],s_suff[0][a]/s_suff[1][a]);
      #endif  
      return(0);
    }

  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_g_a: Starting sample_g_a\n");
  #endif

  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
  S = G_ZERO;
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      if(s_suff[1][a]>0)
	S += s_suff[1][a]*(s_g_a[a]-s_suff[0][a]/s_suff[1][a])*
			    (s_g_a[a]-s_suff[0][a]/s_suff[1][a]);
    }
  logprior = G_ZERO;
  for(i=0;i<3;i++)
    logprior += -(aprior[i]-G_ONE)*log(i_D_g_a->g_a_par[i])-bprior[i]*i_D_g_a->g_a_par[i];


  // Initialization for sampling many times
  N_a = CALLOC(i_D_g_a->ncat,double);
  N_a_beta = Mmatrix_2d(0,i_D_g_a->ncat-1,0,3,sizeof(double),1);
  sum_by_cat = CALLOC(i_D_g_a->ncat,double);
  sum_beta1_by_cat = CALLOC(i_D_g_a->ncat,double);
  suff = Mmatrix_2d(0,1,0,i_D_g_a->ncat-1,sizeof(double),1);
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      N_a[a] = G_ZERO;
      sum_by_cat[a] = G_ZERO;
      N_a_beta[a][0] = G_ZERO;
      N_a_beta[a][1] = G_ZERO;
      N_a_beta[a][2] = G_ZERO;
      N_a_beta[a][3] = G_ZERO;
      sum_beta1_by_cat[a] = G_ZERO;
    }
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(i=0;i<glm->nxcov-1;i++)
	{
	  beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	  beta[i] -= i_length->par->eff[0][i][0][0];
	}
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  N_h_a = (double) i_D_lga->Ages[h][a];
	  N_a[a] += N_h_a;
	  sum_by_cat[a] += i_D_lga->sum_by_cat[h][a];
          N_a_beta[a][0] += N_h_a*beta[0];
          N_a_beta[a][1] += N_h_a*beta[1];
          N_a_beta[a][2] += N_h_a*beta[0]*beta[1];
          N_a_beta[a][3] += N_h_a*beta[1]*beta[1];
	  sum_beta1_by_cat[a] += beta[1]*i_D_lga->sum_by_cat[h][a];
	}
    }

  /* Run Metropolis-Hastings for parameters */


  #ifdef DEBUG_G_A
  printf("c=%lg,theta=%lg,gamma=%lg,S=%lg\n",
         s_par_opt[0],s_par_opt[1],
         s_par_opt[2],S);
  printf("g_a=\n");
  for(a=0;a<i_D_g_a->ncat;a++)
    printf("%lf %lf\n",s_g_a[a],s_suff[0][a]/s_suff[1][a]);
  #endif

  acc_par = CALLOC(3,int);
  max_c = 1.0;
  N_c = 200;
  fac = 1.1;
  sig = 0.3;
  i = (int) (i_D_g_a->g_a_par[2] * (double) N_c/max_c);
  for(k=0;k<num;k++)
    {
      //err= suff_g_a(i_age,i_length,i_D_age,i_D_lga,0,i_suff);
      beta[0] = i_length->par->eff[0][0][0][0];
      beta[1] = i_length->par->eff[0][1][0][0];
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_suff[0][a] = sum_beta1_by_cat[a]-N_a_beta[a][2]+
	    beta[1]*(sum_by_cat[a]-N_a_beta[a][0])-
	    beta[0]*N_a_beta[a][1]-
	    beta[0]*beta[1]*N_a[a];
	  i_suff[1][a] = beta[1]*beta[1]*N_a[a]+G_TWO*beta[1]*N_a_beta[a][1]+
	    N_a_beta[a][3];
	}
      // Sample parameters i g-function
      for(p=0;p<3;p++)
	{
          #ifdef DEBUG_G_A2
	  printf("old:c=%lg,theta=%lg,gamma=%lg,S=%lf\n",
		 i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],i_D_g_a->g_a_par[2],S);
          #endif
	  if(p==0) /* sample c */
	    {
	      s_par_new[0] = scale_proposal(i_D_g_a->g_a_par[0],fac,NULL);
	      s_par_new[1] = i_D_g_a->g_a_par[1];
	      i_new = i;
	      s_par_new[2] = i_D_g_a->g_a_par[2];
	    }
	  else if(p==1) /* sample theta */
	    {
	      s_par_new[0] = i_D_g_a->g_a_par[0];
	      s_par_new[1] = scale_proposal(i_D_g_a->g_a_par[1],fac,NULL);
	      //if(s_par_new[1]>1e13)
	      //	s_par_new[1] = 1e13;
	      i_new = i;
	      s_par_new[2] = i_D_g_a->g_a_par[2];
	    }
	  else /* sample gamma */
	    {
	      s_par_new[0] = i_D_g_a->g_a_par[0];
	      s_par_new[1] = i_D_g_a->g_a_par[1];
	      //i_new = ignuin(max(0,i-2),min(N_c,i+2));
	      //s_par_new[2] = (double) i_new * max_c / (double) N_c;
	      s_par_new[2] = scale_proposal(i_D_g_a->g_a_par[2],fac,NULL);
	    }
	  S_new = G_ZERO;
	  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,s_par_new,s_g_a);
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      if(s_suff[1][a]>0)
	      S_new += s_suff[1][a]*(s_g_a[a]-s_suff[0][a]/s_suff[1][a])*
	        (s_g_a[a]-s_suff[0][a]/s_suff[1][a]);
	    }
	  logprior_new = G_ZERO;
	  for(i=0;i<3;i++)
	    logprior_new += -(aprior[i]-G_ONE)*log(s_par_new[i])-bprior[i]*s_par_new[i];
	  //d = -G_HALF*tau*(S_new-S)+logprior_new-logprior; 
	  d = -G_HALF*tau*(S_new-S);
	  u = genunf(G_ZERO,G_ONE);
	  if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
	    {
	      for(j=0;j<s_npar;j++)
		i_D_g_a->g_a_par[j] = s_par_new[j];
	      i = i_new;
	      acc_par[p]++;
	      S = S_new;
	      logprior = logprior_new;
	    }
          #ifdef DEBUG_G_A2
	  printf("new:c=%lg,theta=%lg,gamma=%lg,S=%lf\n",
		 s_par_new[0],s_par_new[1],s_par_new[2],S_new);
	  printf("g_a=");
	  for(a=0;a<i_D_g_a->ncat;a++)
	    printf("%lf, ",s_g_a[a]);
	  printf("\n");
        #endif
	}
      //Sampling intercept and slope
      if(1)
	{
	  mu[0] = G_ZERO;
	  mu[1] = G_ZERO;
	  Q[0][0] = G_ZERO;
	  Q[0][1] = G_ZERO;
	  Q[1][0] = G_ZERO;
	  Q[1][1] = G_ZERO;
	  bvec[0] = G_ZERO;
	  bvec[1] = G_ZERO;
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      bvec[0] += sum_by_cat[a] - N_a_beta[a][0] - N_a_beta[a][1]*s_g_a[a];
	      bvec[1] += (sum_by_cat[a] - N_a_beta[a][0] - N_a_beta[a][1]*s_g_a[a]) * s_g_a[a];
	      Q[0][0] += N_a[a];
	      Q[0][1] += N_a[a] * s_g_a[a];
	      Q[1][1] += N_a[a] * s_g_a[a] * s_g_a[a];
	    }
	  Q[0][0] *= tau;
	  Q[0][1] *= tau;
	  Q[1][1] *= tau;
	  Q[1][0] = Q[0][1];
	  bvec[0] *= tau;
	  bvec[1] *= tau;
	  for(i=0;i<2;i++)
	    eps[i] = gennor(G_ZERO,G_ONE);
	  err = choldc0(Q,2);
	  cholsl0(Q,2,bvec,mu);
	  chollTl0(Q,2,eps,eps2);
	  i_length->par->eff[0][0][0][0] = mu[0] + eps2[0];
	  i_length->par->eff[0][1][0][0] = mu[1] + eps2[1];
	}
    }

  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);


  for(a=0;a<i_D_g_a->ncat;a++)
    i_D_g_a->g_a[a] = s_g_a[a];
  if(0)
    {
      for(j=0;j<3;j++)
	i_D_g_a->g_a_par[j] = s_par_opt[j];
    }


  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_g_a: End of sample_g_a\n\n");
  #endif

  #ifdef DEBUG_G_A
  printf("c=%lg,theta=%lg,gamma=%lg,S=%lg\n",
         s_par_opt[0],s_par_opt[1],
         s_par_opt[2],S);
  printf("g_a=\n");
  for(a=0;a<i_D_g_a->ncat;a++)
    printf("%lf %lf\n",s_g_a[a],s_suff[0][a]/s_suff[1][a]);
  #endif

  #ifdef DEBUG_G_A2
  printf("c=%lg,theta=%lg,gamma=%lg,S=%lf acc=(%lf %lf %lf)\n",
         i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],
         i_D_g_a->g_a_par[2],S,
         (double) acc_par[0]/(double) num,
         (double) acc_par[1]/(double) num,
         (double) acc_par[2]/(double) num);
  printf("g_a=");
  for(a=0;a<i_D_g_a->ncat;a++)
    printf("%lf, ",s_g_a[a]);
  printf("\n");
  #endif

  Fmatrix_2d(&Q[0][0],&Q[0]);
  FREE(acc_par);
  FREE(N_a);
  Fmatrix_2d(&N_a_beta[0][0],&N_a_beta[0]);
  Fmatrix_2d(&suff[0][0],&suff[0]);
  FREE(sum_by_cat);
  FREE(sum_beta1_by_cat);


  return(0);
}               /* end of sample_g_a */




/*L:sample_g_a_S_R*
________________________________________________________________

		sample_g_a_S_R
________________________________________________________________

Name:		sample_g_a_S_R
Syntax:		
Description:    Samples the function g(a) describing the non-linear
                link between age and length, assuming g(a)=log(1-exp(-K(a-a_0)))
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int sample_g_a_S_R(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		   Data_g_a *i_D_g_a,double **i_suff)
{
  int       a,a0,aA,i,k,acc_par,acc_par2,num=1000,err;
  double    t,df=3.0;
  double    S, S_a,S_opt,g_opt,d,u,v,tau;
  double    S_new,f1,f2,f3,f_old,f_new,prob,sigma;
  Data_glm   *glm;

  s_suff = i_suff;
  s_ncat = i_D_g_a->ncat;
  s_a_vec = i_D_g_a->a_vec;
  a0 = s_a_vec[0];
  aA = s_a_vec[s_ncat-1];
  glm = i_D_lga->glm;
  tau = i_length->par->tau_obs;

  /* Optimize S */
  for(i=0;i<s_npar;i++)
    s_par_opt[i] = i_D_g_a->g_a_par[i];
  err = find_opt_g_a(s_par_opt,&g_opt);
  if(err)
    {
      write_warning("sample_g_a_S_R:Error calling find_opt_g_a\n");
      return(1);
    }
  S_opt = G_ZERO;
  
  //printf("g_opt=");
  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,s_par_opt,s_g_a);
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      if(s_suff[1][a]>0)
       {
	S_a = s_suff[1][a]*(s_g_a[a]-s_suff[0][a]/s_suff[1][a])*
                           (s_g_a[a]-s_suff[0][a]/s_suff[1][a]);
        S_opt += S_a;
       }
      //printf("%lf, %lf, ",s_g_a[a],S_a);
    } 
  //printf("\n");
  //for(i=0;i<s_npar;i++)
  //   fprintf(stderr,"par[%d]=%lg ",i,exp(s_par_opt[i]));
  //fprintf(stderr,"S_opt=%lf, g_opt=%lf\n",S_opt,g_opt);

  /* Run Metropolis-Hastings for parameters */

  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
  S = G_ZERO;
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      if(s_suff[1][a]>0)
	S += s_suff[1][a]*(s_g_a[a]-s_suff[0][a]/s_suff[1][a])*
                          (s_g_a[a]-s_suff[0][a]/s_suff[1][a]);
    }
  acc_par = 0;
  acc_par2 = 0;
  prob = 0.01;
  sigma = 0.01;
  for(k=0;k<num;k++)
    {
      v = ranf();
      if(v<prob)
	{
          for(i=0;i<s_npar;i++)
	    {
	      t = gennor(G_ZERO,sigma)*df/genchi(df);
	      s_par_new[i] = i_D_g_a->g_a_par[i] + t;
	    }
	}
      else
	{
          for(i=0;i<s_npar;i++)
	    {
	      t = gennor(G_ZERO,sigma)*df/genchi(df);
	      s_par_new[i] = s_par_opt[i] + t;
	    }
	}
      f1 = G_ONE;
      f2 = G_ONE;
      f3 = G_ONE;
      for(i=0;i<s_npar;i++)
	{
          f1 *= dt(s_par_new[i],i_D_g_a->g_a_par[i],sigma,df);
          f2 *= dt(s_par_new[i],s_par_opt[i],sigma,df);
          f3 *= dt(i_D_g_a->g_a_par[i],s_par_opt[i],sigma,df);
	}
      f_new = (prob*f1+(G_ONE-prob)*f2);
      f_old = (prob*f1+(G_ONE-prob)*f3);
      //p_new = s_prior_g_a_par(s_npar,s_par_new);
      //f_new *= p_new;
      //p_old = s_prior_g_a_par(s_npar,i_age->g_a_par);
      //f_old *= p_old;
      S_new = G_ZERO;
      err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,s_par_new,s_g_a);
      for(a=0;a<i_D_g_a->ncat;a++)
	{
          if(s_suff[1][a]>0)
	    S_new += s_suff[1][a]*(s_g_a[a]-s_suff[0][a]/s_suff[1][a])*
                                  (s_g_a[a]-s_suff[0][a]/s_suff[1][a]);
	}
      d = -G_HALF*tau*(S_new-S)+log(f_old)-log(f_new);
      u = genunf(G_ZERO,G_ONE);
      if(d > -1.0e-32 && d < 1.0e32 && log(u) < d)
	{
          for(i=0;i<s_npar;i++)
            i_D_g_a->g_a_par[i] = s_par_new[i];
          acc_par++;
          if(v<prob)
            acc_par2++;
          S = S_new;
	}
    }
  /*
  fprintf(stderr,"theta=%lg, gamma=%lg,S=%lf (%lf,%lf)\n",exp(i_age->g_a_par[0]),
	  exp(i_age->g_a_par[1]),S,(double) acc_par/(double) num,
          (double) acc_par2/(double) num);
  */  

  //printf("g_opt=");
  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
  for(a=0;a<i_D_g_a->ncat;a++)
    {
      i_D_g_a->g_a[a] = s_g_a[a];
      //printf("%lf, ",s_g_a[a]);
      //fprintf(stderr,"%d %lg %lg %lg\n",a,s_g_a[a],s_suff[1][a],s_suff[0][a]/s_suff[1][a]);
    }
  //printf("\n");
  //printf("g_par=(%lf,%lf)\n",i_age->g_a_par[0],i_age->g_a_par[1]);
  //fprintf(stderr,"\n");

  return(0);
}               /* end of sample_g_a_S_R */




/*L:calc_g_a*
________________________________________________________________

		calc_g_a
________________________________________________________________

Name:		calc_g_a
Syntax:		
Description:    Calculates g_a function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
int calc_g_a(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  return(s_calc_g_a(i_ncat,i_a_vec,i_par,o_g));

}		/* end of calc_g_a */



/*L:calc_g_a_S_R*
________________________________________________________________

		calc_g_a_S_R
________________________________________________________________

Name:		calc_g_a_S_R
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_S_R(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  int    i;
  double a,amin,amax,g,g_min,g_max,r;
  double theta,gamma,c;

  amin = i_a_vec[0];
  amax = i_a_vec[i_ncat-1];


  c = i_par[0];
  theta = i_par[1];
  gamma = i_par[2];
  //theta = exp(log(i_par[1])*c);
  //gamma = exp(log(i_par[2])*c);

  if(c < 0.00000001)
    {
      // log(a)
      g_min = log(amin);
      g_max = log(amax);
      r = g_max-g_min;
      for(i=0;i<i_ncat;i++)
	{
	  a = i_a_vec[i];
	  o_g[i] = (log(a)-g_min)/r;
	}
    }
  else
    {
      g_min = -log(theta+exp(-gamma * exp(c * log(amin))));
      g_max = -log(theta+exp(-gamma * exp(c * log(amax))));
      r = g_max-g_min;
      for(i=0;i<i_ncat;i++)
	{
	  //a = ((double) i_a_vec[i]-amin)/(amax-amin);
	  a = i_a_vec[i];
	  g =  -log(theta+exp(-gamma * exp(c * log(a))));
          o_g[i] = (g-g_min)/r;
	}
    }

  return(0);
}		/* end of calc_g_a_S_R */



/*L:calc_g_a_S_R_unorm*
________________________________________________________________

		calc_g_a_S_R_unorm
________________________________________________________________

Name:		calc_g_a_S_R_unorm
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_S_R_unorm(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  int    i;
  double a,g;
  double theta,gamma,c;

  c = i_par[0];
  theta = i_par[1];
  gamma = i_par[2];
  //theta = exp(log(i_par[1])*c);
  //gamma = exp(log(i_par[2])*c);
  if(c < s_parmin[0])
    {
      // log(a)
      for(i=0;i<i_ncat;i++)
	{
	  a = i_a_vec[i];
	  o_g[i] = log(a);
	}
    }
  else
    {
      for(i=0;i<i_ncat;i++)
	{
	  //a = ((double) i_a_vec[i]-amin)/(amax-amin);
	  a = i_a_vec[i];
	  g = (log(theta+exp(-gamma)-
               log(theta+exp(-gamma * exp(c * log(a)))))) *
              (theta+exp(-gamma))/
	      (c*theta*gamma*exp(-gamma));

          o_g[i] = g;;
	}
    }
  

  return(0);
}		/* end of calc_g_a_S_R_unorm */



/*L:prior_g_a_S_R*
________________________________________________________________

		prior_g_a_S_R
________________________________________________________________

Name:		prior_g_a_S_R
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int prior_g_a_S_R(int i_npar,double *i_par)
{
  double p,r,theta,gamma;

  theta = exp(i_par[1]);
  gamma = exp(i_par[0]*theta);
  //prior for theta is assumed exponential
  p = exp(-theta);
  //Prior for gamma is assumed gamma(r,r)
  r = 0.01;
  p *= exp(log(gamma)*(r-G_ONE))*exp(-r*gamma);
  //Transforming to prior for i_par
  p /= (theta*theta*gamma);
  
  return(p);
}		/* end of prior_g_a_S_R */



/*L:calc_g_a_polyn3*
________________________________________________________________

		calc_g_a_polyn3
________________________________________________________________

Name:		calc_g_a_polyn3
Syntax:		
Description:    Calculates polyn3 function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_polyn3(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  int    i;
  double a,a_norm,amin,amax,g,g_min,g_max,r;
  double beta1,beta2,beta3;

  amin = i_a_vec[0];
  amax = i_a_vec[i_ncat-1];
  beta1 = exp(i_par[0]);
  beta3 = i_par[1];
  beta2 = G_ONE+beta1-beta3;
  g_min = exp(G_ZERO);
  g_max = exp(G_ONE);
  r = g_max-g_min;
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      a_norm = (log(a)-log(i_a_vec[0]))/
               (log(i_a_vec[i_ncat-1])-log(i_a_vec[0]));
      g = ((beta3*a_norm+beta2)*a_norm+beta3)*a_norm;
      o_g[i] = g;
    }

  return(0);
}		/* end of calc_g_a_polyn3 */



/*L:calc_g_a_polyn3_unorm*
________________________________________________________________

		calc_g_a_polyn3_unorm
________________________________________________________________

Name:		calc_g_a_polyn3_unorm
Syntax:		
Description:    Calculates polyn3 function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_polyn3_unorm(int i_ncat,double *i_a_vec,
                                 double *i_par,double *o_g)
{
  int    i;
  double a,a_norm,amin,amax,g;
  double beta1,beta2,beta3;

  amin = i_a_vec[0];
  amax = i_a_vec[i_ncat-1];
  beta1 = exp(i_par[0]);
  beta3 = i_par[1];
  beta2 = G_ONE+beta1-beta3;
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      a_norm = (a-i_a_vec[0])/(i_a_vec[i_ncat-1]-i_a_vec[0]);
      g = exp(((-beta1*a_norm+beta2)*a_norm+beta3)*a_norm);
      o_g[i] = g;
    }

  return(0);
}		/* end of calc_g_a_polyn3_unorm */



/*L:calc_S_unorm*
________________________________________________________________

		calc_S_unorm
________________________________________________________________

Name:		calc_S_unorm
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double calc_S_unorm(double *i_par)
{
  int    a,err,j;
  double res,S;
  double beta0,beta1;
  double sum_w,mean_u,mean_g,mean_ug,mean_gg;

  for(j=0;j<s_npar;j++)
    s_par[j] = s_parmin[j]+(s_parmax[j]-s_parmin[j])*
      exp(i_par[j+1])/(G_ONE+exp(i_par[j+1]));

  S = G_ZERO;
  err = s_calc_g_a_unorm(s_ncat,s_a_vec,s_par,s_g_a);
  sum_w =  G_ZERO;
  mean_u = G_ZERO;
  mean_g = G_ZERO;
  mean_ug= G_ZERO;
  mean_gg= G_ZERO;
  for(a=0;a<s_ncat;a++)
    {
      sum_w += s_suff[1][a];
      mean_u  += s_suff[0][a];
      mean_g  += s_suff[1][a]*s_g_a[a];
      mean_ug += s_suff[0][a]*s_g_a[a];
      mean_gg += s_suff[1][a]*s_g_a[a]*s_g_a[a];
    }
  mean_u  /= sum_w;
  mean_g  /= sum_w;
  mean_ug /= sum_w;
  mean_gg /= sum_w;
  beta1 = (mean_ug-mean_u*mean_g)/(mean_gg-mean_g*mean_g);
  beta0 = mean_u-beta1*mean_g;
  for(a=0;a<s_ncat;a++)
    {
      if(s_suff[1][a]>0)
	{
	  res = beta0+beta1*s_g_a[a]-s_suff[0][a]/s_suff[1][a];
	  S += s_suff[1][a]*res*res;
	}
    }
  if(!(S > -1 && S < 99999999999999999.99))
    S = 99999999999999999.99;
  #ifdef DEBUG_G_A2
  printf("par[0]=%lg,par[1]=%lg,par[2]=%lf,S=%lf\n",s_par[0],s_par[1],s_par[2],S);
  #endif
  return(S);
}		/* end of calc_S_unorm */



/*L:calc_S_unorm0*
________________________________________________________________

		calc_S_unorm0
________________________________________________________________

Name:		calc_S_unorm0
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double calc_S_unorm0(double *i_par)
{
  int    a,err,j;
  double res,S;
  double beta0,beta1;
  double sum_w,mean_u,mean_g,mean_ug,mean_gg;

  S = G_ZERO;
  for(j=0;j<s_npar;j++)
    {
      s_par[j] = i_par[1+j];
      if(s_par[j] > 10000.00)
	return(99999999999999999.99);
    }
  err = s_calc_g_a_unorm(s_ncat,s_a_vec,s_par,s_g_a);
  sum_w =  G_ZERO;
  mean_u = G_ZERO;
  mean_g = G_ZERO;
  mean_ug= G_ZERO;
  mean_gg= G_ZERO;
  for(a=0;a<s_ncat;a++)
    {
      sum_w += s_suff[1][a];
      mean_u  += s_suff[0][a];
      mean_g  += s_suff[1][a]*s_g_a[a];
      mean_ug += s_suff[0][a]*s_g_a[a];
      mean_gg += s_suff[1][a]*s_g_a[a]*s_g_a[a];
    }
  mean_u  /= sum_w;
  mean_g  /= sum_w;
  mean_ug /= sum_w;
  mean_gg /= sum_w;
  beta1 = (mean_ug-mean_u*mean_g)/(mean_gg-mean_g*mean_g);
  beta0 = mean_u-beta1*mean_g;
  for(a=0;a<s_ncat;a++)
    {
      if(s_suff[1][a]>0)
	{
	  res = beta0+beta1*s_g_a[a]-s_suff[0][a]/s_suff[1][a];
	  S += s_suff[1][a]*res*res;
	}
    }
  if(!(S > -1 && S < 99999999999999999.99))
    S = 99999999999999999.99;
  #ifdef DEBUG_G_A2
  printf("par[0]=%lg,par[1]=%lg,par[2]=%lf,S=%lf\n",s_par[0],s_par[1],s_par[2],S);
  #endif
  return(S);
}		/* end of calc_S_unorm0 */



/*L:calc_S_beta*
________________________________________________________________

		calc_S_beta
________________________________________________________________

Name:		calc_S_beta
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double calc_S_beta(double *i_par)
{
  int    a,err,j;
  double sum_w,mean_u,mean_g,mean_ug,mean_gg;
  double beta0,beta1,S,res;

  for(j=0;j<s_npar;j++)
    s_par[j] = exp(i_par[1+j]);

  if(s_par[0] < 0.0000001 || s_par[0] > 7.0)
    return(99999999999999999.99);
  for(j=1;j<s_npar;j++)
    {
      if(s_par[j] > 10000.00)
	return(99999999999999999.99);
    }
  err = s_calc_g_a_unorm(s_ncat,s_a_vec,s_par,s_g_a);

  sum_w =  G_ZERO;
  mean_u = G_ZERO;
  mean_g = G_ZERO;
  mean_ug= G_ZERO;
  mean_gg= G_ZERO;
  for(a=0;a<s_ncat;a++)
    {
      sum_w += s_suff[1][a];
      mean_u += s_suff[0][a];
      mean_g += s_suff[1][a]*s_g_a[a];
      mean_ug += s_g_a[a]*s_suff[0][a];
      mean_gg += s_suff[1][a]*s_g_a[a]*s_g_a[a];
    }
  mean_u  /= sum_w;
  mean_g  /= sum_w;
  mean_ug /= sum_w;
  mean_gg /= sum_w;
  if(mean_gg < mean_g*mean_g)
	return(99999999999999999.99);

  beta1 = (mean_ug-mean_u*mean_g)/(mean_gg-mean_g*mean_g);
  beta0 = mean_u-beta1*mean_g;

  S = G_ZERO;

  for(a=0;a<s_ncat;a++)
    {
      if(s_suff[1][a]>0)
	{
	  res = beta0+beta1*s_g_a[a]-s_suff[0][a]/s_suff[1][a];
	  S += s_suff[1][a]*res*res;
	}
    }
  if(!(S > -1 && S < 99999999999999999.99))
    S = 99999999999999999.99;
  #ifdef DEBUG_G_A2
  //printf("beta0=%lg,beta1=%lg\n",beta0,beta1);
  //printf("par[0]=%lg,par[1]=%lg,par[2]=%lf,S=%lf\n",s_par[0],s_par[1],s_par[2],S);
  printf("%12.10lg %12.10lg %12.10lg %12.10lf %12.10lg %lg\n",beta0,beta1,s_par[0],s_par[1],s_par[2],S);
  #endif
  return(S);
}		/* end of calc_S_beta */



/*L:calc_S_beta0*
________________________________________________________________

		calc_S_beta0
________________________________________________________________

Name:		calc_S_beta0
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static double calc_S_beta0(double *i_par)
{
  int    a,err,j;
  double sum_w,mean_u,mean_g,mean_ug,mean_gg;
  double beta0,beta1,S,res;

  for(j=0;j<s_npar;j++)
    {
      s_par[j] = i_par[1+j];
      if(s_par[j] > 10000.00)
	return(99999999999999999.99);
    }
  err = s_calc_g_a_unorm(s_ncat,s_a_vec,s_par,s_g_a);

  sum_w =  G_ZERO;
  mean_u = G_ZERO;
  mean_g = G_ZERO;
  mean_ug= G_ZERO;
  mean_gg= G_ZERO;
  for(a=0;a<s_ncat;a++)
    {
      sum_w += s_suff[1][a];
      mean_u += s_suff[0][a];
      mean_g += s_suff[1][a]*s_g_a[a];
      mean_ug += s_g_a[a]*s_suff[0][a];
      mean_gg += s_suff[1][a]*s_g_a[a]*s_g_a[a];
    }
  mean_u  /= sum_w;
  mean_g  /= sum_w;
  mean_ug /= sum_w;
  mean_gg /= sum_w;
  beta1 = (mean_ug-mean_u*mean_g)/(mean_gg-mean_g*mean_g);
  beta0 = mean_u-beta1*mean_g;

  S = G_ZERO;

  for(a=0;a<s_ncat;a++)
    {
      if(s_suff[1][a]>0)
	{
	  res = beta0+beta1*s_g_a[a]-s_suff[0][a]/s_suff[1][a];
	  S += s_suff[1][a]*res*res;
	}
    }
  if(!(S > -1 && S < 99999999999999999.99))
    S = 99999999999999999.99;
  #ifdef DEBUG_G_A2
  //printf("beta0=%lg,beta1=%lg\n",beta0,beta1);
  //printf("par[0]=%lg,par[1]=%lg,par[2]=%lf,S=%lf\n",s_par[0],s_par[1],s_par[2],S);
  printf("%12.10lg %12.10lg %12.10lg %12.10lf %12.10lg %lg\n",beta0,beta1,s_par[0],s_par[1],s_par[2],S);
  #endif
  return(S);
}		/* end of calc_S_beta0 */



/*L:find_opt_g_a*
________________________________________________________________

		find_opt_g_a
________________________________________________________________

Name:		find_opt_g_a
Syntax:		
Description:    optimize g_a
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int find_opt_g_a(double *o_par_opt,double *o_g)
{
  int    err,i,j,p,iter;
  double e,lambda,ftol,Smin,S0;

  
  p = s_npar;

  lambda = 1.0;
  for(j=1;j<=p;j++)
  for(i=1;i<=p;i++)
    s_xi[j][i] = G_ZERO;
  for(j=1;j<=p;j++)
    s_xi[j][j] = 1.0;

  ftol = 0.000001;

  // First try log model
  s_p[1] = 0.0000000000000001;
  for(i=2;i<=s_npar;i++)
    s_p[i] = G_ONE;
  S0 = calc_S_unorm0(s_p);
  #ifdef DEBUG_G_A
  printf("S0=%lf\n",S0);
  #endif

  // Then non-linear model
  s_p[1] = -3.0;
  for(i=0;i<s_npar;i++)
    {
      e = (3.1-s_parmin[i])/(s_parmax[i]-s_parmin[i]);
      s_p[i+1] = log(e/(G_ONE-e));
    }

  err = powell(s_p,s_xi,p,ftol,&iter,&Smin,calc_S_unorm);
  if(err)
    {
      write_warning("find_opt_g_a:Error calling powell\n");
      return(1);
    }
  if(Smin>S0)
    { // linear model
      o_par_opt[0] = G_ZERO;
     for(i=1;i<s_npar;i++)
       o_par_opt[i] = G_ONE;
     Smin = S0;
    }
  else
    {
      for(i=0;i<s_npar;i++)
      o_par_opt[i] = s_parmin[i]+
	(s_parmax[i]-s_parmin[i])*exp(s_p[i+1])/(G_ONE+exp(s_p[i+1]));
    }

  *o_g = Smin;
  return(0);
}		/* end of find_opt_g_a */
