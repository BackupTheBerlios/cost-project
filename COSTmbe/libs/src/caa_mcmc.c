/*!
  \file caa_main.c
  \brief Containing the main routines for MCMC simulation
  \author Geir Storvik, Hanne Rognebakke

*/

/*Include Files:*/
#include "caa.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif


/*!
 \author Geir Storvik
  \brief Main part for lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Calculating sufficient statistics for the lga model
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating precision parameters
  - Calculating likelihood
*/
int MCMC_it_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h)
{
  int err;

  /* Calculate sufficient statistics */
  err = make_suff_lga(i_D_lga,i_D_g_a,i_start_h);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
      #endif
      write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
      return(err);
    }

  /* Sample effects */
  err = sample_gauss_eff(i_length->gr_str,i_length->par,i_D_lga->glm,i_start_h);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_gauss_eff\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_gauss_eff\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("MCMC_it_lga:\n");
     printf("Int=%lf\n",i_length->par->eff[0][0][0][0]);
     printf("Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Int=%lf\n",i_length->par->eff[0][0][0][0]);
  fprintf(g_caa_log,"Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif
  //i_length->par->eff[0][0][0][0] = 1.154;
  //i_length->par->eff[0][1][0][0] = 3.256;
  /* Sample precision parameters */
  //i_length->par->tau[0][1] = 100.0;
  //i_length->par->tau_obs = 100.0;
  err = sample_precision_lin(i_start_h,i_length->par,i_D_lga->glm);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_precision_lin\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_precision_lin\n");
      return(err);
    } 
  err = sample_precision(i_start_h,i_length->par,i_D_lga->glm);
  if(err)
    {
      write_warning("MCMC_it_lga:Error calling sample_precision\n");
      return(err);
    }
  //i_length->par->tau_obs = 45.0;
  #ifdef DEBUG_PROG
  int i,j,isp;
  Data_cov *xcov;
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      isp = xcov->ispat;
      for(j=0;j<xcov->n_cov;j++)
	if(!xcov->fix[j] && j != isp && j != xcov->icell)
	  printf("tau[%d][%d]=%lf\n",i,j,i_length->par->tau[i][j]);
    }
  printf("tau_obs=%lf\n",i_length->par->tau_obs);
  #endif

  err = calc_lik_lin(i_length,i_D_lga,&(i_length->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_lin\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("Loglik_lga=%lf\n",i_length->par->loglik);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_lga=%lf\n",i_length->par->loglik);
  #endif

  return(0);
}		/* end of MCMC_it_lga */



/*!
  \author Geir Storvik
  \brief Nonlinear connection in the lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Calculate sufficiednt statistics for the non-linear model
  - Simulating parameters in g-function
*/
int MCMC_it_g_a(Data_age *i_D_age,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		LW_struct *i_length,int i_start_h,int i_it)
{
  int err;

  err = suff_g_a(i_length,i_D_age,i_D_lga,i_D_g_a,i_start_h,i_D_g_a->suff);
  if(err)
    {
      write_warning("MCMC_it_g_a:Error calling suff_g_a\n");
      return(err);
    }

  if(i_D_g_a->g_a_model>0)
    {
      #ifdef DEBUG_PROG
        printf("\nSampling g_a-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling g_a-parameters\n");
      #endif
      /* Sample nonlinear function g(a) */
      err = sample_g_a(i_length,i_D_age,i_D_lga,i_D_g_a,i_D_g_a->suff,i_it);
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling g_a-parameters: %d \n",err);
      #endif
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_g_a:Error calling sample_g_a\n");
          #endif
	  write_warning("MCMC_it_g_a:Error calling sample_g_a\n");
	  return(err);
	}
    }
  else if(i_D_g_a->g_a_model!=0)
    {
      printf("MCMC_it_g_a:Unknown age-length relation\n");
      return(err);
    }

  return(0);
}		/* end of MCMC_it_g_a */



/*!
  \author Geir Storvik
  \brief Main part for lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating precision parameters
*/
int MCMC_it_wgl(Data_lin *i_D_wgl,LW_struct *i_weight,int start_h)
{
  int err;

  /* Sample effects */
  err = sample_gauss_eff(i_weight->gr_str,i_weight->par,i_D_wgl->glm,start_h);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_gauss_eff\n");
      return(err);
    }

  /* Sample precision parameters */
	  
  err = sample_precision_lin(start_h,i_weight->par,i_D_wgl->glm);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_precision_lin\n");
      return(err);
    }
  err = sample_precision(start_h,i_weight->par,i_D_wgl->glm);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_precision\n");
      return(err);
    }

  return(0);
}		/* end of MCMC_it_wgl */

