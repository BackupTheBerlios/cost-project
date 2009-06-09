#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "caa.h"

/*!
  \brief Writes the containts of an LW-struct to file "filename" in the current directory.

  Only to be used for testing.
  \author Hanne Rognebakke
*/
int write_LW_struct(LW_struct *i_lin,Data_glm *i_glm,char *i_filename)
{
  int i,j,k,a;
  int err;
  Data_cov   *xcov;
  FILE *fp;
  
  fp = fopen(i_filename,"w");

  fprintf(fp,"LW_struct\n\n");
  fprintf(fp,"Eff_str:\n");  
  fprintf(fp,"ncat=%d\n",i_glm->ncat);
  for(a=0;a<i_glm->ncat;a++)
    {
      fprintf(fp,"nxcov=%d\n",i_glm->nxcov);
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];
	  fprintf(fp,"  n_cov=%d\n",xcov->n_cov);
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      for(k=0;k<xcov->n_fac[j];k++)
		fprintf(fp,"  eff[%d][%d][%d][%d]=%f\n",a,i,j,k,i_lin->par->eff[a][i][j][k]);
	      fprintf(fp,"  ssq[%d][%d][%d]=%f\n",a,i,j,i_lin->par->ssq[a][i][j]);
	      fprintf(fp,"  n_ssq[%d][%d][%d]=%f\n",a,i,j,i_lin->par->n_ssq[a][i][j]);
	    }
	}
    }
  for(i=0;i<i_glm->nxcov;i++)
    {
      fprintf(fp,"ar[%d]=%f\n",i,i_lin->par->ar[i]);
      xcov = i_glm->xcov[i];
      for(j=0;j<2;j++)
	fprintf(fp,"prior_ar[%d][%d]=%f\n",i,j,i_lin->par->prior_ar[i][j]);
      for(j=0;j<xcov->n_cov;j++)
	{
	  for(k=0;k<xcov->n_fac[j];k++)
	    fprintf(fp,"prior_mean[%d][%d][%d]=%f\n",i,j,k,i_lin->par->prior_mean[i][j][k]);
	  for(k=0;k<2;k++)
	    fprintf(fp,"prior_prec[%d][%d][%d]=%f\n",i,j,k,i_lin->par->prior_prec[i][j][k]);
	  fprintf(fp,"tau[%d][%d]=%f\n",i,j,i_lin->par->tau[i][j]);
	}
    }
  for(i=0;i<2;i++)
    fprintf(fp,"prior_prec_obs[%d]=%f\n",i,i_lin->par->prior_prec_obs[i]);

  fprintf(fp,"tau_obs=%f\n",i_lin->par->tau_obs);
  fprintf(fp,"loglik=%f\n",i_lin->par->loglik);
  fprintf(fp,"num_var=%d\n",i_lin->par->num_var);

  //fprintf(fp,"cell[0][0]=%f\n",i_lin->par->cell[0][0]);
  fprintf(fp,"mcmc[0]=%f\n",i_lin->par->mcmc[0]);
  
  fprintf(fp,"\nGraph_str:\n");  
  for(i=0;i<i_glm->nxcov;i++)  
    {
      for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	fprintf(fp,"in_gr[%d][%d]=%d\n",i,j,i_lin->gr_str->in_gr[i][j]);
    }
  err = GMRFLib_print_graph(fp,i_lin->gr_str->graph);


  fprintf(fp,"\ncens_model=%d\n",i_lin->cens_model);
  if(i_lin->cens_model)
    {
      fprintf(fp,"cens_k=%f,cens_m=%f,cens_r=%f,cens_Nlim=%f\n",
	      i_lin->cens_k,i_lin->cens_m,i_lin->cens_r,i_lin->cens_Nlim);
    }
  fprintf(fp,"\nfixed_model=%d\n",i_lin->fixed_model);

  fclose(fp);

  return(0);
}           /* end of write_LW_struct */

/*!
  \brief Writes the containts of Data_lin to file "filename" in the current directory.

  Only to be used for testing.
  \author Hanne Rognebakke
*/
int write_Data_lin(Data_lin *i_D_lga,int i_ncat,char *i_filename)
{
  int h,a;
  FILE *fp;
  fp=fopen(i_filename,"w");

  fprintf(fp,"Data_lin\n\n");

  fprintf(fp,"Data_glm:\n");
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      if(i_D_lga->glm->nxcov==3)
	fprintf(fp,"beta_hat[%d][0][0]=%f,beta_hat[%d][0][1]=%f,beta_hat[%d][0][2]=%f\n",
		h,i_D_lga->glm->beta_hat[h][0][0],h,i_D_lga->glm->beta_hat[h][0][1],
		h,i_D_lga->glm->beta_hat[h][0][2]);
      else
	fprintf(fp,"beta_hat[%d][0][0]=%f,beta_hat[%d][0][1]=%f\n",
		h,i_D_lga->glm->beta_hat[h][0][0],h,i_D_lga->glm->beta_hat[h][0][1]);
    }
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      if(i_D_lga->glm->nxcov==3)
	{
	  fprintf(fp,"suff[%d][0][0]=%f,suff[%d][0][1]=%f,suff[%d][0][2]=%f\n",
		  h,i_D_lga->glm->suff[h][0][0],h,i_D_lga->glm->suff[h][0][1],
		  h,i_D_lga->glm->suff[h][0][2]);
	  fprintf(fp,"suff[%d][1][0]=%f,suff[%d][1][1]=%f,suff[%d][1][2]=%f\n",
		  h,i_D_lga->glm->suff[h][1][0],h,i_D_lga->glm->suff[h][1][1],
		  h,i_D_lga->glm->suff[h][1][2]);
	  fprintf(fp,"suff[%d][2][0]=%f,suff[%d][2][1]=%f,suff[%d][2][2]=%f\n",
		  h,i_D_lga->glm->suff[h][2][0],h,i_D_lga->glm->suff[h][2][1],
		  h,i_D_lga->glm->suff[h][2][2]);
	}
      else
	{
	  fprintf(fp,"suff[%d][0][0]=%f,suff[%d][0][1]=%f\n",
		  h,i_D_lga->glm->suff[h][0][0],h,i_D_lga->glm->suff[h][0][1]);
	  fprintf(fp,"suff[%d][1][0]=%f,suff[%d][1][1]=%f\n",
		  h,i_D_lga->glm->suff[h][1][0],h,i_D_lga->glm->suff[h][1][1]);
	}
      fprintf(fp,"ssq[%d]=%f\n",h,i_D_lga->glm->ssq[h]);
    }

  fprintf(fp,"\n");
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<i_ncat;a++)
	{
	  fprintf(fp,"Ages[%d][%d]=%d,sum_by_cat[%d][%d]=%f,sqsum_by_cat[%d][%d]=%f\n",
		  h,a,i_D_lga->Ages[h][a],h,a,i_D_lga->sum_by_cat[h][a],
		  h,a,i_D_lga->sqsum_by_cat[h][a]);
	}
    }
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    fprintf(fp,"haulweight[%d]=%f\n",i_D_lga->haulweight[h]);
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<i_ncat;a++)
	{
	  fprintf(fp,"Ages_fix[%d][%d]=%d,sum_by_cat_fix[%d][%d]=%f,sqsum_by_cat_fix[%d][%d]=%f\n",
		  h,a,i_D_lga->Ages_fix[h][a],h,a,i_D_lga->sum_by_cat_fix[h][a],
		  h,a,i_D_lga->sqsum_by_cat_fix[h][a]);
	}
    }

  fclose(fp);

  return(0);
}            /* end of write_Data_lin */
