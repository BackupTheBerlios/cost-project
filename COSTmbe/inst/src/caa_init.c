/*!
  \file caa_innit.c
  \brief Routines for initialization of variables and structs
  \author Geir Storvik

*/
#include "caa.h"


/*!
  \author Geir Storvik
  \brief Allocate space and initialize age struct and struct for age-data
*/
int init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_coastal_cod,
	     double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
	     double *i_pri_age_prec_par,double *i_pri_age_ar,
             Age_struct **o_age,Age_struct **o_age_mean)
{
  int         a,a2,A,err,h,i,j,ind,N,ind1,ind2,ind3,ind4,ncov,k;
  Age_struct *age;

  err = alloc_age(i_D_age,o_age);
  if(err)
    {
      write_warning("init_age:Error calling alloc_age\n");
      return(err);
    }
  age = *o_age;

  /* Include all effects except haul effect in linear graph */
  for(i=0;i<i_D_age->glm->nxcov;i++)  
    {
      for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
	age->gr_str_f->in_gr[i][j] = 1;
    }

  for(i=0;i<i_D_age->glm->nxcov;i++)  /* Only fixed effects are included in linear graph */
  for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
    age->gr_str_r->in_gr[i][j] = 1-i_D_age->glm->xcov[i]->fix[j];

  A = i_D_age->glm->ncat-1;

  /* Specify that one observation is available on alpha */
  /* Initial values age */
  ind1 = 0;
  ind2 = 0;
  ind3 = 0;
  ind4 = 0;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i==0)
	ncov = i_D_age->glm->xcov[i]->n_cov-1;
      else
	ncov = i_D_age->glm->xcov[i]->n_cov;
      for(j=0;j<ncov;j++)  //Not haul-effect here
	{
	  if(i_D_age->glm->xcov[i]->fix[j])
	    {
	      age->par->tau[i][j] = i_pri_age_eff_prec[ind1];
	      ind1++;
	      for(k=0;k<i_D_age->glm->xcov[i]->n_fac[j];k++)
		{
		  age->par->prior_mean[i][j][k] = i_pri_age_eff_mean[ind2];
		  ind2++;
		}
	    }
	  else
	    {
	      age->par->prior_prec[i][j][0] = i_pri_age_prec_par[ind3];
	      age->par->prior_prec[i][j][1] = i_pri_age_prec_par[ind3+1];
	      ind3 += 2;
	      age->par->tau[i][j] = age->par->prior_prec[i][j][0]/age->par->prior_prec[i][j][1];
	      for(k=0;k<i_D_age->glm->xcov[i]->n_fac[j];k++)
		age->par->prior_mean[i][j][k] = G_ZERO;
	    }
	}
      if(i_D_age->glm->xcov[i]->ispat >= 0)
	{
	  age->par->prior_ar[i][0] = i_pri_age_ar[ind4];
	  age->par->prior_ar[i][1] = i_pri_age_ar[ind4+1];
	  ind4 +=2;
	  age->par->ar[i] = age->par->prior_ar[i][0]/(age->par->prior_ar[i][0]+age->par->prior_ar[i][1]);
	}
    }
  age->par->prior_prec_obs[0] = i_pri_age_prec_par[ind3];
  age->par->prior_prec_obs[1] = i_pri_age_prec_par[ind3+1];
  age->par->tau_obs = age->par->prior_prec_obs[0]/age->par->prior_prec_obs[1];
  age->par->loglik = G_ZERO;

  i_D_age->glm->suff = Mmatrix_3d(0,i_D_age->glm->nHaul-1,   // Free ok
                                  0,i_D_age->glm->nxcov-1,
                                  0,i_D_age->glm->nxcov-1,sizeof(double),1);

  i_D_age->glm->beta_hat = Mmatrix_3d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat-1,0,i_D_age->glm->nxcov-1,
				      sizeof(double),1); // Free ok

  i_D_age->glm->ssq = CALLOC(i_D_age->glm->nHaul,double);  // Free ok

  /* Specify that one observation is available on alpha */
  for(h=0;h<i_D_age->glm->nHaul;h++)  
      i_D_age->glm->suff[h][0][0] = G_ONE;


  for(h=0;h<i_D_age->glm->nHaul;h++)  
      i_D_age->glm->ssq[h] = G_ZERO;

  err = alloc_age(i_D_age,o_age_mean);
  if(err)
    {
      write_warning("init_age:Error calling alloc_age\n");
      return(err);
    }

  int ncat;
  i_D_age->type_age = CALLOC(i_D_age->glm->nHaul,int);   // Free ok
  if(i_age_errors)
    {
      age->A2A = Mmatrix_2d(0,i_D_age->glm->ncat-1,      
                            0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok
      ind = 0;
      if(i_coastal_cod)//Coastal cod
	ncat = i_D_age->glm->ncat/2;
      else
	ncat = i_D_age->glm->ncat;
      for(a=0;a<ncat;a++)
	{
	for(a2=0;a2<ncat;a2++)
	  {
	    age->A2A[a][a2] = i_A2A[ind];
	    if(i_coastal_cod)//Coastal cod
	      age->A2A[ncat+a][ncat+a2] = i_A2A[ind];
	    ind++;
	  }
	}
      /* Find "neighbor" ages */
      age->A_Nneigh = CALLOC(i_D_age->glm->ncat,int);   // Free ok
      age->A_neigh = CALLOC(i_D_age->glm->ncat,int *);  // Free ok
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
          /* First count number of neighbors */
	  N = 0;
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    {
	      if(age->A2A[a2][a]> 0.000000001)
		N++;
	    }
	  age->A_Nneigh[a] = N;
	  age->A_neigh[a] = CALLOC(N,int);              // Free ok
          N = 0;
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    {
	      if(age->A2A[a2][a]> 0.000000001)
		{
		  age->A_neigh[a][N] = a2;
		  N++;
		}
	    }
	}

      for(h=0;h<i_D_age->glm->nHaul;h++)
	i_D_age->type_age[h] = 1;
    }
  else
    {
      for(h=0;h<i_D_age->glm->nHaul;h++)
	i_D_age->type_age[h] = 0;
    }


  age = *o_age_mean;
  age->par->loglik = G_ZERO;



  return(0);
}		/* end of init_age */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated by init_age
*/
int re_init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,
             Age_struct **o_age,Age_struct **o_age_mean)
{
  int         a,err;
  Age_struct *age;

  err = re_alloc_age(i_D_age,o_age);
  if(err)
    {
      write_warning("re_init_age:Error calling re_alloc_age\n");
      return(err);
    }
  age = *o_age;

  Fmatrix_3d(&i_D_age->glm->suff[0][0][0],&i_D_age->glm->suff[0][0],
             &i_D_age->glm->suff[0]);

  Fmatrix_3d(&i_D_age->glm->beta_hat[0][0][0],&i_D_age->glm->beta_hat[0][0],
             &i_D_age->glm->beta_hat[0]);
  FREE(i_D_age->glm->ssq);

  err = re_alloc_age(i_D_age,o_age_mean);
  if(err)
    {
      write_warning("re_init_age:Error calling re_alloc_age\n");
      return(err);
    }

  FREE(i_D_age->type_age);
  if(i_age_errors)
    {
      Fmatrix_2d(&age->A2A[0][0],&age->A2A[0]);
      /* Find "neighbor" ages */
      FREE(age->A_Nneigh);
      for(a=0;a<i_D_age->glm->ncat;a++)
	  FREE(age->A_neigh[a]);
      FREE(age->A_neigh);

    }

  return(0);
}		/* end of re_init_age */


/*!
  \author Geir Storvik
  \brief Allocating space for age-struct
*/
int alloc_age(Data_age *i_D_age,Age_struct **o_age)
{
  int         err,i,j;
  Age_struct *age;

  age = CALLOC(1,Age_struct);                  // Free ok
  age->gr_str_f = CALLOC(1,Graph_str);         // Free ok
  age->gr_str_f->in_gr = CALLOC(i_D_age->glm->nxcov,int *);  //Free ok

  for(i=0;i<i_D_age->glm->nxcov;i++)  
      age->gr_str_f->in_gr[i] = 
	CALLOC(i_D_age->glm->xcov[i]->n_cov,int);  // Free ok

  age->gr_str_r = CALLOC(1,Graph2_str);            // Free ok
  age->gr_str_r->in_gr = CALLOC(i_D_age->glm->nxcov,int *); // Free ok
  for(i=0;i<i_D_age->glm->nxcov;i++)  /* Only fixed effects are included in linear graph */
    {
      age->gr_str_r->in_gr[i] = 
	CALLOC(i_D_age->glm->xcov[i]->n_cov,int);   // Free ok
      for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
	age->gr_str_r->in_gr[i][j] = 1-i_D_age->glm->xcov[i]->fix[j];
    }
  age->alpha = Mmatrix_2d(0,i_D_age->glm->nHaul-1,  // Free ok
		   	    0,i_D_age->glm->ncat-1,sizeof(double),1);

  age->mu = CALLOC(i_D_age->glm->nHaul,double);       // Free ok
  err = alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par);
  if(err)
    {
      write_warning("alloc_age:Error calling alloc_Eff_str\n");
      return(err);
    }
  err = alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par_new);
  if(err)
    {
      write_warning("alloc_age:Error calling alloc_Eff_str\n");
      return(err);
    }

  *o_age = age;
  return(0);
}		/* end of alloc_age */


/*!
  \author Geir Storvik
  \brief Re-allocating space allocated by alloc_age
*/
int re_alloc_age(Data_age *i_D_age,Age_struct **o_age)
{
  int         err,i;
  Age_struct *age;

  age = *o_age;

  for(i=0;i<i_D_age->glm->nxcov;i++)  
      FREE(age->gr_str_f->in_gr[i]);
  FREE(age->gr_str_f->in_gr);
  FREE(age->gr_str_f);


  for(i=0;i<i_D_age->glm->nxcov;i++)  
      FREE(age->gr_str_r->in_gr[i]);
  FREE(age->gr_str_r->in_gr);
  FREE(age->gr_str_r);

  Fmatrix_2d(&age->alpha[0][0],&age->alpha[0]);

  FREE(age->mu);
  err = re_alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par);
  if(err)
    {
      write_warning("re_alloc_age:Error calling re_alloc_Eff_str\n");
      return(err);
    }
  err = re_alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par_new);
  if(err)
    {
      write_warning("re_alloc_age:Error calling re_alloc_Eff_str\n");
      return(err);
    }

  FREE(age);

  return(0);
}		/* end of re_alloc_age */


/*!
  \author Geir Storvik
  \brief Allocate space and initialize linear struct and data for a linear model
*/
int init_lin(Data_lin *i_D_lin,
	     double *i_pri_lin_eff_mean,double *i_pri_lin_eff_prec,
	     double *i_pri_lin_prec_par,double *i_pri_lin_ar,
	     LW_struct **o_lin,LW_struct **o_lin_mean)
{
  int        err,i,j,ind1,ind2,ind3,ind4,k;
  LW_struct *lin;

  err = alloc_lin(i_D_lin,o_lin);
  if(err)
    {
      write_warning("init_lin:Error calling alloc_lin\n");
      return(err);
    }
  /* Initial values */
  ind1 = 0;
  ind2 = 0;
  ind3 = 0;
  ind4 = 0;
  lin = *o_lin;
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      for(j=0;j<i_D_lin->glm->xcov[i]->n_cov;j++)
	{
	  if(i_D_lin->glm->xcov[i]->fix[j])
	    {
	      lin->par->tau[i][j] = i_pri_lin_eff_prec[ind1];
	      ind1++;
	      for(k=0;k<i_D_lin->glm->xcov[i]->n_fac[j];k++)
		{
		  lin->par->prior_mean[i][j][k] = i_pri_lin_eff_mean[ind2];
		  ind2++;
		}
	    }
	  else
	    {
	      lin->par->prior_prec[i][j][0] = i_pri_lin_prec_par[ind3];
	      lin->par->prior_prec[i][j][1] = i_pri_lin_prec_par[ind3+1];
	      ind3 += 2;
	      lin->par->tau[i][j] = lin->par->prior_prec[i][j][0]/lin->par->prior_prec[i][j][1];
	      for(k=0;k<i_D_lin->glm->xcov[i]->n_fac[j];k++)
		lin->par->prior_mean[i][j][k] = G_ZERO;

	      for(k=0;k<i_D_lin->glm->xcov[i]->n_fac[j];k++)
		lin->par->prior_mean[i][j][k] = G_ZERO;
	    }
	}
      if(i_D_lin->glm->xcov[i]->ispat >= 0)
	{
	  lin->par->prior_ar[i][0] = i_pri_lin_ar[ind4];
	  lin->par->prior_ar[i][1] = i_pri_lin_ar[ind4+1];
	  ind4 +=2;
	  lin->par->ar[i] = lin->par->prior_ar[i][0]/
	    (lin->par->prior_ar[i][0]+lin->par->prior_ar[i][1]);
	}
    }
  lin->par->prior_prec_obs[0] = i_pri_lin_prec_par[ind3];
  lin->par->prior_prec_obs[1] = i_pri_lin_prec_par[ind3+1];
  lin->par->tau_obs = lin->par->prior_prec_obs[0]/lin->par->prior_prec_obs[1];
  lin->par->loglik   = G_ZERO;


  err = alloc_lin(i_D_lin,o_lin_mean);
  if(err)
    {
      write_warning("init_lin:Error calling alloc_lin\n");
      return(err);
    }
  lin = *o_lin_mean;

  lin->par->loglik   = G_ZERO;


  return(0);
}		/* end of init_lin */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated by init_lin
*/
int re_init_lin(Data_lin *i_D_lin,
             LW_struct **o_lin,LW_struct **o_lin_mean)
{
  int        err;

  err = re_alloc_lin(i_D_lin,o_lin);
  if(err)
    {
      write_warning("re_init_lin:Error calling re_alloc_lin\n");
      return(err);
    }

  err = re_alloc_lin(i_D_lin,o_lin_mean);
  if(err)
    {
      write_warning("re_init_lin:Error calling alloc_lin\n");
      return(err);
    }

  return(0);
}		/* end of re_init_lin */


/*!
  \author Geir Storvik
  \brief Allocating space for linear struct 
*/
int alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin)
{
  int        err,i,j;
  LW_struct *lin;

  lin = CALLOC(1,LW_struct);             // Free ok
  lin->gr_str = CALLOC(1,Graph_str);     // Free ok
  err = alloc_Eff_str(1,i_D_lin->glm->nxcov,i_D_lin->glm->xcov,&lin->par);
  if(err)
    {
      write_warning("alloc_lin:Error calling alloc_Eff_str\n");
      return(err);
    }
  
  lin->gr_str->in_gr = CALLOC(i_D_lin->glm->nxcov,int *);    // Free ok
  for(i=0;i<i_D_lin->glm->nxcov;i++)  
    {
      lin->gr_str->in_gr[i] = CALLOC(i_D_lin->glm->xcov[i]->n_cov,int);   // Free ok
     for(j=0;j<i_D_lin->glm->xcov[i]->n_cov;j++)
      lin->gr_str->in_gr[i][j] = 1;/* All effects are included in linear graph */
    }

  *o_lin = lin;

  return(0);     /* end of alloc_lin */
}



/*!
  \author Geir Storvik
  \brief Allocating space for linear struct 
*/
int re_alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin)
{
  int        err,i;
  LW_struct *lin;

  lin = *o_lin;

  err = re_alloc_Eff_str(1,i_D_lin->glm->nxcov,i_D_lin->glm->xcov,&lin->par);
  if(err)
    {
      write_warning("re_alloc_lin:Error calling alloc_Eff_str\n");
      return(err);
    }
  
  for(i=0;i<i_D_lin->glm->nxcov;i++)  
     FREE(lin->gr_str->in_gr[i]);
  FREE(lin->gr_str->in_gr);

  FREE(lin->gr_str);
  FREE(lin);

  return(0);     /* end of re_alloc_lin */
}


/*!
  \author Geir Storvik
  \brief Initialize parameters in lga model
*/
int init_lga_par(LW_struct *i_length,Data_lin *i_D_lga,
		 double *i_lga_fixed_int,double *i_lga_fixed_slp,double *i_lga_fixed_tau)
{
  int     h;
  double  Int,Slp,W;

  Int = G_ZERO;
  Slp = G_ZERO;
  W = G_ZERO;

  if(i_length->fixed_model == 0)
    {
      for(h=0;h<i_D_lga->glm->nHaul;h++)
	{
	  Int += i_D_lga->glm->beta_hat[h][0][0]*i_D_lga->glm->suff[h][0][0];
	  Slp += i_D_lga->glm->beta_hat[h][0][1]*i_D_lga->glm->suff[h][0][0];
	  W += i_D_lga->glm->suff[h][0][0];
	}
      Int /= W;
      Slp /= W;
    }
  else
    {
      Int = i_lga_fixed_int[0];
      Slp = i_lga_fixed_slp[0];
      i_length->par->tau_obs = i_lga_fixed_tau[0];
    }

  /* Assume first covariate is constant term */
  if(i_D_lga->glm->xcov[0]->n_fac[0]==1)
      i_length->par->eff[0][0][0][0] = Int;
  else
    {
      write_warning("init_lga_par:First covariate should be constant term\n");
      return(1);
    }

  if(i_D_lga->glm->xcov[1]->n_fac[0]==1)
      i_length->par->eff[0][1][0][0] = Slp;
  else
    {
      write_warning("init_lga_par:First covariate should be constant term\n");
      return(1);
    }
  #ifdef DEBUG_PROG
  printf("Init_lga_par:\n");
  printf("Int = %f\n",i_length->par->eff[0][0][0][0]);
  printf("Slp = %f\n",i_length->par->eff[0][1][0][0]);
  #endif


  return(0);    /* End of init_lga_par */
}



/*F:init_glm_sim*

________________________________________________________________

		init_glm_sim
________________________________________________________________

Name:		init_glm_sim
Syntax:		
Description: Initialize age struct
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_init.c,v 1.1 2009/06/09 10:30:47 mmerzere Exp $
________________________________________________________________
*/
int init_glm_sim(Data_glm *i_glm,Data_glm **o_glm_sim)
{
  int         h;
  Data_glm   *glm_sim;

  glm_sim = CALLOC(1,Data_glm);             // Free ok
  glm_sim->suff = Mmatrix_3d(0,i_glm->nHaul-1,0,1,0,1,   // Free ok
			 sizeof(double),1);
  for(h=0;h<i_glm->nHaul;h++)
    glm_sim->suff[h][0][0] = i_glm->suff[h][0][0];
  glm_sim->nindex = i_glm->nindex;
  glm_sim->index = i_glm->index;
  glm_sim->nxcov = i_glm->nxcov;
  glm_sim->xcov = i_glm->xcov;
  glm_sim->nHaul =i_glm->nHaul;
  glm_sim->beta_hat = Mmatrix_3d(0,i_glm->nHaul-1,0,0,0,1,sizeof(double),1);    // Free ok
  glm_sim->ssq = CALLOC(i_glm->nHaul,double);           // Free ok

  *o_glm_sim = glm_sim;

  return(0);
}		/* end of init_glm_sim */



/*F:re_init_glm_sim*

________________________________________________________________

		re_init_glm_sim
________________________________________________________________

Name:		re_init_glm_sim
Syntax:		
Description: Initialize age struct
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_init.c,v 1.1 2009/06/09 10:30:47 mmerzere Exp $
________________________________________________________________
*/
int re_init_glm_sim(Data_glm *i_glm,Data_glm **o_glm_sim)
{
  Data_glm   *glm_sim;

  *o_glm_sim = glm_sim;
  Fmatrix_3d(&glm_sim->suff[0][0][0],&glm_sim->suff[0][0],&glm_sim->suff[0]);

  Fmatrix_3d(&glm_sim->beta_hat[0][0][0],&glm_sim->beta_hat[0][0],&glm_sim->beta_hat[0]);
  FREE(glm_sim->ssq);

  FREE(glm_sim);

  return(0);
}		/* end of re_init_glm_sim */

