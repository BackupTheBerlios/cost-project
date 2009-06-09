/*!
  \file caa_routines.c
  \brief Routines for converting data into approperiate structures and some other general
  routines
*/
#include "caa.h"

int make_cell_constr_age(Data_age *i_D_age,int *i_icell,
			 double *i_age_int_Sigma_cell,double *i_age_int_constr_cell,
			 int i_age_int_nconstr_cell,
			 double *i_age_hsz_Sigma_cell,double *i_age_hsz_constr_cell,
			 int i_age_hsz_nconstr_cell)
{
  int         i,j,k,ind,ncell,err;
  double     *constr_cell, *Sigma_cell, **W, det;
  Data_cov   *xcov;

  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
#ifdef DEBUG_CELL
      fprintf(stderr,"age icell=%d,ncat=%d\n",i_icell[i]-1,i_D_age->glm->ncat);
#endif
      xcov = i_D_age->glm->xcov[i];
      xcov->icell = i_icell[i]-1;
      if(xcov->icell>=0)
	{
	  ncell = xcov->n_fac[xcov->icell]*i_D_age->glm->ncat;
	  W = CALLOC2_d(ncell,ncell);
	  xcov->Q_cell = CALLOC2_d(ncell,ncell);
	  xcov->Qchol_cell = CALLOC2_d(ncell,ncell);
	  xcov->cell_vec = CALLOC(ncell,double);
	  if(i==0)
	    Sigma_cell = i_age_int_Sigma_cell;
	  else
	    Sigma_cell = i_age_hsz_Sigma_cell;
	  ind=0;
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  W[j][k] = Sigma_cell[ind];
		  ind++;
		}
	    }
	  err = cholinv0(W,ncell,xcov->Q_cell,&det);
	  for(j=0;j<ncell;j++)
	  for(k=0;k<ncell;k++)
	    xcov->Qchol_cell[j][k] = xcov->Q_cell[j][k];
	  err = choldc0(xcov->Qchol_cell,ncell);
	  FREE2_d(W,ncell);

#ifdef DEBUG_CELL
	  FILE *unit;
	  unit = fopen("Q_cell_age.dat","w");
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%lf ",xcov->Q_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
#endif
	  if(i==0)
	    {
	      xcov->n_constr_cell = i_age_int_nconstr_cell;
	      constr_cell = i_age_int_constr_cell;
	    }
	  else
	    {
	      xcov->n_constr_cell = i_age_hsz_nconstr_cell;
	      constr_cell = i_age_hsz_constr_cell;
	    }
	  xcov->constr_cell = CALLOC2_d(xcov->n_constr_cell,ncell);
	  ind = 0;
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  xcov->constr_cell[j][k] = constr_cell[ind];
		  ind++;
		}
	    }
#ifdef DEBUG_CELL
	  unit = fopen("constr_cell_age.dat","w");
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%lf ",xcov->constr_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
#endif
	}
    }
  return(0);
}

int make_cell_constr_lin(Data_lin *i_D_lin,int *i_icell,
			 double *i_lin_int_Sigma_cell,double *i_lin_int_constr_cell,
			 int i_lin_int_nconstr_cell,
			 double *i_lin_slp_Sigma_cell,double *i_lin_slp_constr_cell,
			 int i_lin_slp_nconstr_cell,
			 double *i_lin_hsz_Sigma_cell,double *i_lin_hsz_constr_cell,
			 int i_lin_hsz_nconstr_cell)
{
  int         i,j,k,ind,ncell,err;
  double     *constr_cell, det;
  double     *Sigma_cell, **W;
  Data_cov   *xcov;

  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      xcov = i_D_lin->glm->xcov[i];
      xcov->icell = i_icell[i]-1;
      if(xcov->icell>=0)
	{
	  ncell = xcov->n_fac[xcov->icell];
#ifdef DEBUG_CELL
	  fprintf(stderr,"lin icell=%d,ncell=%d,nconstr=%d\n",
		  i_icell[i]-1,ncell,i_lin_int_nconstr_cell);
#endif
	  W = CALLOC2_d(ncell,ncell);
	  xcov->Q_cell = CALLOC2_d(ncell,ncell);
	  xcov->Qchol_cell = CALLOC2_d(ncell,ncell);
	  xcov->cell_vec = CALLOC(ncell,double);
	  if(i==0)
	    Sigma_cell = i_lin_int_Sigma_cell;
	  else if(i==1)
	    Sigma_cell = i_lin_slp_Sigma_cell;
	  else
	    Sigma_cell = i_lin_hsz_Sigma_cell;
	  ind=0;
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  W[j][k] = Sigma_cell[ind];
		  ind++;
		}
	    }
	  err = cholinv0(W,ncell,xcov->Q_cell,&det);
	  for(j=0;j<ncell;j++)
	  for(k=0;k<ncell;k++)
	    xcov->Qchol_cell[j][k] = xcov->Q_cell[j][k];
	  err = choldc0(xcov->Qchol_cell,ncell);
	  FREE2_d(W,ncell);

#ifdef DEBUG_CELL
	  FILE *unit;
	  unit = fopen("Q_cell_lin.dat","w");
	  for(j=0;j<ncell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%f ",xcov->Q_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
#endif
	  if(i==0)
	    {
	      xcov->n_constr_cell = i_lin_int_nconstr_cell;
	      constr_cell = i_lin_int_constr_cell;
	    }
	  else if(i==1)
	    {
	      xcov->n_constr_cell = i_lin_slp_nconstr_cell;
	      constr_cell = i_lin_slp_constr_cell;
	    }
	  else
	    {
	      xcov->n_constr_cell = i_lin_hsz_nconstr_cell;
	      constr_cell = i_lin_hsz_constr_cell;
	    }
	  xcov->constr_cell = CALLOC2_d(xcov->n_constr_cell,ncell);
	  ind = 0;
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		{
		  xcov->constr_cell[j][k] = constr_cell[ind];
		  ind++;
		}
	    }
#ifdef DEBUG_CELL
	  unit = fopen("constr_cell_lin.dat","w");
	  for(j=0;j<xcov->n_constr_cell;j++)
	    {
	      for(k=0;k<ncell;k++)
		fprintf(unit,"%lf ",xcov->constr_cell[j][k]);
	      fprintf(unit,"\n");
	    }
	  fclose(unit);
#endif
	}
    }
  return(0);
}

/*!
  \author Geir Storvik
  \brief  Put cell distributions into structs
*/  
int makedata_cell_dist(SEXP i_dist_cell,int *i_icell,
		       Age_struct *i_age,Data_age *i_D_age,
		       LW_struct *i_length,Data_lin *i_D_lga,
		       LW_struct *i_weight,Data_lin *i_D_wgl)
{
  int   i,j,k,ind;
  int *num_cell_o, *num_cell_u;
  double *E=NULL,*C=NULL;
  Data_cov *xcov;

  num_cell_o =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.o")));
  num_cell_u =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.u")));

  i_age->par->cell = CALLOC(i_D_age->glm->nxcov,double *);
  i_length->par->cell = CALLOC(i_D_lga->glm->nxcov,double *);
  i_weight->par->cell = CALLOC(i_D_wgl->glm->nxcov,double *);
  ind = 0;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      xcov = i_D_age->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_age->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"age.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"age.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.C"));
		}
	  
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }
  ind = 2;

  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_length->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }
  ind = 5;
  
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_weight->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }
  return(0);
}

int re_makedata_cell_dist(SEXP i_dist_cell,int *i_icell,
			  Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl)
{
  int    i;
  Data_cov *xcov;

  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      xcov = i_D_age->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_age->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_age->par->cell);

  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_length->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_length->par->cell);

  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_weight->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_weight->par->cell);

  return(0);

}

/*!
  \author Geir Storvik
  \brief  Put cell distributions into structs
*/  
int makedata_cell_dist_CC(SEXP i_dist_cell,int *i_icell,
			  Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl,
			  LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			  LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC)
{
  int   i,j,k,ind;
  int *num_cell_o, *num_cell_u;
  double *E=NULL,*C=NULL;
  Data_cov *xcov;

  num_cell_o =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.o")));
  num_cell_u =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.u")));

  i_age->par->cell = CALLOC(i_D_age->glm->nxcov,double *);
  i_length->par->cell = CALLOC(i_D_lga->glm->nxcov,double *);
  i_weight->par->cell = CALLOC(i_D_wgl->glm->nxcov,double *);
  i_length_CC->par->cell = CALLOC(i_D_lga_CC->glm->nxcov,double *);
  i_weight_CC->par->cell = CALLOC(i_D_wgl_CC->glm->nxcov,double *);
  ind = 0;
  /* Age */
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      xcov = i_D_age->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_age->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"age.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"age.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.C"));
		}
	  
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }

  ind = 2;
  /* lga skrei */
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_length->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }

  ind = 5;
  /* wgl skrei */
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_weight->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }

  ind = 8;
  /* lga Costal cod */  
  for(i=0;i<i_D_lga_CC->glm->nxcov;i++)
    {
      xcov = i_D_lga_CC->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_length_CC->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"lga.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }

  ind = 11;
  /* wgl Coastal cod */
  for(i=0;i<i_D_wgl_CC->glm->nxcov;i++)
    {
      xcov = i_D_wgl_CC->glm->xcov[i];
      xcov->icell = i_icell[ind]-1;
      if(i_icell[ind]>=0)
	{
	  xcov->cell_dist = CALLOC(1,Cell_dist);
	  xcov->cell_dist->n_o = num_cell_o[ind];
	  xcov->cell_dist->n_u = num_cell_u[ind];
	  i_weight_CC->par->cell[i] = CALLOC(xcov->cell_dist->n_o+xcov->cell_dist->n_u,double);
	  if(xcov->cell_dist->n_u>0)
	    {
	      if(i==0)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.int.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.C"));
		}
	      else if(i==1)
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.slp.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.C"));
		}
	      else
		{
		  xcov->cell_dist->n_C = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.hsz.nC"));
		  E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.E"));
		  C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.C"));
		}
	      xcov->cell_dist->E = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_o);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_o;k++)
		  xcov->cell_dist->E[j][k] = E[j*xcov->cell_dist->n_o+k];
	      
	      xcov->cell_dist->C = CALLOC2_d(xcov->cell_dist->n_u,xcov->cell_dist->n_C);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		for(k=0;k<xcov->cell_dist->n_C;k++)
		  xcov->cell_dist->C[j][k] = C[j*xcov->cell_dist->n_C+k];
	    }
	}
      ind++;
    }
  return(0);
}

int re_makedata_cell_dist_CC(SEXP i_dist_cell,int *i_icell,
			     Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC)
{
  int    i;
  Data_cov *xcov;

  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      xcov = i_D_age->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_age->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_age->par->cell);

  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_length->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_length->par->cell);

  for(i=0;i<i_D_lga_CC->glm->nxcov;i++)
    {
      xcov = i_D_lga_CC->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_length_CC->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_length_CC->par->cell);

  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_weight->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_weight->par->cell);

  for(i=0;i<i_D_wgl_CC->glm->nxcov;i++)
    {
      xcov = i_D_wgl_CC->glm->xcov[i];
      if(xcov->icell>=0 && 
	 xcov->cell_dist->n_u>0)
	{
	  FREE(i_weight_CC->par->cell[i]);
	  FREE2_d(xcov->cell_dist->E,xcov->cell_dist->n_u);
	  FREE2_d(xcov->cell_dist->C,xcov->cell_dist->n_u);
	  FREE(xcov->cell_dist);
	}
    }
  FREE(i_weight_CC->par->cell);

  return(0);

}

/*!
  \author Geir Storvik
  \brief  Simulate unobserved cell effects
*/  
int simulate_cell_effects(Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl)
{
  
  int       a,i,j,k,n,ncat;
  double   *eps, *cell, sd;
  Data_cov *xcov;
  
  n = 1;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_age->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_lga->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_wgl->glm->xcov[i]->cell_dist->n_C);
    }
  eps = CALLOC(n,double);
  
  ncat = i_D_age->glm->ncat;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_age->glm->xcov[i];
	  cell = i_age->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_age->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(a=0;a<ncat;a++)
	    {
	      for(k=0;k<n;k++)
		cell[a*n+k] = i_age->par->eff[a][i][xcov->icell][k];
	    }
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[ncat*n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[ncat*n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[ncat*n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga->glm->xcov[i];
	  cell = i_length->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  cell = i_weight->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  FREE(eps);

  return(0);
}

/*!
  \author Geir Storvik
  \brief  Simulate unobserved cell effects
*/  
int simulate_cell_effects_CC(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC)
{
  
  int       a,i,j,k,n,ncat;
  double   *eps, *cell, sd;
  Data_cov *xcov;
  
  n = 1;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_age->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_lga->glm->xcov[i]->cell_dist->n_C);
    }
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>=0)
	n = MAX(n,i_D_wgl->glm->xcov[i]->cell_dist->n_C);
    }
  eps = CALLOC(n,double);
  
  /* age */
  ncat = i_D_age->glm->ncat;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i_D_age->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_age->glm->xcov[i];
	  cell = i_age->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_age->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(a=0;a<ncat;a++)
	    {
	      for(k=0;k<n;k++)
		cell[a*n+k] = i_age->par->eff[a][i][xcov->icell][k];
	    }
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[ncat*n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[ncat*n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[ncat*n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  /* lga */
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      if(i_D_lga->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga->glm->xcov[i];
	  cell = i_length->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* lga - coastal cod */
  for(i=0;i<i_D_lga_CC->glm->nxcov;i++)
    {
      if(i_D_lga_CC->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_lga_CC->glm->xcov[i];
	  cell = i_length_CC->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_length_CC->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_length_CC->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* wgl */
  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      if(i_D_wgl->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl->glm->xcov[i];
	  cell = i_weight->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }
  /* wgl - coastal cod */
  for(i=0;i<i_D_wgl_CC->glm->nxcov;i++)
    {
      if(i_D_wgl_CC->glm->xcov[i]->icell>0)
	{
	  xcov = i_D_wgl_CC->glm->xcov[i];
	  cell = i_weight_CC->par->cell[i];
	  n = xcov->n_fac[xcov->icell];
	  sd = G_ONE/sqrt(i_weight_CC->par->tau[i][i_D_age->glm->xcov[i]->icell]);
	  sd = G_ZERO;
	  // First put observed cell effects into cell vector
	  for(k=0;k<n;k++)
	    cell[k] = i_weight_CC->par->eff[0][i][xcov->icell][k];
	  if(xcov->cell_dist->n_u>0)
	    {
	      // Mean part
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  cell[n+j] = G_ZERO;
		  for(k=0;k<xcov->cell_dist->n_o;k++)
		    cell[n+j] += xcov->cell_dist->E[j][k]*cell[k];
		}
	      // Random part
	      for(k=0;k<xcov->cell_dist->n_C;k++)
		eps[k] = sd*gennor(G_ZERO,G_ONE);
	      for(j=0;j<xcov->cell_dist->n_u;j++)
		{
		  for(k=0;k<xcov->cell_dist->n_C;k++)
		    cell[n+j] += xcov->cell_dist->C[j][k]*eps[k];
		}
	    }
	}
    }

  FREE(eps);

  return(0);
}
