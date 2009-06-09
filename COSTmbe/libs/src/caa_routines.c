/*!
  \file caa_routines.c
  \brief Routines for converting data into approperiate structures and some other general
  routines
*/
#include "caa.h"

static int make_c_cov(int i_n_cov,int i_ispat,int *i_fix,int *i_x_cov,
                      int i_nHaul,Data_cov **o_xcov);
static int re_make_c_cov(int i_n_cov,int i_ispat,int *i_fix,int *i_x_cov,
			 int i_nHaul,Data_cov **o_xcov);
static int make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov);
static int re_make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov);
static int convert_cov(int i_nHaul,int *i_nFac,Data_cov *i_xcov);
static int re_convert_cov(int i_nHaul,int *i_nFac,Data_cov *i_xcov);
static int convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff);
static int re_convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			      Data_cov *o_xcov,int i_neff);
int compare(int *i_x,int *i_y);
static int update_average_par(int i_n,int i_ncat,Eff_str *i_par,Eff_str *x_par_mean,
                              int i_nxcov,Data_cov **i_xcov);

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif

    

/*!
  \brief Writes input for model1 to file caa_input_model1.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik, Hanne Rognebakke
*/
int  write_input_model1_new(SEXP i_mcmc_par,SEXP i_constr,SEXP i_seed,SEXP i_num_par,SEXP i_nBoats, 
			    SEXP i_common_par,SEXP i_dataList,SEXP i_ageList,SEXP i_lgaList,
			    SEXP i_priorList)
{
  SEXP      elmt = R_NilValue;
  int       a,h,i,ind,ind2,isp,n,nFish;
  int       nBoats,nAges,coastal_cod;
  int      *mcmc_par,*num_par,*nFishBoat,*start_noAge,*num_noAge;
  int      *totage,*replength,*totseason,*tottype;
  int      *num_adj_area,*adj_area;
  double   *int_len_lim,*totlength,*haulweight,*A2A;

  int      *a_vec,*age_n_cov,*age_ispat,*hsz_quad;
  int      *age_int_nFac,*age_int_fix,*age_int_c_cov;
  int      *age_hsz_nFac,*age_hsz_fix,*age_hsz_c_cov;

  int      *lga_n_cov,*lga_ispat;
  int      *lga_int_nFac,*lga_int_fix,*lga_int_c_cov;
  int      *lga_slp_nFac,*lga_slp_fix,*lga_slp_c_cov;
  int      *lga_hsz_nFac,*lga_hsz_fix,*lga_hsz_c_cov;

  int       lga_g_a_model,lga_g_a_ncat,lga_fixed_model,lga_cens_model;
  int      *lga_g_a_a2Age_vec;
  double   *lga_g_a_avec,*lga_g_a_par_init;
  double   *lga_fixed_int,*lga_fixed_slp,*lga_fixed_tau;
  double   *lga_fixed_c,*lga_fixed_theta,*lga_fixed_gamma;
  double   *lga_cens;

  double   *pri_eff_mean,*pri_eff_prec,*pri_prec_par,*pri_ar;
  FILE     *caa_input;

  caa_input = fopen("caa_input_model1.txt","w");

  mcmc_par = INTEGER_POINTER(AS_INTEGER(i_mcmc_par));
  fprintf(caa_input,"mcmc_par=%d %d %d\n",mcmc_par[0],mcmc_par[1],mcmc_par[2]);
  fprintf(caa_input,"constr=%d,seed=%d\n",INTEGER_VALUE(i_constr),INTEGER_VALUE(i_seed));
  nBoats = INTEGER_VALUE(i_nBoats);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "nAges")))
    nAges = INTEGER_VALUE(elmt);
  fprintf(caa_input,"nBoats=%d,nAges=%d\n",nBoats,nAges);

  if(!Rf_isNull(elmt = getListElement(i_dataList, "nFishBoat")))
    nFishBoat = INTEGER_POINTER(AS_INTEGER(elmt)); 
  if(!Rf_isNull(elmt = getListElement(i_dataList, "start_noAge")))
    start_noAge = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_dataList, "num_noAge")))
    num_noAge = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_dataList, "totseason")))
    totseason = INTEGER_POINTER(AS_INTEGER(elmt));
  nFish = 0;
  for(h=0;h<nBoats;h++)
    {
      fprintf(caa_input,"i=%d,nFishBoat=%d,start_noAge=%d,num_noAge=%d,season=%d\n",
	      h,nFishBoat[h],start_noAge[h],num_noAge[h],totseason[h]);
      nFish += nFishBoat[h];
    }
  if(!Rf_isNull(elmt = getListElement(i_dataList, "n_int_len")))
    n = INTEGER_VALUE(elmt); 
  if(!Rf_isNull(elmt = getListElement(i_dataList, "int_len_lim")))
    int_len_lim = NUMERIC_POINTER(elmt);
  fprintf(caa_input,"n_int_len_lim=%d\n",n);
  for(i=0;i<n;i++)
    fprintf(caa_input,"%f\n",int_len_lim[i]);

  if(!Rf_isNull(elmt = getListElement(i_dataList, "totage")))
    totage = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_dataList, "totlength")))
    totlength = NUMERIC_POINTER(elmt); 
  if(!Rf_isNull(elmt = getListElement(i_dataList, "replength")))
    replength = INTEGER_POINTER(AS_INTEGER(elmt));
  coastal_cod = INTEGER_VALUE(getListElement(i_dataList, "coastal_cod"));
  if(coastal_cod)
    {
      if(!Rf_isNull(elmt = getListElement(i_dataList, "tottype")))
	tottype = INTEGER_POINTER(AS_INTEGER(elmt));
    }
  if(!Rf_isNull(elmt = getListElement(i_ageList, "a_vec")))
    a_vec = INTEGER_POINTER(AS_INTEGER(elmt));
  fprintf(caa_input,"n=%d,totage[i],totlength[i],replength[i]",nFish);
  if(coastal_cod)
    fprintf(caa_input,",tottype[i]:\n");
  else
    fprintf(caa_input,":\n");
  for(i=0;i<nFish;i++)
    {
      fprintf(caa_input,"i=%d,%d,%f,%d",i,totage[i],totlength[i],replength[i]);
      if(coastal_cod)
	fprintf(caa_input,",%d\n",tottype[i]);
      else
	fprintf(caa_input,"\n");
    }
  fprintf(caa_input,"nAges=%d\n",nAges);
  for(a=0;a<nAges;a++)
    fprintf(caa_input,"a_vec[%d]=%d\n",a,a_vec[a]);

  if(!Rf_isNull(elmt = getListElement(i_ageList, "n_cov")))
    age_n_cov = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ispat")))
    age_ispat = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_ageList, "int_nFac")))
    age_int_nFac = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_ageList, "int_fix")))
    age_int_fix = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_ageList, "int_c_cov")))
    age_int_c_cov = INTEGER_POINTER(AS_INTEGER(elmt));
  fprintf(caa_input,"age_int_n_cov=%d,age_int_ispat=%d\n",age_n_cov[0],age_ispat[0]);
  for(a=0;a<(age_n_cov[0]);a++)
    fprintf(caa_input,"a_nFac[%d]=%d,a_fix[%d]=%d\n",a,age_int_nFac[a],a,age_int_fix[a]);
  fprintf(caa_input,"age_int_c_cov:\n");
  ind = 0;
  for(h=0;h<nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(age_n_cov[0]);i++)
	{
	  fprintf(caa_input,"%d ",age_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(age_n_cov[1]>0)
    {
      fprintf(caa_input,"haulsize\n");
      if(!Rf_isNull(elmt = getListElement(i_dataList, "haulweight")))
	haulweight = NUMERIC_POINTER(elmt); 
      for(h=0;h<nBoats;h++)
	fprintf(caa_input,"i=%d,haulsize=%f\n",h,haulweight[h]);
      if(!Rf_isNull(elmt = getListElement(i_ageList, "hsz_nFac")))
	age_hsz_nFac = INTEGER_POINTER(AS_INTEGER(elmt));
      if(!Rf_isNull(elmt = getListElement(i_ageList, "hsz_fix")))
	age_hsz_fix = INTEGER_POINTER(AS_INTEGER(elmt));
      if(!Rf_isNull(elmt = getListElement(i_ageList, "hsz_c_cov")))
	age_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(elmt));
      fprintf(caa_input,"age_hsz_n_cov=%d,age_hsz_ispat=%d\n",
	      age_n_cov[1],age_ispat[1]);
      for(a=0;a<(age_n_cov[1]);a++)
	fprintf(caa_input,"a_nFac[%d]=%d,a_fix[%d]=%d\n",
		a,age_hsz_nFac[a],a,age_hsz_fix[a]);
      fprintf(caa_input,"age_hsz_c_cov:\n");
      ind = 0;
      for(h=0;h<nBoats;h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(age_n_cov[1]);i++)
	    {
	      fprintf(caa_input,"%d ",age_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  if(!Rf_isNull(elmt = getListElement(i_ageList, "age_errors")))
    n = INTEGER_VALUE(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "A2A")))
    A2A = NUMERIC_POINTER(elmt);
  fprintf(caa_input,"age_errors=%d\n",n);
  if(n)
    {
      fprintf(caa_input,"A2A:\n");
      ind=0;
      for(a=0;a<nAges;a++)
	{
	  for(a=0;a<nAges;a++)
	    {
	      fprintf(caa_input,"%f ",A2A[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  if(!Rf_isNull(elmt = getListElement(i_lgaList, "n_cov")))
    lga_n_cov = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "ispat")))
    lga_ispat = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "int_nFac")))
    lga_int_nFac = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "int_fix")))
    lga_int_fix = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "int_c_cov")))
    lga_int_c_cov = INTEGER_POINTER(AS_INTEGER(elmt));
  fprintf(caa_input,"lga_int_n_cov=%d,lga_int_ispat=%d\n",
	  lga_n_cov[0],lga_ispat[0]);
  for(a=0;a<lga_n_cov[0];a++)
    fprintf(caa_input,"lga_int_nFac[%d]=%d,lga_int_fix[%d]=%d\n",
	    a,lga_int_nFac[a],a,lga_int_fix[a]);
  fprintf(caa_input,"lga_int_c_cov:\n");
  ind = 0;
  for(h=0;h<nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(lga_n_cov[0]);i++)
	{
	  fprintf(caa_input,"%d ",lga_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(!Rf_isNull(elmt = getListElement(i_lgaList, "slp_nFac")))
    lga_slp_nFac = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "slp_fix")))
    lga_slp_fix = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "slp_c_cov")))
    lga_slp_c_cov = INTEGER_POINTER(AS_INTEGER(elmt));
  fprintf(caa_input,"lga_slp_n_cov=%ld,lga_slp_ispat=%d\n",lga_n_cov[1],lga_ispat[1]);
  for(a=0;a<lga_n_cov[1];a++)
    fprintf(caa_input,"lga_slp_nFac[%d]=%d,lga_slp_fix[%d]=%d\n",
	    a,lga_slp_nFac[a],a,lga_slp_fix[a]);
  fprintf(caa_input,"lga_slp_c_cov\n");
  ind = 0;
  for(h=0;h<nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(lga_n_cov[1]);i++)
	{
	  fprintf(caa_input,"%d ",lga_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(lga_n_cov[2]>0)
    {
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "hsz_nFac")))
	lga_hsz_nFac = INTEGER_POINTER(AS_INTEGER(elmt));
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "hsz_fix")))
	lga_hsz_fix = INTEGER_POINTER(AS_INTEGER(elmt));
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "hsz_c_cov")))
	lga_hsz_c_cov = INTEGER_POINTER(AS_INTEGER(elmt));
      fprintf(caa_input,"lga_hsz_n_cov=%d,lga_hsz_ispat=%d\n",lga_n_cov[2],lga_ispat[2]);
      for(a=0;a<lga_n_cov[2];a++)
	fprintf(caa_input,"lga_hsz_nFac[%d]=%d,lga_hsz_fix[%d]=%d\n",
		a,lga_hsz_nFac[a],a,lga_hsz_fix[a]);
      fprintf(caa_input,"lga_hsz_c_cov\n");
      ind = 0;
      for(h=0;h<nBoats;h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(lga_n_cov[2]);i++)
	    {
	      fprintf(caa_input,"%d",lga_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  lga_g_a_model = INTEGER_VALUE(getListElement(i_lgaList, "g_a_model"));
  lga_g_a_ncat = INTEGER_VALUE(getListElement(i_lgaList,"g_a_ncat"));
  lga_g_a_a2Age_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_lgaList,"g_a_a2Age_vec")));
  lga_g_a_avec = NUMERIC_POINTER(getListElement(i_lgaList,"g_a_avec"));
  fprintf(caa_input,"g_a_model=%d\n",lga_g_a_model);
  for(a=0;a<lga_g_a_ncat;a++)
    fprintf(caa_input,"lga_g_a_a_vec[%d]=%f\n",a,lga_g_a_avec[a]);
  if(coastal_cod)
    for(a=0;a<(int)nAges/2;a++)
      fprintf(caa_input,"lga_g_a_a2Age_vec[%d]=%d\n",a,lga_g_a_a2Age_vec[a]);
  else
    for(a=0;a<nAges;a++)
      fprintf(caa_input,"lga_g_a_a2Age_vec[%d]=%d\n",a,lga_g_a_a2Age_vec[a]);
  if(lga_g_a_model == 1)
    {
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "g_a_par_init")))
	lga_g_a_par_init = NUMERIC_POINTER(elmt);
      fprintf(caa_input,"g_a_par_init\n");
      fprintf(caa_input,"c=%f,theta=%f,gamma=%f\n",lga_g_a_par_init[0],
	      lga_g_a_par_init[1],lga_g_a_par_init[2]);
    }

  if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_model")))
    n = INTEGER_VALUE(elmt);
  fprintf(caa_input,"lga_fixed_model=%d\n",n);
  if(n != 0)
    {
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_int")))
	lga_fixed_int = NUMERIC_POINTER(elmt);
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_slp")))
	lga_fixed_slp = NUMERIC_POINTER(elmt);
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_tau")))
	lga_fixed_tau = NUMERIC_POINTER(elmt);
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_g_a_c")))
	lga_fixed_c = NUMERIC_POINTER(elmt);
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_g_a_theta")))
	lga_fixed_theta = NUMERIC_POINTER(elmt);
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "fixed_g_a_gamma")))
	lga_fixed_gamma = NUMERIC_POINTER(elmt);
      fprintf(caa_input,"lga_fixed_int[0]=%f,lga_fixed_slp[0]=%f,lga_fixed_tau[0]=%f\n",
	      lga_fixed_int[0],lga_fixed_slp[0],lga_fixed_tau[0]);
      fprintf(caa_input,"lga_fixed_g_a_c[0]=%f,lga_fixed_g_a_theta[0]=%f,lga_fixed_g_a_gamma[0]=%f\n",
	      lga_fixed_c[0],lga_fixed_theta[0],lga_fixed_gamma[0]);
    }

  n = INTEGER_VALUE(getListElement(i_ageList, "hsz_quad"));
  fprintf(caa_input,"quad:%d\n",n);

  if(INTEGER_VALUE(getListElement(i_lgaList, "cens_model")))
    {
      if(!Rf_isNull(elmt = getListElement(i_lgaList, "cens")))
	lga_cens = NUMERIC_POINTER(elmt);
      fprintf(caa_input,"cens:k=%f,m=%f,r=%f,Nlim=%f\n",lga_cens[0],lga_cens[1],lga_cens[2],lga_cens[3]);
    }

  if(!Rf_isNull(elmt = getListElement(i_lgaList, "num_adj_area")))
    num_adj_area = INTEGER_POINTER(AS_INTEGER(elmt));
  if(!Rf_isNull(elmt = getListElement(i_lgaList, "adj_area")))
    adj_area = INTEGER_POINTER(AS_INTEGER(elmt));
  isp = age_ispat[0];
  if(isp>=0)
    {
      fprintf(caa_input,"age_int_ispat=%d,nFac=%ld\n",isp,age_int_nFac[isp-1]);
      ind = 0;
      for(a=0;a<age_int_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"age_int_ispat:i=%d,num=%d,adj=",a,num_adj_area[a]);
	  for(i=0;i<num_adj_area[a];i++)
	    {
	      fprintf(caa_input,"%d ",adj_area[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(age_n_cov[1]>0)
    {
      isp = age_ispat[1];
      if(isp>=0)
	{
	  fprintf(caa_input,"age_hsz_ispat=%d,nFac=%d\n",isp,age_hsz_nFac[isp-1]);
	  for(a=0;a<age_hsz_nFac[isp-1];a++)
	    fprintf(caa_input,"age_hsz_ispat:i=%d,num=%d,\n",a,num_adj_area[a]);
	}
    }
  isp = lga_ispat[0];
  if(isp>=0)
    {
      fprintf(caa_input,"lga_int_ispat=%d,nFac=%d\n",isp,lga_int_nFac[isp-1]);
      ind = 0;
      for(a=0;a<lga_int_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"lga_int_ispat:i=%d,num=%d,\n",a,num_adj_area[a]);
	  for(i=0;i<num_adj_area[a];i++)
	    {
	      fprintf(caa_input,"%d ",adj_area[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  isp = lga_ispat[1];
  if(isp>=0)
    {
      ind = 0;
      fprintf(caa_input,"lga_slp_ispat=%d,nFac=%d\n",isp,lga_slp_nFac[isp-1]);
      for(a=0;a<lga_slp_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"lga_slp_ispat:i=%d,num=%d, adj=",a,num_adj_area[a]);
	  for(i=0;i<num_adj_area[a];i++)
	    {
	      fprintf(caa_input,"%d ",adj_area[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(lga_n_cov[2]>0)
    {
      isp = lga_ispat[2];
      if(isp>=0)
	fprintf(caa_input,"lga_hsz_ispat=%d,nFac=%d\n",isp,lga_hsz_nFac[isp-1]);
    }

  num_par = INTEGER_POINTER(AS_INTEGER(i_num_par));
  for(a=0;a<length(i_num_par);a++)
      fprintf(caa_input,"i=%d,num_par=%d,\n",a,num_par[a]);

  if(!Rf_isNull(elmt = getListElement(i_priorList, "age_eff_mean")))
    pri_eff_mean = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "age_eff_prec")))
    pri_eff_prec = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "age_prec_par")))
    pri_prec_par = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "age_ar")))
    pri_ar = NUMERIC_POINTER(elmt);
  ind = 0;
  ind2 = 0;
  fprintf(caa_input,"age prior prec+mean\n");
  for(i=0;i<age_n_cov[0];i++)
    {
      if(age_int_fix[i])
	{
	  fprintf(caa_input,"int %d %f ",i,pri_eff_prec[ind2]);
	  ind2++;
	  for(n=0;n<age_int_nFac[i];n++)
	    {
	      fprintf(caa_input,"%f ",pri_eff_mean[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(age_n_cov[1]>0)
    {
      for(i=0;i<age_n_cov[1];i++)
	{
	  if(age_hsz_fix[i])
	    {
	      fprintf(caa_input,"hsz %d %f ",i,pri_eff_prec[ind2]);
	      ind2++;
	      for(n=0;n<age_hsz_nFac[i];n++)
		{
		  fprintf(caa_input,"%f ",pri_eff_mean[ind]);
		  ind++;
		}
	      fprintf(caa_input,"\n");
	    }
	}
    }
  ind = 0;
  fprintf(caa_input,"age prior prec\n");
  for(i=0;i<age_n_cov[0];i++)
    {
      if(!age_int_fix[i])
	{
	  fprintf(caa_input,"int %d %f,%f\n",i,pri_prec_par[ind],pri_prec_par[ind+1]);
	  ind +=2;
	}
    }
  if(age_n_cov[1]>0)
    {
      for(i=0;i<age_n_cov[1];i++)
	{
	  if(!age_hsz_fix[i])
	    {
	      fprintf(caa_input,"hsz %d %f,%f ",i,pri_prec_par[ind],pri_prec_par[ind+1]);
	      ind +=2;
	    }
	}
      fprintf(caa_input,"\n");
    }
  fprintf(caa_input,"age prior ar\n");
  ind = 0;
  if(age_ispat[0]>= 0)
    {
      fprintf(caa_input,"int %f %f\n",pri_ar[ind],pri_ar[ind+1]);
      ind+=2;
    }
  if(age_n_cov[1]>0 && age_ispat[1]>=0)
    fprintf(caa_input,"hsz %f %f\n",pri_ar[ind],pri_ar[ind+1]);

  if(!Rf_isNull(elmt = getListElement(i_priorList, "lga_eff_mean")))
    pri_eff_mean = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "lga_eff_prec")))
    pri_eff_prec = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "lga_prec_par")))
    pri_prec_par = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_priorList, "lga_ar")))
    pri_ar = NUMERIC_POINTER(elmt);
  ind = 0;
  ind2 = 0;
  fprintf(caa_input,"lga prior prec+mean\n");
  for(i=0;i<lga_n_cov[0];i++)
    {
      if(lga_int_fix[i])
	{
	  fprintf(caa_input,"int %d %f ",i,pri_eff_prec[ind2]);
	  ind2++;
	  for(n=0;n<lga_int_nFac[i];n++)
	    {
	      fprintf(caa_input,"%f ",pri_eff_mean[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  for(i=0;i<lga_n_cov[1];i++)
    {
      if(lga_slp_fix[i])
	{
	  fprintf(caa_input,"slp %d %f ",i,pri_eff_prec[ind2]);
	  ind2++;
	  for(n=0;n<lga_slp_nFac[i];n++)
	    {
	      fprintf(caa_input,"%f ",pri_eff_mean[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(lga_n_cov[2]>0)
    {
      for(i=0;i<lga_n_cov[2];i++)
	{
	  if(lga_hsz_fix[i])
	    {
	      fprintf(caa_input,"hsz %d %f ",i,pri_eff_prec[ind2]);
	      ind2++;
	      for(n=0;n<lga_hsz_nFac[i];n++)
		{
		  fprintf(caa_input,"%f ",pri_eff_mean[ind]);
		  ind++;
		}
	      fprintf(caa_input,"\n");
	    }
	}
    }
  ind = 0;
  fprintf(caa_input,"lga prior prec\n");
  for(i=0;i<lga_n_cov[0];i++)
    {
      if(!lga_int_fix[i])
	{
	  fprintf(caa_input,"int %d %f,%f\n",i,pri_prec_par[ind],pri_prec_par[ind+1]);
	  ind +=2;
	}
    }
  for(i=0;i<lga_n_cov[1];i++)
    {
      if(!lga_slp_fix[i])
	{
	  fprintf(caa_input,"slp %d %f,%f\n",i,pri_prec_par[ind],pri_prec_par[ind+1]);
	  ind +=2;
	}
    }
  if(lga_n_cov[2]>0)
    {
      for(i=0;i<lga_n_cov[2];i++)
	{
	  if(!lga_hsz_fix[i])
	    {
	      fprintf(caa_input,"hsz %d %f,%f ",i,pri_prec_par[ind],pri_prec_par[ind+1]);
	      ind +=2;
	    }
	}
      fprintf(caa_input,"\n");
    }
  fprintf(caa_input,"lga prior ar\n");
  ind = 0;
  if(lga_ispat[0]>= 0)
    {
      fprintf(caa_input,"int %f %f\n",pri_ar[ind],pri_ar[ind+1]);
      ind+=2;
    }
  if(lga_ispat[1]>= 0)
    {
      fprintf(caa_input,"slp %f %f\n",pri_ar[ind],pri_ar[ind+1]);
      ind+=2;
    }
  if(lga_n_cov[2]>0 && lga_ispat[2]>=0)
    fprintf(caa_input,"hsz %f %f\n",pri_ar[ind],pri_ar[ind+1]);


  fclose(caa_input);

  return(0);
}		/* end of caa_write_input_model1_new */


    

/*!
  \brief Writes input for model1 to file caa_input_model1.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik
*/
int  write_input_model1(int *i_mcmc_par,int *i_constr,int *i_seed,
			int *nBoats,
			int *i_totage,double *i_totlength,double *i_haulweight,
			int *i_nFishBoat,int *i_replength,
			int *i_start_noAge,int *i_num_noAge,
			int *i_n_int_len,double *i_int_len_lim,
			int *nAges,int *a_vec,int *i_n_cov,int *i_ispat,
                        int *i_age_int_nFac,int *i_age_int_fix,int *i_age_int_c_cov,
                        int *i_age_hsz_nFac,int *i_age_hsz_fix,int *i_age_hsz_c_cov,
			int *i_age_errors,double *i_A2A,
			int *i_lga_int_nFac,int *i_lga_int_fix,int *i_lga_int_c_cov,
			int *i_lga_slp_nFac,int *i_lga_slp_fix,int *i_lga_slp_c_cov,
			int *i_lga_hsz_nFac,int *i_lga_hsz_fix,int *i_lga_hsz_c_cov,
			int *i_lga_g_a_model,double *i_lga_g_a_par_init,
			int *i_lga_fixed_model,
			double *i_lga_fixed_int,double *i_lga_fixed_slp,double *i_lga_fixed_tau,
			double *i_lga_fixed_g_a_c,double *i_lga_fixed_g_a_theta,double *i_lga_fixed_g_a_gamma,
			int *i_quad_hsz,
			int *i_lga_cens_model,double *i_lga_cens,
			int *i_num_adj_area,int *i_adj_area,
			int *i_num_par,
			double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
			double *i_pri_age_prec_par,double *i_pri_age_ar,
			double *i_pri_lga_eff_mean,double *i_pri_lga_eff_prec,
			double *i_pri_lga_prec_par,double *i_pri_lga_ar)
{
  int        a,h,i,ind,ind2,j,isp,k,n;
  FILE      *caa_input;

  caa_input = fopen("caa_input_model1.txt","w");
  fprintf(caa_input,"mcmc_par=%ld %ld %ld\n",i_mcmc_par[0],i_mcmc_par[1],i_mcmc_par[2]);
  fprintf(caa_input,"constr=%ld,seed=%ld\n",*i_constr,*i_seed);
  fprintf(caa_input,"nBoats=%ld,nAges=%ld\n",
	  *nBoats,*nAges);
  n = 0;
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"i=%d,nFishBoat=%ld,start_noAge=%ld,num_noAge=%ld\n",
	      h,i_nFishBoat[h],i_start_noAge[h],i_num_noAge[h]);
      n += i_nFishBoat[h];
    }
  fprintf(caa_input,"n_int_len_lim=%ld\n",*i_n_int_len);
  for(i=0;i<(*i_n_int_len);i++)
    fprintf(caa_input,"%lf\n",i_int_len_lim[i]);

  fprintf(caa_input,"n=%d,totage[i],totlength[i],replength[i]:\n",n);
  for(i=0;i<n;i++)
    fprintf(caa_input,"i=%d,%ld,%lf,%ld\n",i,i_totage[i],i_totlength[i],i_replength[i]);
  fprintf(caa_input,"nAges=%ld\n",*nAges);
  for(a=0;a<(*nAges);a++)
    fprintf(caa_input,"a_vec[%d]=%ld\n",a,a_vec[a]);

  fprintf(caa_input,"age_int_n_cov=%ld,i_age_int_ispat=%ld\n",i_n_cov[0],i_ispat[0]);
  for(a=0;a<(i_n_cov[0]);a++)
    fprintf(caa_input,"a_nFac[%d]=%ld,a_fix[%d]=%ld\n",a,i_age_int_nFac[a],a,i_age_int_fix[a]);
  fprintf(caa_input,"age_int_c_cov:\n");
  ind = 0;
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[0]);i++)
	{
	  fprintf(caa_input,"%ld ",i_age_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(i_n_cov[1]>0)
    {
      fprintf(caa_input,"haulsize\n");
      for(h=0;h<(*nBoats);h++)
	fprintf(caa_input,"i=%d,haulsize=%lf\n",h,i_haulweight[h]);
      fprintf(caa_input,"age_hsz_n_cov=%ld,i_age_hsz_ispat=%ld\n",
	      i_n_cov[1],i_ispat[1]);
      for(a=0;a<(i_n_cov[1]);a++)
	fprintf(caa_input,"a_nFac[%d]=%ld,a_fix[%d]=%ld\n",
		a,i_age_hsz_nFac[a],a,i_age_hsz_fix[a]);
      fprintf(caa_input,"age_hsz_c_cov:\n");
      ind = 0;
      for(h=0;h<(*nBoats);h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(i_n_cov[1]);i++)
	    {
	      fprintf(caa_input,"%ld ",i_age_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  fprintf(caa_input,"age_errors=%ld\n",*i_age_errors);
  if(*i_age_errors)
    {
      fprintf(caa_input,"A2A:\n");
      ind=0;
      for(a=0;a<(*nAges);a++)
	{
	  for(a=0;a<(*nAges);a++)
	    {
	      fprintf(caa_input,"%lf ",i_A2A[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  fprintf(caa_input,"lga_int_n_cov=%ld,lga_int_ispat=%ld\n",
	  i_n_cov[2],i_ispat[2]);
  for(a=0;a< i_n_cov[2];a++)
    fprintf(caa_input,"lga_int_nFac[%d]=%ld,lga_int_fix[%d]=%ld\n",
	    a,i_lga_int_nFac[a],a,i_lga_int_fix[a]);
  fprintf(caa_input,"lga_int_c_cov:\n");
  ind = 0;
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[2]);i++)
	{
	  fprintf(caa_input,"%ld ",i_lga_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_slp_n_cov=%ld,lga_slp_ispat=%ld\n",i_n_cov[3],i_ispat[3]);
  for(a=0;a< i_n_cov[3];a++)
    fprintf(caa_input,"lga_slp_nFac[%d]=%ld,lga_slp_fix[%d]=%ld\n",
	    a,i_lga_slp_nFac[a],a,i_lga_slp_fix[a]);
  fprintf(caa_input,"lga_slp_c_cov\n");
  ind = 0;
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[3]);i++)
	{
	  fprintf(caa_input,"%ld ",i_lga_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(i_n_cov[4]>0)
    {
      fprintf(caa_input,"lga_hsz_n_cov=%ld,lga_hsz_ispat=%ld\n",i_n_cov[4],i_ispat[4]);
      for(a=0;a< i_n_cov[4];a++)
	fprintf(caa_input,"lga_hsz_nFac[%d]=%ld,lga_hsz_fix[%d]=%ld\n",
		a,i_lga_hsz_nFac[a],a,i_lga_hsz_fix[a]);
      fprintf(caa_input,"lga_hsz_c_cov\n");
      ind = 0;
      for(h=0;h<(*nBoats);h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(i_n_cov[4]);i++)
	    {
	      fprintf(caa_input,"%ld",i_lga_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  fprintf(caa_input,"g_a_model=%ld\n",*i_lga_g_a_model);
  if(*i_lga_g_a_model == 1)
    {
      fprintf(caa_input,"g_a_par_init\n");
      fprintf(caa_input,"c=%lf,theta=%lf,gamma=%lf\n",i_lga_g_a_par_init[0],i_lga_g_a_par_init[1],i_lga_g_a_par_init[2]);
    }
  fprintf(caa_input,"lga_fixed_model=%ld\n",*i_lga_fixed_model);
  if(*i_lga_fixed_model != 0)
    {
      fprintf(caa_input,"lga_fixed_int[0]=%lf,lga_fixed_slp[0]=%lf,lga_fixed_tau[0]=%lf\n",
	      i_lga_fixed_int[0],i_lga_fixed_slp[0],i_lga_fixed_tau[0]);
      fprintf(caa_input,"lga_fixed_g_a_c[0]=%lf,lga_fixed_g_a_theta[0]=%lf,lga_fixed_g_a_gamma[0]=%lf\n",
	      i_lga_fixed_g_a_c[0],i_lga_fixed_g_a_theta[0],i_lga_fixed_g_a_gamma[0]);
    }
  fprintf(caa_input,"quad:%ld\n",*i_quad_hsz);
  if(*i_lga_cens_model)
    fprintf(caa_input,"cens:k=%lf,m=%lf,r=%lf,Nlim=%lf\n",i_lga_cens[0],i_lga_cens[1],i_lga_cens[2],i_lga_cens[3]);
  isp = i_ispat[0];
  if(isp>=0)
    {
      fprintf(caa_input,"age_int_ispat=%d,nFac=%ld\n",isp,i_age_int_nFac[isp-1]);
      i = 0;
      for(a=0;a<i_age_int_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"age_int_ispat:i=%d,num=%ld,adj=",a,i_num_adj_area[a]);
	  for(j=0;j<i_num_adj_area[a];j++)
	    {
	      fprintf(caa_input,"%ld ",i_adj_area[i]);
	      i++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(i_n_cov[1]>0)
    {
      isp = i_ispat[1];
      if(isp>=0)
	{
	  fprintf(caa_input,"age_hsz_ispat=%d,nFac=%ld\n",isp,i_age_hsz_nFac[isp-1]);
	  for(a=0;a<i_age_hsz_nFac[isp-1];a++)
	    fprintf(caa_input,"age_hsz_ispat:i=%d,num=%ld,\n",a,i_num_adj_area[a]);
	}
    }
  isp = i_ispat[2];
  if(isp>=0)
    {
      fprintf(caa_input,"lga_int_ispat=%d,nFac=%ld\n",isp,i_lga_int_nFac[isp-1]);
      for(a=0;a<i_lga_int_nFac[isp-1];a++)
	fprintf(caa_input,"lga_int_ispat:i=%d,num=%ld,\n",a,i_num_adj_area[a]);
    }
  isp = i_ispat[3];
  if(isp>=0)
    {
      i = 0;
      fprintf(caa_input,"lga_slp_ispat=%d,nFac=%ld\n",isp,i_lga_slp_nFac[isp-1]);
      for(a=0;a<i_lga_slp_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"lga_slp_ispat:i=%d,num=%ld, adj=",a,i_num_adj_area[a]);
	  for(j=0;j<i_num_adj_area[a];j++)
	    {
	      fprintf(caa_input,"%ld ",i_adj_area[i]);
	      i++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(i_n_cov[4]>0)
    {
      isp = i_ispat[4];
      if(isp>=0)
	fprintf(caa_input,"lga_hsz_ispat=%d,nFac=%ld\n",isp,i_lga_hsz_nFac[isp-1]);
    }
  for(a=0;a<3;a++)
      fprintf(caa_input,"i=%d,num_par=%ld,\n",a,i_num_par[a]);

  ind = 0;
  ind2 = 0;
  fprintf(caa_input,"age prior prec+mean\n");
  for(j=0;j<i_n_cov[0];j++)
    {
      if(i_age_int_fix[j])
	{
	  fprintf(caa_input,"int %d %lf ",j,i_pri_age_eff_prec[ind2]);
	  ind2++;
	  for(k=0;k<i_age_int_nFac[j];k++)
	    {
	      fprintf(caa_input,"%lf ",i_pri_age_eff_mean[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(i_n_cov[1]>0)
    {
      for(j=0;j<i_n_cov[1];j++)
	{
	  if(i_age_hsz_fix[j])
	    {
	      fprintf(caa_input,"hsz %d %lf ",j,i_pri_age_eff_prec[ind2]);
	      ind2++;
	      for(k=0;k<i_age_hsz_nFac[j];k++)
		{
		  fprintf(caa_input,"%lf ",i_pri_age_eff_mean[ind]);
		  ind++;
		}
	      fprintf(caa_input,"\n");
	    }
	}
    }
  ind = 0;
  fprintf(caa_input,"age prior prec\n");
  for(j=0;j<i_n_cov[0];j++)
    {
      if(!i_age_int_fix[j])
	{
	  fprintf(caa_input,"int %d %lf,%lf\n",j,i_pri_age_prec_par[ind],i_pri_age_prec_par[ind+1]);
	  ind +=2;
	}
    }
  if(i_n_cov[1]>0)
    {
      for(j=0;j<i_n_cov[1];j++)
	{
	  if(!i_age_hsz_fix[j])
	    {
	      fprintf(caa_input,"hsz %d %lf,%lf ",j,i_pri_age_prec_par[ind],i_pri_age_prec_par[ind+1]);
	      ind +=2;
	    }
	}
      fprintf(caa_input,"\n");
    }
  fprintf(caa_input,"age prior ar\n");
  ind = 0;
  if(i_ispat[0]>= 0)
    {
      fprintf(caa_input,"int %lf %lf\n",i_pri_age_ar[ind],i_pri_age_ar[ind+1]);
      ind+=2;
    }
  if(i_n_cov[1]>0 && i_ispat[1]>=0)
    fprintf(caa_input,"hsz %lf %lf\n",i_pri_age_ar[ind],i_pri_age_ar[ind+1]);
  fclose(caa_input);
  

  return(0);
}		/* end of caa_write_input_model1 */

    

/*!
  \brief Writes input for model2 to file caa_input_model2.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik
*/
int write_input_model2(int *i_mcmc_par,
		       int *i_constr,int *i_seed,
		       int *i_wgl_nBoats,double *i_totlength,double *i_totweight,double *i_haulweight,
		       int *i_replength,int *i_nFishBoat,int *i_n_cov,int *i_ispat,
		       int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
		       int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
		       int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
		       int *i_num_adj_area,int *i_adj_area,
		       int *i_num_par)
{
  int        a,h,i,isp,j,n;

  FILE      *caa_input;

  caa_input = fopen("caa_input_model2.txt","w");
  fprintf(caa_input,"mcmc_par = %d %d %d\n",i_mcmc_par[0],i_mcmc_par[1],i_mcmc_par[2]);
  fprintf(caa_input,"constr=%d, seed=%d\n",*i_constr,*i_seed);
  fprintf(caa_input,"nBoats=%d\n",*i_wgl_nBoats);
  n = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"i=%d,nFishBoat=%d\n",h,i_nFishBoat[h]);
      n += i_nFishBoat[h];
    }
  for(i=0;i<n;i++)
    fprintf(caa_input,"%f %f %d\n",i_totlength[i],i_totweight[i],i_replength[i]);
  fprintf(caa_input,"wgl_nBoats=%d\n",*i_wgl_nBoats);

  fprintf(caa_input,"wgl_int_n_cov=%d,wgl_int_ispat=%d\n",i_n_cov[0],i_ispat[0]);
  for(a=0;a< i_n_cov[0];a++)
    fprintf(caa_input,"wgl_int_nFac[%d]=%d,wgl_int_fix[%d]=%d\n",a,i_wgl_int_nFac[a],a,i_wgl_int_fix[a]);
  for(h=0;h<(*i_wgl_nBoats);h++)
    fprintf(caa_input,"wgl_int_c_cov[%d]=%d\n",h,i_wgl_int_c_cov[h]);

  fprintf(caa_input,"wgl_slp_n_cov=%d,wgl_slp_ispat=%d\n",i_n_cov[1],i_ispat[1]);
  for(a=0;a< i_n_cov[1];a++)
    fprintf(caa_input,"wgl_slp_nFac[%d]=%d,wgl_slp_fix[%d]=%d\n",a,i_wgl_slp_nFac[a],a,i_wgl_slp_fix[a]);
  for(h=0;h<(*i_wgl_nBoats);h++)
    fprintf(caa_input,"wgl_c_cov[%d]=%d\n",h,i_wgl_slp_c_cov[h]);

  if(i_n_cov[2]>0)
    {
      fprintf(caa_input,"haulsize\n");
      for(h=0;h<(*i_wgl_nBoats);h++)
	fprintf(caa_input,"i=%d,haulsize=%f\n",h,i_haulweight[h]);
      fprintf(caa_input,"wgl_hsz_n_cov=%d,wgl_hsz_ispat=%d\n",i_n_cov[2],i_ispat[2]);
      for(a=0;a< i_n_cov[2];a++)
	fprintf(caa_input,"wgl_hsz_nFac[%d]=%d,wgl_hsz_fix[%d]=%d\n",a,i_wgl_hsz_nFac[a],a,i_wgl_hsz_fix[a]);
      for(h=0;h<(*i_wgl_nBoats);h++)
	fprintf(caa_input,"wgl_c_cov[%d]=%d\n",h,i_wgl_hsz_c_cov[h]);
    }
  h = 0;
  isp = i_ispat[0];
  if(isp>=0)
    {
      fprintf(caa_input,"int_ispat=%d,nFac=%d\n",isp,i_wgl_int_nFac[isp-1]);
      i = 0;
      for(a=0;a<i_wgl_int_nFac[isp-1];a++)
	{
	  fprintf(caa_input,"int_ispat:i=%d,num=%d, adj=",a,i_num_adj_area[a]);
	  for(j=0;j<i_num_adj_area[a];j++)
	    {
	      fprintf(caa_input,"%d ",i_adj_area[i]);
	      i++;
	    }
	  fprintf(caa_input,"\n");
	}
    }
  if(i_n_cov[2]>0)
    {
      isp = i_ispat[1];
      if(isp>=0)
	fprintf(caa_input,"int_hsz_ispat=%d,nFac=%d\n",isp,i_wgl_hsz_nFac[isp-1]);
    }

  fclose(caa_input);
  
  return(0);
}


/*!
  \brief Writes input for prediction to file caa_input_predict.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik
*/
/*F:*

________________________________________________________________

		
________________________________________________________________

Name:		 write_input_predict
Syntax:	 See call below
Description:    Writes input for model1 to file caa_input_predict.txt in the
                current directory. Only to be used for testing.
Side effects:   A file caa_input_model1.txt is made
Return value:   0
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_routines.c,v 1.1 2009/06/09 09:27:06 mmerzere Exp $
________________________________________________________________
*/
int write_input_predict(int *i_Nmcmc,int *i_burnin,double *i_mcmc1,double *i_mcmc2,
			int *i_num_par1,int *i_num_par2,
			int *nBoats,int *nAges,
			int *a_vec,
			int *i_n_cov,int *i_ispat,
			int *i_age_int_nFac,int *i_age_int_fix,int *i_age_int_c_cov,
			int *i_age_hsz_nFac,int *i_age_hsz_fix,int *i_age_hsz_c_cov,
			int *i_lga_nBoats,
			int *i_lga_int_nFac,int *i_lga_int_fix,int *i_lga_int_c_cov,
			int *i_lga_slp_nFac,int *i_lga_slp_fix,int *i_lga_slp_c_cov,
			int *i_lga_hsz_nFac,int *i_lga_hsz_fix,int *i_lga_hsz_c_cov,
			int *i_lga_g_a_model,
			int *i_lga_cens_model,double *i_lga_cens,
			int *i_wgl_nBoats,
			int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
			int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
			int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
			int *i_num_adj_area,int *i_adj_area,
			int *i_tot_nCell,int *i_tot_nFactors,
			int *i_tot_fac_age_int,int *i_tot_fac_age_hsz,
			int *i_tot_fac_lga_int,int *i_tot_fac_lga_slp,int *i_tot_fac_lga_hsz,
			int *i_tot_fac_wgl_int,int *i_tot_fac_wgl_slp,int *i_tot_fac_wgl_hsz,
			int *i_inc_haul,int *i_tot_factors,double *i_tot_catch,double *i_par_haulsize,
			int *i_N_l_int,double *i_l_int,int *i_nMC)
{
  int        a,c,h,i,j,ind;
  FILE      *caa_input;

  caa_input = fopen("caa_input_predict.txt","w");

  fprintf(caa_input,"Nmcmc = %d burnin= %d\n",*i_Nmcmc,*i_burnin);

  fprintf(caa_input,"nBoats_age=%d,\n",*nBoats);
  for(a=0;a<(*nAges);a++)
    fprintf(caa_input,"i=%d,avec=%d\n",a,a_vec[a]);

  fprintf(caa_input,"age_int_ncov=%d,i_age_int_ispat=%d\n",i_n_cov[0],i_ispat[0]);
  for(i=0;i<i_n_cov[0];i++)
    fprintf(caa_input,"i=%d,i_age_int_nFac=%d,i_age_int_fix=%d\n",i,i_age_int_nFac[i],i_age_int_fix[i]);
  ind = 0;
  fprintf(caa_input,"age_int_c_cov\n");
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<i_n_cov[0];i++)
	{
	  fprintf(caa_input,"%d ",i_age_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"age_hsz_ncov=%d,i_age_hsz_ispat=%d\n",i_n_cov[1],i_ispat[1]);
  for(i=0;i<i_n_cov[1];i++)
    fprintf(caa_input,"i=%d,i_age_hsz_nFac=%d,i_age_hsz_fix=%d\n",i,i_age_hsz_nFac[i],i_age_hsz_fix[i]);
  ind = 0;
  fprintf(caa_input,"age_hsz_c_cov\n");
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<i_n_cov[1];i++)
	{
	  fprintf(caa_input,"%d ",i_age_hsz_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }


  fprintf(caa_input,"i_lga_nBoats=%d\n",*i_lga_nBoats);
  fprintf(caa_input,"lga_int_ncov=%d,lga_int_ispat=%d\n",i_n_cov[2],i_ispat[2]);
  for(i=0;i<(i_n_cov[2]);i++)
    fprintf(caa_input,"i=%d,lga_int_nFac=%d,lga_int_fix=%d\n",i,i_lga_int_nFac[i],i_lga_int_fix[i]);
  fprintf(caa_input,"lga_int_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[2]);i++)
	{
	  fprintf(caa_input,"%d ",i_lga_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_slp_ncov=%d,lga_slp_ispat=%d\n",i_n_cov[3],i_ispat[3]);
  for(i=0;i<(i_n_cov[3]);i++)
    fprintf(caa_input,"i=%d,lga_slp_nFac=%d,lga_slp_fix=%d\n",i,i_lga_slp_nFac[i],i_lga_slp_fix[i]);
  fprintf(caa_input,"lga_slp_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[3]);i++)
	{
	  fprintf(caa_input,"%d ",i_lga_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_hsz_ncov=%d,lga_hsz_ispat=%d\n",i_n_cov[4],i_ispat[4]);
  for(i=0;i<(i_n_cov[4]);i++)
    fprintf(caa_input,"i=%d,lga_hsz_nFac=%d,lga_hsz_fix=%d\n",i,i_lga_hsz_nFac[i],i_lga_hsz_fix[i]);
  fprintf(caa_input,"lga_hsz_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[4]);i++)
	{
	  fprintf(caa_input,"%d ",i_lga_hsz_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"g_a_model=%d\n",*i_lga_g_a_model);
  if(*i_lga_cens_model)
    fprintf(caa_input,"cens:k=%f,m=%f,r=%f,Nlim=%f\n",i_lga_cens[0],i_lga_cens[1],i_lga_cens[2],i_lga_cens[3]);


  fprintf(caa_input,"i_wgl_nBoats=%d\n",
	  *i_wgl_nBoats);
  fprintf(caa_input,"wgl_int_ncov=%d,wgl_int_ispat=%d\n",i_n_cov[5],i_ispat[5]);

  for(i=0;i<(i_n_cov[5]);i++)
    fprintf(caa_input,"i=%d,wgl_int_nFac=%d,wgl_int_fix=%d\n",i,i_wgl_int_nFac[i],i_wgl_int_fix[i]);
  fprintf(caa_input,"wgl_int_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[5]);i++)
	{
	  fprintf(caa_input,"%d ",i_wgl_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_slp_ncov=%d,wgl_slp_ispat=%d\n",i_n_cov[6],i_ispat[6]);
  for(i=0;i<(i_n_cov[6]);i++)
    fprintf(caa_input,"i=%d,wgl_slp_nFac=%d,wgl_slp_fix=%d\n",i,i_wgl_slp_nFac[i],i_wgl_slp_fix[i]);
  fprintf(caa_input,"wgl_slp_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[6]);i++)
	{
	  fprintf(caa_input,"%d ",i_wgl_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_hsz_ncov=%d,wgl_hsz_ispat=%d\n",i_n_cov[7],i_ispat[7]);
  for(i=0;i<(i_n_cov[7]);i++)
    fprintf(caa_input,"i=%d,wgl_hsz_nFac=%d,wgl_hsz_fix=%d\n",i,i_wgl_hsz_nFac[i],i_wgl_hsz_fix[i]);
  fprintf(caa_input,"wgl_hsz_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[7]);i++)
	{
	  fprintf(caa_input,"%d ",i_wgl_hsz_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(i_ispat[0] > -1 && i_age_int_nFac[i_ispat[0]]>0)
    {
      fprintf(caa_input,"Neighborhood\n"); 
      ind = 0;
      for(i=0;i<(i_age_int_nFac[i_ispat[0]-1]);i++)
	 {
	   fprintf(caa_input,"%d ",i_num_adj_area[i]);
	   for(j=0;j<i_num_adj_area[i];j++)
	     {
	       fprintf(caa_input,"%d ",i_adj_area[ind]);
	       ind++;
	     }
	   fprintf(caa_input,"\n");
	 }
    }

  fprintf(caa_input,"tot_nCell=%d, tot_nFactors=%d\n",*i_tot_nCell,*i_tot_nFactors);

  fprintf(caa_input,"tot_fac_age_int:");
  for(i=0;i<i_n_cov[0];i++)
    fprintf(caa_input,"%d ",i_tot_fac_age_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_age_hsz:");
  for(i=0;i<i_n_cov[1];i++)
    fprintf(caa_input,"%d ",i_tot_fac_age_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"tot_fac_lga_int:");
  for(i=0;i<i_n_cov[2];i++)
    fprintf(caa_input,"%d ",i_tot_fac_lga_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_lga_slp:");
  for(i=0;i<i_n_cov[3];i++)
    fprintf(caa_input,"%d ",i_tot_fac_lga_slp[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_lga_hsz:");
  for(i=0;i<i_n_cov[4];i++)
    fprintf(caa_input,"%d ",i_tot_fac_lga_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"tot_fac_wgl_int:");
  for(i=0;i<i_n_cov[2];i++)
    fprintf(caa_input,"%d ",i_tot_fac_wgl_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_wgl_slp:");
  for(i=0;i<i_n_cov[3];i++)
    fprintf(caa_input,"%d ",i_tot_fac_wgl_slp[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_wgl_hsz:");
  for(i=0;i<i_n_cov[4];i++)
    fprintf(caa_input,"%d ",i_tot_fac_wgl_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"inc_haul= %d %d %d %d %d %d %d\n",
	  i_inc_haul[0],i_inc_haul[1],i_inc_haul[2],i_inc_haul[3],
	  i_inc_haul[4],i_inc_haul[5],i_inc_haul[6]);
  fprintf(caa_input,"tot_factors catch\n");
  ind = 0;
  for(c=0;c<(*i_tot_nCell);c++)
    {
      fprintf(caa_input,"%d ",c);
      for(i=0;i<(*i_tot_nFactors);i++)
	{
	  fprintf(caa_input,"%d ",i_tot_factors[ind]);
	  ind++;
	}
      fprintf(caa_input,"%f\n",i_tot_catch[c]);
    }
  if((i_n_cov[1]+i_n_cov[4]+i_n_cov[7])>0)
    {
      fprintf(caa_input,"parameters for haulsize\n");
      for(c=0;c<(*i_tot_nCell);c++)
	fprintf(caa_input,"%f %f\n",i_par_haulsize[2*c],i_par_haulsize[2*c+1]);
    }
     
  fprintf(caa_input,"N_l_int=%d\n",*i_N_l_int);
  fprintf(caa_input,"l_int\n");
  for(i=0;i<(*i_N_l_int);i++)
    fprintf(caa_input,"%f\n",i_l_int[i]);

  fprintf(caa_input,"nMC=%d\n",*i_nMC);
  fclose(caa_input);

  return(0);
}		/* end of caa_write_input_predict */


/*!
  \brief Writes input for prediction to file caa_input_predict.txt in the current directory. 

  Only to be used for testing.
  \author Geir Storvik
*/
/*F:*

________________________________________________________________

		
________________________________________________________________

Name:		 write_input_predict_new
Syntax:	 See call below
Description:    Writes input for model1 to file caa_input_predict.txt in the
                current directory. Only to be used for testing.
Side effects:   A file caa_input_model1.txt is made
Return value:   0
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_routines.c,v 1.1 2009/06/09 09:27:06 mmerzere Exp $
________________________________________________________________
*/
int write_input_predict_new(SEXP i_mcmc_samp,SEXP i_common_par,
			    SEXP i_data_age,SEXP i_data_lga,SEXP i_data_wgl,
			    SEXP i_data_catch,
			    SEXP i_par_haulsize,SEXP i_dist_cell,
			    SEXP i_N_l_int,SEXP i_l_int,SEXP i_nMC)
{
  int        a,c,h,i,j,ind;
  FILE      *caa_input;

  /* Variables connect to input data */
  int     nMCMC, burnin, nMC;  
  double *mcmc1=NULL, *mcmc2=NULL;
  int *num_par1, *num_par2;
  int  age_nBoats, nAges, *a_vec;
  int *n_cov, *ispat, *icell;
  int *age_int_nFac, *age_int_fix, *age_int_c_cov;
  int *age_hsz_nFac, *age_hsz_fix, *age_hsz_c_cov;
  int  lga_nBoats;
  int *lga_int_nFac, *lga_int_fix, *lga_int_c_cov;
  int *lga_slp_nFac, *lga_slp_fix, *lga_slp_c_cov;
  int *lga_hsz_nFac, *lga_hsz_fix, *lga_hsz_c_cov;
  int  lga_g_a_model, lga_cens_model;
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
  int *num_cell_o, *num_cell_u;
  int  age_int_cell_nC,age_hsz_cell_nC;
  int  lga_int_cell_nC,lga_slp_cell_nC,lga_hsz_cell_nC;
  int  wgl_int_cell_nC,wgl_slp_cell_nC,wgl_hsz_cell_nC;
  double *age_int_cell_E=NULL, *age_hsz_cell_E=NULL;
  double *lga_int_cell_E=NULL, *lga_slp_cell_E=NULL, *lga_hsz_cell_E=NULL;
  double *wgl_int_cell_E=NULL, *wgl_slp_cell_E=NULL, *wgl_hsz_cell_E=NULL;
  double *age_int_cell_C=NULL, *age_hsz_cell_C=NULL;
  double *lga_int_cell_C=NULL, *lga_slp_cell_C=NULL, *lga_hsz_cell_C=NULL;
  double *wgl_int_cell_C=NULL, *wgl_slp_cell_C=NULL, *wgl_hsz_cell_C=NULL;
  int N_l_int,  n_MC;
  int  *inc_haul;
  double *l_int;

  caa_input = fopen("caa_input_predict.txt","w");

  nMCMC = INTEGER_VALUE(getListElement(i_mcmc_samp,"nMCMC"));
  burnin = INTEGER_VALUE(getListElement(i_mcmc_samp,"burnin"));
  mcmc1 = NUMERIC_POINTER(getListElement(i_mcmc_samp,"samples1"));
  mcmc2 = NUMERIC_POINTER(getListElement(i_mcmc_samp,"samples2"));
  num_par1 = INTEGER_POINTER(AS_INTEGER(getListElement(i_mcmc_samp,"numpar1")));
  num_par2 = INTEGER_POINTER(AS_INTEGER(getListElement(i_mcmc_samp,"numpar2")));

  n_cov = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"ncov")));
  ispat = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"ispat")));
  icell = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"icell")));
  num_adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"neigh.num")));
  adj_area = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"neigh.adj")));
  inc_haul = INTEGER_POINTER(AS_INTEGER(getListElement(i_common_par,"inchaul")));

  age_nBoats = INTEGER_VALUE(getListElement(i_data_age,"nBoats"));
  nAges = INTEGER_VALUE(getListElement(i_data_age,"nAges"));
  a_vec = INTEGER_POINTER(AS_INTEGER(getListElement(i_data_age,"avec")));
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
  
  num_cell_o =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.o")));
  num_cell_u =  INTEGER_POINTER(AS_INTEGER(getListElement(i_dist_cell,"num.cell.u")));
  age_int_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.E"));
  age_hsz_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.E"));
  lga_int_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.E"));
  lga_slp_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.E"));
  lga_hsz_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.E"));
  wgl_int_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.E"));
  wgl_slp_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.E"));
  wgl_hsz_cell_E = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.E"));
  age_int_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"age.int.nC"));
  age_hsz_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"age.hsz.nC"));
  lga_int_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"lga.int.nC"));
  lga_slp_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"lga.slp.nC"));
  lga_hsz_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"lga.hsz.nC"));
  wgl_int_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.int.nC"));
  wgl_slp_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.slp.nC"));
  wgl_hsz_cell_nC = INTEGER_VALUE(getListElement(i_dist_cell,"wgl.hsz.nC"));
  age_int_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.int.C"));
  age_hsz_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"age.hsz.C"));
  lga_int_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.int.C"));
  lga_slp_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.slp.C"));
  lga_hsz_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"lga.hsz.C"));
  wgl_int_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.int.C"));
  wgl_slp_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.slp.C"));
  wgl_hsz_cell_C = NUMERIC_POINTER(getListElement(i_dist_cell,"wgl.hsz.C"));

  nMC = INTEGER_VALUE(i_nMC);
  N_l_int = INTEGER_VALUE(i_N_l_int);
  l_int = NUMERIC_POINTER(i_l_int);


  fprintf(caa_input,"Nmcmc = %d burnin= %d\n",nMCMC,burnin);

  fprintf(caa_input,"nBoats_age=%d,\n",age_nBoats);
  for(a=0;a<nAges;a++)
    fprintf(caa_input,"i=%d,avec=%d\n",a,a_vec[a]);

  fprintf(caa_input,"age_int_ncov=%d,age_int_ispat=%d\n",n_cov[0],ispat[0]);
  for(i=0;i<n_cov[0];i++)
    fprintf(caa_input,"i=%d,age_int_nFac=%d,age_int_fix=%d\n",i,age_int_nFac[i],age_int_fix[i]);
  ind = 0;
  fprintf(caa_input,"age_int_c_cov\n");
  for(h=0;h<age_nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<n_cov[0];i++)
	{
	  fprintf(caa_input,"%d ",age_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"age_hsz_ncov=%d,age_hsz_ispat=%d\n",n_cov[1],ispat[1]);
  for(i=0;i<n_cov[1];i++)
    fprintf(caa_input,"i=%d,age_hsz_nFac=%d,age_hsz_fix=%d\n",i,age_hsz_nFac[i],age_hsz_fix[i]);
  if(n_cov[1]>0)
    {
      ind = 0;
      fprintf(caa_input,"age_hsz_c_cov\n");
      for(h=0;h<age_nBoats;h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<n_cov[1];i++)
	    {
	      fprintf(caa_input,"%d ",age_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }


  fprintf(caa_input,"lga_nBoats=%d\n",lga_nBoats);
  fprintf(caa_input,"lga_int_ncov=%d,lga_int_ispat=%d\n",n_cov[2],ispat[2]);
  for(i=0;i<(n_cov[2]);i++)
    fprintf(caa_input,"i=%d,lga_int_nFac=%d,lga_int_fix=%d\n",i,lga_int_nFac[i],lga_int_fix[i]);
  fprintf(caa_input,"lga_int_c_cov\n");
  ind = 0;
  for(h=0;h<lga_nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(n_cov[2]);i++)
	{
	  fprintf(caa_input,"%d ",lga_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_slp_ncov=%d,lga_slp_ispat=%d\n",n_cov[3],ispat[3]);
  for(i=0;i<(n_cov[3]);i++)
    fprintf(caa_input,"i=%d,lga_slp_nFac=%d,lga_slp_fix=%d\n",i,lga_slp_nFac[i],lga_slp_fix[i]);
  fprintf(caa_input,"lga_slp_c_cov\n");
  ind = 0;
  for(h=0;h<lga_nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(n_cov[3]);i++)
	{
	  fprintf(caa_input,"%d ",lga_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_hsz_ncov=%d,lga_hsz_ispat=%d\n",n_cov[4],ispat[4]);
  for(i=0;i<(n_cov[4]);i++)
    fprintf(caa_input,"i=%d,lga_hsz_nFac=%d,lga_hsz_fix=%d\n",i,lga_hsz_nFac[i],lga_hsz_fix[i]);
  if(n_cov[4]>0)
    {
      fprintf(caa_input,"lga_hsz_c_cov\n");
      ind = 0;
      for(h=0;h<lga_nBoats;h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(n_cov[4]);i++)
	    {
	      fprintf(caa_input,"%d ",lga_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  fprintf(caa_input,"g_a_model=%d\n",lga_g_a_model);
  if(lga_cens_model)
    fprintf(caa_input,"cens:k=%f,m=%f,r=%f,Nlim=%f\n",lga_cens[0],lga_cens[1],lga_cens[2],lga_cens[3]);


  fprintf(caa_input,"wgl_nBoats=%d\n",wgl_nBoats);
  fprintf(caa_input,"wgl_int_ncov=%d,wgl_int_ispat=%d\n",n_cov[5],ispat[5]);

  for(i=0;i<(n_cov[5]);i++)
    fprintf(caa_input,"i=%d,wgl_int_nFac=%d,wgl_int_fix=%d\n",i,wgl_int_nFac[i],wgl_int_fix[i]);
  fprintf(caa_input,"wgl_int_c_cov\n");
  ind = 0;
  for(h=0;h<wgl_nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(n_cov[5]);i++)
	{
	  fprintf(caa_input,"%d ",wgl_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_slp_ncov=%d,wgl_slp_ispat=%d\n",n_cov[6],ispat[6]);
  for(i=0;i<(n_cov[6]);i++)
    fprintf(caa_input,"i=%d,wgl_slp_nFac=%d,wgl_slp_fix=%d\n",i,wgl_slp_nFac[i],wgl_slp_fix[i]);
  fprintf(caa_input,"wgl_slp_c_cov\n");
  ind = 0;
  for(h=0;h<wgl_nBoats;h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(n_cov[6]);i++)
	{
	  fprintf(caa_input,"%d ",wgl_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_hsz_ncov=%d,wgl_hsz_ispat=%d\n",n_cov[7],ispat[7]);
  for(i=0;i<(n_cov[7]);i++)
    fprintf(caa_input,"i=%d,wgl_hsz_nFac=%d,wgl_hsz_fix=%d\n",i,wgl_hsz_nFac[i],wgl_hsz_fix[i]);
  if(n_cov[7]>0)
    {
      fprintf(caa_input,"wgl_hsz_c_cov\n");
      ind = 0;
      for(h=0;h<wgl_nBoats;h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(n_cov[7]);i++)
	    {
	      fprintf(caa_input,"%d ",wgl_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  if(ispat[0] > -1 && age_int_nFac[ispat[0]]>0)
    {
      fprintf(caa_input,"Neighborhood\n"); 
      ind = 0;
      for(i=0;i<(age_int_nFac[ispat[0]-1]);i++)
	 {
	   fprintf(caa_input,"%d ",num_adj_area[i]);
	   for(j=0;j<num_adj_area[i];j++)
	     {
	       fprintf(caa_input,"%d ",adj_area[ind]);
	       ind++;
	     }
	   fprintf(caa_input,"\n");
	 }
    }

  fprintf(caa_input,"tot_nCell=%d, tot_nFactors=%d\n",tot_nCell,tot_nFactors);

  fprintf(caa_input,"tot_fac_age_int:");
  for(i=0;i<n_cov[0];i++)
    fprintf(caa_input,"%d ",tot_fac_age_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_age_hsz:");
  for(i=0;i<n_cov[1];i++)
    fprintf(caa_input,"%d ",tot_fac_age_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"tot_fac_lga_int:");
  for(i=0;i<n_cov[2];i++)
    fprintf(caa_input,"%d ",tot_fac_lga_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_lga_slp:");
  for(i=0;i<n_cov[3];i++)
    fprintf(caa_input,"%d ",tot_fac_lga_slp[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_lga_hsz:");
  for(i=0;i<n_cov[4];i++)
    fprintf(caa_input,"%d ",tot_fac_lga_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"tot_fac_wgl_int:");
  for(i=0;i<n_cov[2];i++)
    fprintf(caa_input,"%d ",tot_fac_wgl_int[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_wgl_slp:");
  for(i=0;i<n_cov[3];i++)
    fprintf(caa_input,"%d ",tot_fac_wgl_slp[i]);
  fprintf(caa_input,"\n");
  fprintf(caa_input,"tot_fac_wgl_hsz:");
  for(i=0;i<n_cov[4];i++)
    fprintf(caa_input,"%d ",tot_fac_wgl_hsz[i]);
  fprintf(caa_input,"\n");

  fprintf(caa_input,"inc_haul= %d %d %d %d %d %d %d %d\n",
	  inc_haul[0],inc_haul[1],inc_haul[2],inc_haul[3],
	  inc_haul[4],inc_haul[5],inc_haul[6],inc_haul[7]);
  fprintf(caa_input,"tot_factors catch\n");
  ind = 0;
  for(c=0;c<tot_nCell;c++)
    {
      fprintf(caa_input,"%d ",c);
      for(i=0;i<tot_nFactors;i++)
	{
	  fprintf(caa_input,"%d ",tot_factors[ind]);
	  ind++;
	}
      fprintf(caa_input,"%f\n",tot_catch[c]);
    }
  if((n_cov[1]+n_cov[4]+n_cov[7])>0)
    {
      fprintf(caa_input,"parameters for haulsize\n");
      for(c=0;c<tot_nCell;c++)
	fprintf(caa_input,"%f %f\n",par_haulsize[2*c],par_haulsize[2*c+1]);
    }
     
  fprintf(caa_input,"N_l_int=%d\n",N_l_int);
  fprintf(caa_input,"l_int\n");
  for(i=0;i<N_l_int;i++)
    fprintf(caa_input,"%f\n",l_int[i]);

  fprintf(caa_input,"nMC=%d\n",nMC);

  // Distribution for cell effects
  fprintf(caa_input,"age_int_cell_E,num_cell=(%d,%d):\n",num_cell_o[0],num_cell_u[0]);
  for(i=0;i<MIN(10,num_cell_u[0]);i++)
    {
      for(j=0;j<MIN(5,num_cell_o[0]);j++)
	fprintf(caa_input,"%f ",age_int_cell_E[i*num_cell_o[0]+j]);
      fprintf(caa_input,"\n");
    }
  fprintf(caa_input,"age_int_cell_C,n=%d;\n",age_int_cell_nC);
  for(i=0;i<MIN(10,num_cell_u[0]);i++)
    {
      for(j=0;j<MIN(5,age_int_cell_nC);j++)
	fprintf(caa_input,"%f ",age_int_cell_C[i*age_int_cell_nC+j]);
      fprintf(caa_input,"\n");
    }
  fclose(caa_input);

  return(0);
}		/* end of caa_write_input_predict */

/*F:*

________________________________________________________________

		
________________________________________________________________

Name:		 write_input_BF
Syntax:	 See call below
Description:    Writes input for caa_Bayes_factor routine to file 
                caa_input_BF.txt in the
                current directory. Only to be used for testing.
Side effects:   A file caa_input_model1.txt is made
Return value:   0
Global or static variables used: None
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_routines.c,v 1.1 2009/06/09 09:27:06 mmerzere Exp $
________________________________________________________________
*/
int write_input_marg_dens(int *i_Nmcmc,int *i_burnin,double *i_mcmc1,double *i_mcmc2,
			  int *i_num_par1,int *i_num_par2,
			  int *nBoats,int *nAges,int *a_vec,
			  int *i_n_cov,int *i_ispat,
			  int *i_age_int_nFac,int *i_age_int_fix,int *i_age_int_c_cov,
			  int *i_age_hsz_nFac,int *i_age_hsz_fix,int *i_age_hsz_c_cov,
			  int *i_lga_nBoats,
			  int *i_lga_int_nFac,int *i_lga_int_fix,int *i_lga_int_c_cov,
			  int *i_lga_slp_nFac,int *i_lga_slp_fix,int *i_lga_slp_c_cov,
			  int *i_lga_hsz_nFac,int *i_lga_hsz_fix,int *i_lga_hsz_c_cov,
			  int *i_lga_g_a_model,
			  int *i_wgl_nBoats,
			  int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
			  int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
			  int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
			  int *i_num_adj_area,int *i_adj_area)
{
  int        a,h,i,j,ind;
  FILE      *caa_input;

  caa_input = fopen("caa_input_predict.txt","w");

  fprintf(caa_input,"Nmcmc = %ld burnin= %ld\n",*i_Nmcmc,*i_burnin);

  fprintf(caa_input,"nBoats_age=%ld,nAges=%ld\n",
	  *nBoats,*nAges);
  for(a=0;a<(*nAges);a++)
    fprintf(caa_input,"i=%d,avec=%ld\n",a,a_vec[a]);

  fprintf(caa_input,"age_int_ncov=%ld,i_age_int_ispat=%ld\n",i_n_cov[0],i_ispat[0]);
  for(i=0;i<(i_n_cov[0]);i++)
    fprintf(caa_input,"i=%d,i_age_int_nFac=%ld,i_age_int_fix=%ld\n",i,i_age_int_nFac[i],i_age_int_fix[i]);
  ind = 0;
  fprintf(caa_input,"age_int_c_cov\n");
  for(h=0;h<(*nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[0]);i++)
	{
	  fprintf(caa_input,"%ld ",i_age_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(i_n_cov[1]>0)
    {
      fprintf(caa_input,"age_hsz_ncov=%ld,i_age_hsz_ispat=%ld\n",i_n_cov[1],i_ispat[1]);
      for(i=0;i<(i_n_cov[1]);i++)
	fprintf(caa_input,"i=%d,i_age_hsz_nFac=%ld,i_age_hsz_fix=%ld\n",i,i_age_hsz_nFac[i],i_age_hsz_fix[i]);
      ind = 0;
      fprintf(caa_input,"age_hsz_c_cov\n");
      for(h=0;h<(*nBoats);h++)
	{
	  fprintf(caa_input,"%d ",h);
	  for(i=0;i<(i_n_cov[1]);i++)
	    {
	      fprintf(caa_input,"%ld ",i_age_hsz_c_cov[ind]);
	      ind++;
	    }
	  fprintf(caa_input,"\n");
	}
    }

  fprintf(caa_input,"i_lga_nBoats=%ld\n",*i_lga_nBoats);
  fprintf(caa_input,"lga_int_ncov=%ld,lga_int_ispat=%ld\n",i_n_cov[2],i_ispat[2]);
  for(i=0;i<(i_n_cov[2]);i++)
    fprintf(caa_input,"i=%d,lga_int_nFac=%ld,lga_int_fix=%ld\n",i,i_lga_int_nFac[i],i_lga_int_fix[i]);
  fprintf(caa_input,"lga_int_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[2]);i++)
	{
	  fprintf(caa_input,"%ld ",i_lga_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_slp_ncov=%ld,lga_slp_ispat=%ld\n",i_n_cov[3],i_ispat[3]);
  for(i=0;i<(i_n_cov[3]);i++)
    fprintf(caa_input,"i=%d,lga_slp_nFac=%ld,lga_slp_fix=%ld\n",i,i_lga_slp_nFac[i],i_lga_slp_fix[i]);
  fprintf(caa_input,"lga_slp_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[3]);i++)
	{
	  fprintf(caa_input,"%ld ",i_lga_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"lga_hsz_ncov=%ld,lga_hsz_ispat=%ld\n",i_n_cov[4],i_ispat[4]);
  for(i=0;i<(i_n_cov[4]);i++)
    fprintf(caa_input,"i=%d,lga_hsz_nFac=%ld,lga_hsz_fix=%ld\n",i,i_lga_hsz_nFac[i],i_lga_hsz_fix[i]);
  fprintf(caa_input,"lga_hsz_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_lga_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[4]);i++)
	{
	  fprintf(caa_input,"%ld ",i_lga_hsz_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"i_wgl_nBoats=%ld\n",*i_wgl_nBoats);
  fprintf(caa_input,"wgl_int_ncov=%ld,wgl_int_ispat=%ld\n",i_n_cov[5],i_ispat[5]);
  for(i=0;i<(i_n_cov[5]);i++)
    fprintf(caa_input,"i=%d,wgl_int_nFac=%ld,wgl_int_fix=%ld\n",i,i_wgl_int_nFac[i],i_wgl_int_fix[i]);
  fprintf(caa_input,"wgl_int_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[5]);i++)
	{
	  fprintf(caa_input,"%ld ",i_wgl_int_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_slp_ncov=%ld,wgl_slp_ispat=%ld\n",i_n_cov[6],i_ispat[6]);
  for(i=0;i<(i_n_cov[6]);i++)
    fprintf(caa_input,"i=%d,wgl_slp_nFac=%ld,wgl_slp_fix=%ld\n",i,i_wgl_slp_nFac[i],i_wgl_slp_fix[i]);
  fprintf(caa_input,"wgl_slp_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[6]);i++)
	{
	  fprintf(caa_input,"%ld ",i_wgl_slp_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  fprintf(caa_input,"wgl_hsz_ncov=%ld,wgl_hsz_ispat=%ld\n",i_n_cov[7],i_ispat[7]);
  for(i=0;i<(i_n_cov[7]);i++)
    fprintf(caa_input,"i=%d,wgl_hsz_nFac=%ld,wgl_hsz_fix=%ld\n",i,i_wgl_hsz_nFac[i],i_wgl_hsz_fix[i]);
  fprintf(caa_input,"wgl_hsz_c_cov\n");
  ind = 0;
  for(h=0;h<(*i_wgl_nBoats);h++)
    {
      fprintf(caa_input,"%d ",h);
      for(i=0;i<(i_n_cov[7]);i++)
	{
	  fprintf(caa_input,"%ld ",i_wgl_hsz_c_cov[ind]);
	  ind++;
	}
      fprintf(caa_input,"\n");
    }

  if(i_age_int_nFac[i_ispat[0]]>0)
    {
      fprintf(caa_input,"Neighborhood\n"); 
      ind = 0;
      for(i=0;i<(i_age_int_nFac[i_ispat[0]-1]);i++)
	 {
	   fprintf(caa_input,"%ld ",i_num_adj_area[i]);
	   for(j=0;j<i_num_adj_area[i];j++)
	     {
	       fprintf(caa_input,"%ld ",i_adj_area[ind]);
	       ind++;
	     }
	   fprintf(caa_input,"\n");
	 }
    }

  return(0);
}		/* end of caa_write_input_marg_dens */

    

/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing original data

  Makes a struct of type Data_orig (see caa.h for definition)

  Space allocated in this routine is reallocated in re_makedata_orig
*/
int makedata_orig(SEXP i_dataList, Data_orig **o_D_orig)
{
  Data_orig  *D_orig;
  int         i;
  SEXP        elmt = R_NilValue;

  D_orig = CALLOC(1,Data_orig);     // Free ok

  /* Picking out data vectors from input list */
  /* NOTE! Variables are read-only!! */

  if(!Rf_isNull(elmt = getListElement(i_dataList, "totage")))
    D_orig->totage = INTEGER_POINTER(AS_INTEGER(elmt)); // vector of ages
 
  if(!Rf_isNull(elmt = getListElement(i_dataList, "totlength")))
    D_orig->totlength = NUMERIC_POINTER(elmt); // vector of lengths

  if(!Rf_isNull(elmt = getListElement(i_dataList, "haulweight")))
    D_orig->haulweight = NUMERIC_POINTER(elmt); // vector of haul sizes

  if(!Rf_isNull(elmt = getListElement(i_dataList, "totseason")))
    D_orig->season = INTEGER_POINTER(AS_INTEGER(elmt)); // vector of seasons or months

  D_orig->coastal_cod = INTEGER_VALUE(getListElement(i_dataList, "coastal_cod"));
  if(D_orig->coastal_cod)
    {
      if(!Rf_isNull(elmt = getListElement(i_dataList, "tottype")))
	D_orig->tottype = INTEGER_POINTER(AS_INTEGER(elmt)); // vector of type (coastal/atlantic cod)
    }

  if(!Rf_isNull(elmt = getListElement(i_dataList, "nFishBoat")))
    D_orig->nFishBoat = INTEGER_POINTER(AS_INTEGER(elmt)); // number of fish pr boat

  if(!Rf_isNull(elmt = getListElement(i_dataList, "replength")))
    D_orig->replength = INTEGER_POINTER(AS_INTEGER(elmt)); // repetitions of lengths

  if(!Rf_isNull(elmt = getListElement(i_dataList, "start_noAge")))
    D_orig->start_noAge = INTEGER_POINTER(AS_INTEGER(elmt)); // index at which no-aged fish start for each haul

  if(!Rf_isNull(elmt = getListElement(i_dataList, "num_noAge")))
    D_orig->num_noAge = INTEGER_POINTER(AS_INTEGER(elmt)); // number of no-aged fish for each haul

  if(!Rf_isNull(elmt = getListElement(i_dataList, "n_int_len")))
    D_orig->n_int_len = INTEGER_VALUE(elmt); // number of intervals for length

  if(!Rf_isNull(elmt = getListElement(i_dataList, "int_len_lim")))
    D_orig->int_len_lim = NUMERIC_POINTER(elmt); // lower limits of length-intervals

  if(!Rf_isNull(elmt = getListElement(i_dataList, "int_len_vec")))
    D_orig->int_len = NUMERIC_POINTER(AS_NUMERIC(elmt)); // length value for intervals


  *o_D_orig = D_orig;

  return(0);
}		/* end of makedata_orig */

    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in makedata_orig
*/
int re_makedata_orig(Data_orig **o_D_orig)
{
  Data_orig   *D_orig;

  D_orig = *o_D_orig;

  FREE(D_orig);
  
  return(0);
}		/* end of re_makedata_orig */

    

/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing parameters for coastal cod

  Makes a struct of type Data_CC (see caa.h for definition)

  Space allocated in this routine is reallocated in re_makedata_CC
*/
int makedata_CC(SEXP i_ageList, Data_CC **o_D_CC)
{
  Data_CC  *D_CC;
  int       i;
  SEXP      elmt = R_NilValue;

  D_CC = CALLOC(1,Data_CC);    

  /* Picking out data vectors from input list */
  /* NOTE! Variables are read-only!! */

  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype1.CC")))
    D_CC->ptype1_CC = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype1.S")))
    D_CC->ptype1_S = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype2.CC")))
    D_CC->ptype2_CC = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype2.S")))
    D_CC->ptype2_S = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype4.CC")))
    D_CC->ptype4_CC = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype4.S")))
    D_CC->ptype4_S = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt = getListElement(i_ageList, "ptype5.CC")))
    D_CC->ptype5_CC = NUMERIC_POINTER(elmt);
  if(!Rf_isNull(elmt  = getListElement(i_ageList, "ptype5.S")))
    D_CC->ptype5_S = NUMERIC_POINTER(elmt);

  *o_D_CC = D_CC;

  return(0);
}		/* end of makedata_CC */

    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in makedata_CC
*/
int re_makedata_CC(Data_CC **o_D_CC)
{
  Data_CC   *D_CC;

  D_CC = *o_D_CC;

  FREE(D_CC);
  
  return(0);
}		/* end of re_makedata_CC */

    

/*!
  \author Geir Storvik
  \brief Makes a struct of type containing model and covariates for age model

  Makes a struct of type Data_age (see caa.h for definition)
  describing the model structure, covariates and number of hauls
  with and without age-data.

  Space allocated in this routine is reallocated in re_makedata_age1 
*/
int makedata_age1(int i_nBoats,int i_nAges,int *i_a_vec,
		  int i_int_n_cov,int *i_int_nFac,int i_int_ispat,
		  int *i_int_fix,int *i_int_x_cov,
		  int *i_int_num,int *i_int_adj_area,
		  int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,
		  int *i_hsz_fix,int *i_hsz_x_cov,
		  int *i_hsz_num,int *i_hsz_adj_area,
		  Data_age **o_D_age)
{
  int         a,err;
  Data_cov   *xcov;
  Data_age   *D_age;

  D_age = CALLOC(1,Data_age);       // Free ok
  D_age->glm = CALLOC(1,Data_glm);  // Free ok

  /* Sizes */
  D_age->glm->nHaul = i_nBoats;
  D_age->glm->ncat=i_nAges;

  D_age->a_vec = CALLOC(D_age->glm->ncat,int);  // Free ok
  for(a=0;a<D_age->glm->ncat;a++)
    D_age->a_vec[a] = i_a_vec[a];

  /* Covariates */
  if(i_hsz_n_cov>0)
    D_age->glm->nxcov = 2;  /* Intercept and haulsize */
  else
    D_age->glm->nxcov = 1;  /* Only Intercept */
  D_age->glm->xcov = CALLOC(D_age->glm->nxcov,Data_cov *);   // Free ok

  // Intercept
  make_c_cov(i_int_n_cov,i_int_ispat,i_int_fix,i_int_x_cov,D_age->glm->nHaul,&xcov);
  D_age->glm->xcov[0] = xcov;

  /* Convert covariates */
  err = convert_cov(D_age->glm->nHaul,i_int_nFac,D_age->glm->xcov[0]);
  if(err)
    {
      write_warning("makedata_age1:Error calling convert_cov\n");
      return(err);
    }

  /* spatial structure */
  if(D_age->glm->xcov[0]->ispat > -1)
      make_spat_struct(i_int_num,i_int_adj_area,xcov);

  // Haulsize
  if(i_hsz_n_cov>0)
    {
      make_c_cov(i_hsz_n_cov,i_hsz_ispat,i_hsz_fix,i_hsz_x_cov,D_age->glm->nHaul,&xcov);
      D_age->glm->xcov[1] = xcov;

      /* Convert covariates */
      err = convert_cov(D_age->glm->nHaul,i_hsz_nFac,D_age->glm->xcov[1]);
      if(err)
	{
	  write_warning("makedata_age1:Error calling convert_cov\n");
	  return(err);
	}

      /* spatial structure */
      if(D_age->glm->xcov[1]->ispat > -1)
	make_spat_struct(i_hsz_num,i_hsz_adj_area,xcov);
    }

  *o_D_age = D_age;

  return(0);
}		/* end of makedata_age1 */

    

/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_age1
*/
int re_makedata_age1(int i_nBoats,int i_nAges,int *i_a_vec,
		     int i_int_n_cov,int *i_int_nFac,int i_int_ispat,
		     int *i_int_fix,int *i_int_x_cov,
		     int *i_int_num,int *i_int_adj_area,
		     int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,
		     int *i_hsz_fix,int *i_hsz_x_cov,
		     int *i_hsz_num,int *i_hsz_adj_area,
		     Data_age **o_D_age)
{
  int         err;
  Data_cov   *xcov;
  Data_age   *D_age;

  D_age = *o_D_age;

  if(i_hsz_n_cov>0)
    {
      xcov = D_age->glm->xcov[1];

      if(D_age->glm->xcov[1]->ispat > -1)
	re_make_spat_struct(i_hsz_num,i_hsz_adj_area,xcov);

      err = re_convert_cov(D_age->glm->nHaul,i_hsz_nFac,D_age->glm->xcov[1]);
      if(err)
	{
	  write_warning("re_makedata_age1:Error calling re_convert_cov\n");
	  return(err);
	}

      re_make_c_cov(i_hsz_n_cov,i_hsz_ispat,i_hsz_fix,i_hsz_x_cov,D_age->glm->nHaul,&xcov);
    }

  xcov = D_age->glm->xcov[0];

  if(D_age->glm->xcov[0]->ispat > -1)
      re_make_spat_struct(i_int_num,i_int_adj_area,xcov);

  err = re_convert_cov(D_age->glm->nHaul,i_int_nFac,D_age->glm->xcov[0]);
  if(err)
    {
      write_warning("re_makedata_age1:Error calling re_convert_cov\n");
      return(err);
    }

  re_make_c_cov(i_int_n_cov,i_int_ispat,i_int_fix,i_int_x_cov,D_age->glm->nHaul,&xcov);

  FREE(D_age->glm->xcov);

  FREE(D_age->glm);

  FREE(D_age->a_vec);

  FREE(D_age);

  return(0);
}		/* end of re_makedata_age1 */

    

/*!
  \author Geir Storvik
  \brief Put observed ages for different hauls into D_age->Ages for Amigo data.

  Note that the first hauls are for length-only data so that
  Memory allocated in this routine is reallocated in re_makedata_age2
*/
int makedata_age2(Data_orig *i_D_orig,Data_age *i_D_age)
{
  int         a,a2,f,h,n,ind;

  /* Ages */
  i_D_age->n_h = CALLOC(i_D_age->glm->nHaul,int);    // Free ok
 
  i_D_age->Ages = CALLOC2_i(i_D_age->glm->nHaul,i_D_age->glm->ncat);

  i_D_age->Ages_fix = CALLOC2_i(i_D_age->glm->nHaul,i_D_age->glm->ncat);

  i_D_age->Ages_disc = Mmatrix_2d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat,
				   sizeof(int),1); //Free ok
  i_D_age->Ages_land = Mmatrix_2d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat,
				   sizeof(int),1); //Free ok

  if(i_D_orig->coastal_cod) //difference between coastal cod and skrei
    {
      a2 = (int) i_D_age->glm->ncat/2;
      ind = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages[h][a] = 0;
	  for(f=0;f<i_D_orig->nFishBoat[h];f++)
	    {
	      if(i_D_orig->totage[ind] > -1000)
		{
		  a = i_D_orig->totage[ind]-i_D_age->a_vec[0];
		  if(a < 0 || a >= i_D_age->a_vec[i_D_age->glm->ncat-1])
		    {
		      fprintf(stderr,"h=%d,f=%d,ind=%d,totage=%ld,a_vec0=%d\n",
			      h,f,ind,i_D_orig->totage[ind],i_D_age->a_vec[0]);
		      write_warning("makedata_age2:Something is wrong\n");
		    }
		  /* Use certain observations for starting values, 
		     even when classification error */
		  if(i_D_orig->tottype[ind] == 1) //certain coastal cod
		    i_D_age->Ages[h][a] += i_D_orig->replength[ind];
		  else if(i_D_orig->tottype[ind] == 5) //certain skrei
		    i_D_age->Ages[h][a2+a] += i_D_orig->replength[ind];
		}
	      ind++;
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages_fix[h][a] = i_D_age->Ages[h][a];
	  n = 0;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    n += i_D_age->Ages[h][a];
	  i_D_age->n_h[h] = n;
	}
    }
  else  //original version
    {
      ind = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages[h][a] = 0;
	  for(f=0;f<i_D_orig->nFishBoat[h];f++)
	    {
	      if(i_D_orig->totage[ind] > -1000)
		{
		  a = i_D_orig->totage[ind]-i_D_age->a_vec[0];
		  if(a < 0 || a >= i_D_age->glm->ncat)
		    {
		      fprintf(stderr,"h=%d,f=%d,ind=%d,totage=%ld,a_vec0=%d\n",
			      h,f,ind,i_D_orig->totage[ind],i_D_age->a_vec[0]);
		      write_warning("makedata_age2:Something is wrong\n");
		    }
		  i_D_age->Ages[h][i_D_orig->totage[ind]-i_D_age->a_vec[0]] += i_D_orig->replength[ind];
		}
	      ind++;
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages_fix[h][a] = i_D_age->Ages[h][a];
	  n = 0;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    n += i_D_age->Ages[h][a];
	  i_D_age->n_h[h] = n;
	}
    }

  return(0);
}		/* end of makedata_age2 */

    

/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_age2
*/
int re_makedata_age2(Data_age *i_D_age)
{
  /* Ages */
  FREE(i_D_age->n_h);
  FREE2_i(i_D_age->Ages,i_D_age->glm->nHaul);
  FREE2_i(i_D_age->Ages_fix,i_D_age->glm->nHaul);
  // Fmatrix_2d(&i_D_age->Ages[0][0],&i_D_age->Ages[0]);
  Fmatrix_2d(&i_D_age->Ages_disc[0][0],&i_D_age->Ages_disc[0]);
  Fmatrix_2d(&i_D_age->Ages_land[0][0],&i_D_age->Ages_land[0]);

  return(0);
}		/* end of re_makedata_age2 */



/*!
  \author Hanne Rognebakke
  \brief Allocate space and initialize struct for g_a-data
*/
int makedata_g_a(int lga_g_a_ncat,int lga_g_a_nSeason,double *i_avec,int *i_a2Age_vec,
		 int i_g_a_model,Data_age *i_D_age,Data_g_a **o_g_a)
{
  int       err,a,A;
  Data_g_a  *g_a;
 

  /* Allocating space */
  g_a = CALLOC(1,Data_g_a);                  // Free ok
  g_a->ncat = lga_g_a_ncat;
  g_a->nSeason = lga_g_a_nSeason;
  g_a->a_vec = CALLOC(g_a->ncat,double);  // Free ok
  for(a=0;a<g_a->ncat;a++)
    g_a->a_vec[a] = i_avec[a];

  g_a->a2Age_vec = CALLOC(i_D_age->glm->ncat,int);  // Free ok
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      g_a->a2Age_vec[a] = i_a2Age_vec[a];
      //printf("g_a_a2Age_vec[%d]=%d\n",a,g_a->a2Age_vec[a]);
    }

  g_a->g_a_model = i_g_a_model;
  /* if g_a->g_a_model==0  log-linear model */
  if(g_a->g_a_model==1) /* Schnute-Richards model */
    {
      g_a->g_a_npar = 3;
      g_a->g_a_par = CALLOC(3,double);            // Free ok
    }
  else if(g_a->g_a_model==2) /* polynomial model, not implemented yet */
    {
      g_a->g_a_npar = 2;
      g_a->g_a_par = CALLOC(2,double);            // Free ok
    }
  else if(g_a->g_a_model!=0)
    {
      write_warning("makedata_g_a:Unknown g_a_model\n");
      return(1);
    }

  g_a->g_a = CALLOC(g_a->ncat,double);   // Free ok
  /* Initialization of nonlinear function between log-length and age */
  A = g_a->ncat-1;
  for(a=0;a<g_a->ncat;a++)
    {
      g_a->g_a[a] = 
	(log(g_a->a_vec[a])-log(g_a->a_vec[0]))/
	(log(g_a->a_vec[A])-log(g_a->a_vec[0]));
      //printf("g_a[%d]=%f\n",a,g_a->g_a[a]);
    }

  /* Allocating space for sufficient statistics for g(a) */
  g_a->suff = Mmatrix_2d(0,1,0,g_a->ncat-1,sizeof(double),1);  // Free ok

  *o_g_a = g_a;

  return(0);
}		/* end of makedata_g_a */



/*!
  \author Hanne Rognebakke
  \brief Re-allocating space allocated by makedata_g_a
*/
int re_makedata_g_a(Data_g_a **o_g_a)
{
  int       err;

  Data_g_a  *g_a;

  g_a = *o_g_a;
  FREE(g_a->a_vec);
  FREE(g_a->a2Age_vec);
  FREE(g_a->g_a);
  if(g_a->g_a_model>0)
      FREE(g_a->g_a_par);
  Fmatrix_2d(&g_a->suff[0][0],&g_a->suff[0]);

  FREE(g_a);

  return(0);
}		/* end of re_makedata_g_a */



/*!
  \author Geir Storvik
  \brief Construct a struct describing the model structure and covariates for linear models

  With and without covariates.
  Space allocated in this routine is reallocated in re_makedata_lin1
  Note that effects are transformed to start at zero.

  Memory allocated in this routines is reallocated in re_makedata_lin1
*/
int makedata_lin1(int i_nBoats,
		  int i_int_n_cov,int *i_int_nFac,int i_int_ispat,
		  int *i_int_fix,int *i_int_x_cov,
		  int *i_int_num,int *i_int_adj_area,
		  int i_slp_n_cov,int *i_slp_nFac,int i_slp_ispat,
		  int *i_slp_fix,int *i_slp_x_cov,
		  int *i_slp_num,int *i_slp_adj_area,
		  int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,
		  int *i_hsz_fix,int *i_hsz_x_cov,
		  int *i_hsz_num,int *i_hsz_adj_area,
		  Data_lin **o_D_lin)
{
  int        err;
  Data_cov  *xcov;
  Data_lin  *D_lin;

  D_lin = CALLOC(1,Data_lin);         // Free ok
  D_lin->glm = CALLOC(1,Data_glm);    // Free ok
  D_lin->glm->ncat = 1;

  /* Sizes */
  D_lin->glm->nHaul = i_nBoats;

  /* Covariates */
  if(i_hsz_n_cov>0)
    D_lin->glm->nxcov = 3;  /* Intercept and slope and haulsize */
  else
    D_lin->glm->nxcov = 2;  /* Intercept and slope */
  D_lin->glm->xcov = CALLOC(D_lin->glm->nxcov,Data_cov *);   // Free ok


  /* Intercept */
  err = make_c_cov(i_int_n_cov,i_int_ispat,i_int_fix,i_int_x_cov,D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("makedata_lin1:Error calling make_c_cov for intercept\n");
      return(err);
    }
  D_lin->glm->xcov[0] = xcov;

  /* Convert covariates */
  err = convert_cov(D_lin->glm->nHaul,i_int_nFac,D_lin->glm->xcov[0]);
  if(err)
    {
      write_warning("makedata_lin1:Error calling convert_cov for intercept\n");
      return(err);
    }

  /* Make spatial structure */

  if(D_lin->glm->xcov[0]->ispat > -1)
    {
      err = make_spat_struct(i_int_num,i_int_adj_area,D_lin->glm->xcov[0]);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling make_spat_struct\n");
	  return(err);
	}
    }
  
  /* Slope */
  err = make_c_cov(i_slp_n_cov,i_slp_ispat,i_slp_fix,i_slp_x_cov,D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("makedata_lin1:Error calling make_c_cov\n");
      return(err);
    }
  D_lin->glm->xcov[1] = xcov;

  /* Convert covariates */
  err = convert_cov(D_lin->glm->nHaul,i_slp_nFac,D_lin->glm->xcov[1]);
  if(err)
    {
      write_warning("makedata_lin1:Error calling convert_cov for slope\n");
      return(err);
    }

  /* Make spatial structure */
  if(D_lin->glm->xcov[1]->ispat > -1)
    {
      err = make_spat_struct(i_slp_num,i_slp_adj_area,D_lin->glm->xcov[1]);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling make_spat_struct\n");
	  return(err);
	}
    }

  /* haulsize */

  if(i_hsz_n_cov>0)
    {
      err = make_c_cov(i_hsz_n_cov,i_hsz_ispat,i_hsz_fix,i_hsz_x_cov,D_lin->glm->nHaul,&xcov);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling make_c_cov for haulsize\n");
	  return(err);
	}
      D_lin->glm->xcov[2] = xcov;

      /* Convert covariates */
      err = convert_cov(D_lin->glm->nHaul,i_hsz_nFac,D_lin->glm->xcov[2]);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling convert_cov\n");
	  return(err);
	}

      /* Make spatial structure */
      if(D_lin->glm->xcov[2]->ispat > -1)
	{
	  err = make_spat_struct(i_hsz_num,i_hsz_adj_area,D_lin->glm->xcov[1]);
	  if(err)
	    {
	      write_warning("makedata_lin1:Error calling make_spat_struct\n");
	      return(err);
	    }
	}
    }


  *o_D_lin = D_lin;

  return(0);
}		/* end of makedata_lin1 */



/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_lin1
*/
int re_makedata_lin1(int i_nBoats,
		     int i_int_n_cov,int *i_int_nFac,int i_int_ispat,int *i_int_fix,int *i_int_x_cov,
		     int *i_int_num,int *i_int_adj_area,
		     int i_slp_n_cov,int *i_slp_nFac,int i_slp_ispat,int *i_slp_fix,int *i_slp_x_cov,
		     int *i_slp_num,int *i_slp_adj_area,
		     int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,int *i_hsz_fix,int *i_hsz_x_cov,
		     int *i_hsz_num,int *i_hsz_adj_area,
		     Data_lin **o_D_lin)
{
  int        err;
  Data_cov  *xcov;
  Data_lin  *D_lin;

  D_lin = *o_D_lin;

  /* Int */
  xcov = D_lin->glm->xcov[0];
  if(xcov->ispat > -1)
    {
      err = re_make_spat_struct(i_int_num,i_int_adj_area,xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_make_spat_struct\n");
	  return(err);
	}
    }
  err = re_convert_cov(D_lin->glm->nHaul,i_int_nFac,xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_convert_cov\n");
      return(err);
    }
  err = re_make_c_cov(i_int_n_cov,i_int_ispat,i_int_fix,i_int_x_cov,D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_make_c_cov\n");
      return(err);
    }


  /* Slope */
  xcov = D_lin->glm->xcov[1];
  if(xcov->ispat > -1)
    {
      err = re_make_spat_struct(i_slp_num,i_slp_adj_area,xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_ make_spat_struct\n");
	  return(err);
	}
    }
  err = re_convert_cov(D_lin->glm->nHaul,i_slp_nFac,xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_convert_cov\n");
      return(err);
    }
  err = re_make_c_cov(i_slp_n_cov,i_slp_ispat,i_slp_fix,i_slp_x_cov,D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_make_c_cov\n");
      return(err);
    }

  if(i_hsz_n_cov>0)
    {
      xcov = D_lin->glm->xcov[2];
      if(xcov->ispat > -1)
	{
	  err = re_make_spat_struct(i_hsz_num,i_hsz_adj_area,xcov);
	  if(err)
	    {
	      write_warning("re_makedata_lin1:Error calling re_ make_spat_struct\n");
	      return(err);
	    }
	}
      err = re_convert_cov(D_lin->glm->nHaul,i_hsz_nFac,xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_convert_cov\n");
	  return(err);
	}
      err = re_make_c_cov(i_hsz_n_cov,i_hsz_ispat,i_hsz_fix,i_hsz_x_cov,D_lin->glm->nHaul,&xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_make_c_cov\n");
	  return(err);
	}
    }

  FREE(D_lin->glm->xcov);
 
  FREE(D_lin->glm);

  FREE(D_lin);

  return(0);
}		/* end of re_makedata_lin1 */



/*!
  \author Geir Storvik
  \brief Calculates sufficient statistics for the lga model 

  Use amigo data and age-stratified by length data.

  This routine is only used in initialization. 
  It calculates summary statistics based on aged fish.
  Updating sufficient statistics based on length-only data is performed
  by the make_suff_lga routine.

  \todo The routine copy summaries to the fix-versions. This should be made
  more clean later on.

  Memory allocated in this routine is reallocated in re_makedata_lga_suff.
*/
int makedata_lga_suff(Data_lin *i_D_lga,Data_orig *i_D_orig,Data_g_a *i_D_g_a)
{
  int        a,f,h,season,err,n,ncat,nmiss,N,ind_f,ind_a;
  double    *x,*y,beta0,beta1,ssq,length,r;


  ncat = i_D_g_a->ncat;

  /* Sufficient Statistics */
  i_D_lga->Ages = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				   0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				     0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->Ages_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				       0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
					 0,ncat-1,sizeof(double),1);    // Free ok
  /* Need size of beta_hat equal to nxcov when finding b-vector in sample_gauss_eff */
  i_D_lga->glm->beta_hat = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
				      0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,
				  0,i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok


  N = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);
  x = CALLOC(N,double);                  // Free ok
  y = CALLOC(N,double);                  // Free ok

  ind_f = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<ncat;a++)
	{ 
	  i_D_lga->Ages[h][a] = 0;
	  i_D_lga->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	}
      i_D_lga->glm->suff[h][0][1] = G_ZERO;
      i_D_lga->glm->suff[h][1][1] = G_ZERO;
      n = 0;
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
          a = i_D_orig->totage[ind_f];
	  season = i_D_orig->season[h];
          length = i_D_orig->totlength[ind_f];
	  if(a > -9000 && length > -9000)
	    {
	      ind_a = i_D_g_a->a2Age_vec[a-(int) i_D_g_a->a_vec[0]]+(season-1);
              x[n] = i_D_g_a->g_a[ind_a];
              y[n] = length;
	      r = (double) i_D_orig->replength[ind_f];
              i_D_lga->Ages[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga->sum_by_cat[h][ind_a] += r*y[n];
	      i_D_lga->sqsum_by_cat[h][ind_a] += r*y[n]*y[n];
	      i_D_lga->glm->suff[h][0][0] += r;
	      i_D_lga->glm->suff[h][0][1] += r*x[n];
	      i_D_lga->glm->suff[h][1][1] += r*x[n]*x[n];
	      n++;
	    }
          ind_f ++;
	}
      i_D_lga->glm->suff[h][1][0] = i_D_lga->glm->suff[h][0][1];

      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat,i_D_g_a->g_a,i_D_lga->sum_by_cat[h],i_D_lga->sqsum_by_cat[h],
                            i_D_lga->glm->suff[h],&beta0,&beta1,&ssq);
	  i_D_lga->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga->glm->ssq[h] = ssq;
	}

      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = i_D_lga->glm->suff[h][0][0] * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff[h][0][1] * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][2] = i_D_lga->glm->suff[h][0][0] * i_D_lga->haulweight[h]*i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff[h][0][2];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff[h][1][2];
	  i_D_lga->glm->beta_hat[h][0][2] = G_ZERO; /* needed in find b-vector in sample_gauss_eff() */
	}
      // Copying summaries to the fix-versions. 
      // This should probably be changed to put all directly into the fix-versions
      for(a=0;a<ncat;a++)
	{
	  i_D_lga->Ages_fix[h][a] = i_D_lga->Ages[h][a];
	  i_D_lga->sum_by_cat_fix[h][a] = i_D_lga->sum_by_cat[h][a];
	  i_D_lga->sqsum_by_cat_fix[h][a] = i_D_lga->sqsum_by_cat[h][a];
	}

    } 


  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_lga_suff */



/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_lga_suff
*/
/*F:re_makedata_lga_suff*

________________________________________________________________

		re_makedata_lga_suff
________________________________________________________________

Name:		re_makedata_lga_suff
Syntax:		
Description:    Reallocates memory allocated in makedata_lga_suff
Side effects:
Return value:
Global or static variables used:
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: $Id: caa_routines.c,v 1.1 2009/06/09 09:27:06 mmerzere Exp $
________________________________________________________________

*/
int re_makedata_lga_suff(Data_lin *i_D_lga)
{
  Fmatrix_2d(&i_D_lga->Ages[0][0],&i_D_lga->Ages[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat[0][0],&i_D_lga->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat[0][0],&i_D_lga->sqsum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->Ages_fix[0][0],&i_D_lga->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat_fix[0][0],&i_D_lga->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat_fix[0][0],&i_D_lga->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat[0][0][0],&i_D_lga->glm->beta_hat[0][0],&i_D_lga->glm->beta_hat[0]);
  FREE(i_D_lga->glm->ssq);
  Fmatrix_3d(&i_D_lga->glm->suff[0][0][0],&i_D_lga->glm->suff[0][0],&i_D_lga->glm->suff[0]);
  
  return(0);
}		/* end of re_makedata_lga_suff */



/*!
  \author Hanne Rognebakke
  \brief Calculates sufficient statistics for the lga model when coastal cod 

  Use amigo data and age-stratified by length data.

  This routine is only used in initialization. 
  It calculates summary statistics based on aged fish.
  Updating sufficient statistics based on length-only data is performed
  by the make_suff_lga routine.

  \todo The routine copy summaries to the fix-versions. This should be made
  more clean later on.

  Memory allocated in this routine is reallocated in re_makedata_lga_suff.
*/
int makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC,Data_orig *i_D_orig,
			 Data_g_a *i_D_g_a,Data_g_a *i_D_g_a_CC)
{
  int        ncat,ncat_CC,a,f,h,season,err,N,ind_f,ind_a;
  double     beta0,beta1,ssq,length,r;


  ncat = i_D_g_a->ncat;
  ncat_CC = i_D_g_a_CC->ncat;

  /* Sufficient Statistics */
  i_D_lga->Ages = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat = 
    Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat = 
    Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->Ages_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul,
				 0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat_fix = 
    Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat_fix = 
    Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  /* Need size of beta_hat equal to nxcov when finding b-vector in sample_gauss_eff */
  i_D_lga->glm->beta_hat =
    Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
	       0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff = 
    Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,0,
	       i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok
  
  i_D_lga_CC->Ages = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul,
				0,ncat_CC-1,sizeof(int),1); // Free ok
  i_D_lga_CC->sum_by_cat = 
    Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->sqsum_by_cat = 
    Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->Ages_fix = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul,
				    0,ncat_CC-1,sizeof(int),1); // Free ok
  i_D_lga_CC->sum_by_cat_fix = 
    Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->sqsum_by_cat_fix = 
    Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  /* Need size of beta_hat equal to nxcov when finding b-vector in sample_gauss_eff */
  i_D_lga_CC->glm->beta_hat =
    Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,0,
	       0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->glm->ssq = CALLOC(i_D_lga_CC->glm->nHaul,double);             // Free ok
  i_D_lga_CC->glm->suff = 
    Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,i_D_lga_CC->glm->nxcov-1,
	       0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);
  ind_f = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<ncat;a++)
	{
	  i_D_lga->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	}
      for(a=0;a<ncat_CC;a++)
	{
	  i_D_lga_CC->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga_CC->sqsum_by_cat[h][a] = G_ZERO;
	}
      i_D_lga->glm->suff[h][0][1] = G_ZERO;
      i_D_lga->glm->suff[h][1][1] = G_ZERO;
      i_D_lga_CC->glm->suff[h][0][1] = G_ZERO;
      i_D_lga_CC->glm->suff[h][1][1] = G_ZERO;
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  a = i_D_orig->totage[ind_f];
	  season = i_D_orig->season[h];
	  length = i_D_orig->totlength[ind_f];
	  if(a > -1000 && length > -1000)
	    {
	      /* Use certain observations for starting values even if classification error */
	      if(i_D_orig->tottype[ind_f]==1) //certain coastal cod
		{
		  ind_a = i_D_g_a_CC->a2Age_vec[a-(int) i_D_g_a_CC->a_vec[0]]+(season-1);
		  r = (double) i_D_orig->replength[ind_f];
		  i_D_lga_CC->Ages[h][ind_a] += r;
		  i_D_lga_CC->sum_by_cat[h][ind_a] += r*length;
		  i_D_lga_CC->sqsum_by_cat[h][ind_a] += r*length*length;
		  i_D_lga_CC->glm->suff[h][0][0] += r;
		  i_D_lga_CC->glm->suff[h][0][1] += r*i_D_g_a_CC->g_a[ind_a];
		  i_D_lga_CC->glm->suff[h][1][1] += r*i_D_g_a_CC->g_a[ind_a]*i_D_g_a_CC->g_a[ind_a];
		}
	      else if(i_D_orig->tottype[ind_f]==5) //certain skrei
		{
		  ind_a = i_D_g_a->a2Age_vec[a-(int) i_D_g_a->a_vec[0]]+(season-1);
		  r = (double) i_D_orig->replength[ind_f];
		  i_D_lga->Ages[h][ind_a] += r;
		  i_D_lga->sum_by_cat[h][ind_a] += r*length;
		  i_D_lga->sqsum_by_cat[h][ind_a] += r*length*length;
		  i_D_lga->glm->suff[h][0][0] += r;
		  i_D_lga->glm->suff[h][0][1] += r*i_D_g_a->g_a[ind_a];
		  i_D_lga->glm->suff[h][1][1] += r*i_D_g_a->g_a[ind_a]*i_D_g_a->g_a[ind_a];
		}
	    }
	  ind_f ++;
	}
      i_D_lga->glm->suff[h][1][0] = i_D_lga->glm->suff[h][0][1];
      i_D_lga_CC->glm->suff[h][1][0] = i_D_lga_CC->glm->suff[h][0][1];
      
      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat,i_D_g_a->g_a,i_D_lga->sum_by_cat[h],i_D_lga->sqsum_by_cat[h],
			    i_D_lga->glm->suff[h],&beta0,&beta1,&ssq);
	  i_D_lga->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga->glm->ssq[h] = ssq;
	}
      if(i_D_lga_CC->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat_CC,i_D_g_a_CC->g_a,i_D_lga_CC->sum_by_cat[h],
			    i_D_lga_CC->sqsum_by_cat[h],
			    i_D_lga_CC->glm->suff[h],&beta0,&beta1,&ssq);
	  i_D_lga_CC->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga_CC->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga_CC->glm->ssq[h] = ssq;
	}     
      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = i_D_lga->glm->suff[h][0][0] * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff[h][0][1] * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][2] = i_D_lga->glm->suff[h][0][0] * i_D_lga->haulweight[h]*i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff[h][0][2];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff[h][1][2];
	  i_D_lga->glm->beta_hat[h][0][2] = G_ZERO; /* needed in find b-vector in sample_gauss_eff() */
	}
      if(i_D_lga_CC->glm->nxcov==3)
	{
	  i_D_lga_CC->glm->suff[h][0][2] = i_D_lga_CC->glm->suff[h][0][0] * i_D_lga_CC->haulweight[h];
	  i_D_lga_CC->glm->suff[h][1][2] = i_D_lga_CC->glm->suff[h][0][1] * i_D_lga_CC->haulweight[h];
	  i_D_lga_CC->glm->suff[h][2][2] = i_D_lga_CC->glm->suff[h][0][0] * i_D_lga_CC->haulweight[h]*i_D_lga_CC->haulweight[h];
	  i_D_lga_CC->glm->suff[h][2][0] = i_D_lga_CC->glm->suff[h][0][2];
	  i_D_lga_CC->glm->suff[h][2][1] = i_D_lga_CC->glm->suff[h][1][2];
	  i_D_lga_CC->glm->beta_hat[h][0][2] = G_ZERO; /* needed in find b-vector in sample_gauss_eff() */
	}

      // Copying summaries to the fix-versions. 
      // This should probably be changed to put all directly into the fix-versions
      for(a=0;a<ncat;a++)
	{
	  i_D_lga->Ages_fix[h][a] = i_D_lga->Ages[h][a];
	  i_D_lga->sum_by_cat_fix[h][a] = i_D_lga->sum_by_cat[h][a];
	  i_D_lga->sqsum_by_cat_fix[h][a] = i_D_lga->sqsum_by_cat[h][a];
	} 
      for(a=0;a<ncat_CC;a++)
	{
	  i_D_lga_CC->Ages_fix[h][a] = i_D_lga_CC->Ages[h][a];
	  i_D_lga_CC->sum_by_cat_fix[h][a] = i_D_lga_CC->sum_by_cat[h][a];
	  i_D_lga_CC->sqsum_by_cat_fix[h][a] = i_D_lga_CC->sqsum_by_cat[h][a];
	} 
    }

  FILE *fp;
  fp=fopen("suff_lga_test_init.txt","w");
  //fprintf(fp,"h beta0,beta1,ssq,suff00,suff01,suff10,suff11\n");
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      fprintf(fp,"%d,%f,%f,%f,%f,%f,%f,%f\n",h,
	      i_D_lga->glm->beta_hat[h][0][0],i_D_lga->glm->beta_hat[h][0][1],
	      i_D_lga->glm->ssq[h],i_D_lga->glm->suff[h][0][0],i_D_lga->glm->suff[h][0][1],
	      i_D_lga->glm->suff[h][1][0],i_D_lga->glm->suff[h][1][1]);
    }
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      fprintf(fp,"%d,%f,%f,%f,%f,%f,%f,%f\n",h,
	      i_D_lga_CC->glm->beta_hat[h][0][0],i_D_lga_CC->glm->beta_hat[h][0][1],
	      i_D_lga_CC->glm->ssq[h],i_D_lga_CC->glm->suff[h][0][0],i_D_lga_CC->glm->suff[h][0][1],
	      i_D_lga_CC->glm->suff[h][1][0],i_D_lga_CC->glm->suff[h][1][1]);
    }
  fclose(fp);

  return(0);
}		/* end of makedata_lga_suff_CC */



/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in makedata_lga_suff_CC
*/
int re_makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC)
{
  Fmatrix_2d(&i_D_lga->Ages[0][0],&i_D_lga->Ages[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat[0][0],&i_D_lga->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat[0][0],&i_D_lga->sqsum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->Ages_fix[0][0],&i_D_lga->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat_fix[0][0],&i_D_lga->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat_fix[0][0],&i_D_lga->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat[0][0][0],&i_D_lga->glm->beta_hat[0][0],&i_D_lga->glm->beta_hat[0]);
  FREE(i_D_lga->glm->ssq);
  Fmatrix_3d(&i_D_lga->glm->suff[0][0][0],&i_D_lga->glm->suff[0][0],&i_D_lga->glm->suff[0]);
  
  Fmatrix_2d(&i_D_lga_CC->Ages[0][0],&i_D_lga_CC->Ages[0]);
  Fmatrix_2d(&i_D_lga_CC->sum_by_cat[0][0],&i_D_lga_CC->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga_CC->sqsum_by_cat[0][0],&i_D_lga_CC->sqsum_by_cat[0]);
  Fmatrix_2d(&i_D_lga_CC->Ages_fix[0][0],&i_D_lga_CC->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga_CC->sum_by_cat_fix[0][0],&i_D_lga_CC->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga_CC->sqsum_by_cat_fix[0][0],&i_D_lga_CC->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga_CC->glm->beta_hat[0][0][0],&i_D_lga_CC->glm->beta_hat[0][0],&i_D_lga_CC->glm->beta_hat[0]);
  FREE(i_D_lga_CC->glm->ssq);
  Fmatrix_3d(&i_D_lga_CC->glm->suff[0][0][0],&i_D_lga_CC->glm->suff[0][0],&i_D_lga_CC->glm->suff[0]);
  

  return(0);
}		/* end of re_makedata_lga_suff_CC */


/*!
  \author Geir Storvik
  \brief Calculates sufficient statistics for the wgl model 

  Memory allocated in this routine is reallocated in re_makedata_wgl_suff.
*/
int makedata_wgl_suff(Data_lin *i_D_lin,
		      int *i_nFishBoat,double *i_totlength,double *i_totweight,
		      int *i_replength,double *i_haulweight)
{
  int        f,h,err,n,N,ind;
  double    *x,*y,beta0,beta1,ssq,l,w,r;


  /* Sufficient Statistics */
  i_D_lin->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    N = max(N,i_nFishBoat[h]);

  x = CALLOC(N,double);               // Free ok
  y = CALLOC(N,double);               // Free ok

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      i_D_lin->glm->suff[h][0][0] = G_ZERO;
      i_D_lin->glm->suff[h][0][1] = G_ZERO;
      i_D_lin->glm->suff[h][1][1] = G_ZERO;
      n = 0;
      for(f=0;f<i_nFishBoat[h];f++)
	{
	  l = i_totlength[ind];
	  w = i_totweight[ind];
	  if(l > -9000 && w > -9000)
	    {
	      r = (double) i_replength[ind];
	      x[n] = l;
	      y[n] = w;
	      i_D_lin->glm->suff[h][0][0] += r;
	      i_D_lin->glm->suff[h][0][1] += r*x[n];
	      i_D_lin->glm->suff[h][1][1] += r*x[n]*x[n];
	      n += 1;
	    }
	  ind++;
	}
      if(n > 0)
	{
	  err = lm_fit(n,x,y,&beta0,&beta1,&ssq); 
	  //Note: does not use replength! Ok for starting values?

          i_D_lin->glm->beta_hat[h][0][0] = beta0;
          i_D_lin->glm->beta_hat[h][0][1] = beta1;
          i_D_lin->glm->ssq[h] = ssq;
	}
      i_D_lin->glm->suff[h][1][0] = 
	i_D_lin->glm->suff[h][0][1];
      if(i_D_lin->glm->nxcov==3)
	{
	  i_D_lin->glm->suff[h][0][2] = i_D_lin->glm->suff[h][0][0] * i_haulweight[h];
	  i_D_lin->glm->suff[h][1][2] = i_D_lin->glm->suff[h][0][1]*i_haulweight[h];
	  i_D_lin->glm->suff[h][2][2] = i_D_lin->glm->suff[h][0][0] * i_haulweight[h]*i_haulweight[h];
	  i_D_lin->glm->suff[h][2][0] = i_D_lin->glm->suff[h][0][2];
	  i_D_lin->glm->suff[h][2][1] = i_D_lin->glm->suff[h][1][2];
          i_D_lin->glm->beta_hat[h][0][2] = G_ZERO;    /* needed in find b-vector in sample_gauss_eff() */
	}
    }
  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_wgl_suff */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_wgl_suff
*/
int re_makedata_wgl_suff(Data_lin *i_D_lin,
			 int *i_nFishBoat,double *i_totlength,
			 double *i_totweight,int *i_replength)
{
  Fmatrix_3d(&i_D_lin->glm->beta_hat[0][0][0],&i_D_lin->glm->beta_hat[0][0],&i_D_lin->glm->beta_hat[0]);
  FREE(i_D_lin->glm->ssq);
  Fmatrix_3d(&i_D_lin->glm->suff[0][0][0],&i_D_lin->glm->suff[0][0],&i_D_lin->glm->suff[0]);

  return(0);
}		/* end of re_makedata_wgl_suff */


/*!
  \author Hanne Rognebakke
  \brief Calculates sufficient statistics for the wgl model when coastal cod

  Memory allocated in this routine is reallocated in re_makedata_wgl_suff_CC.
*/
int makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC,
			 int *i_nFishBoat,double *i_totlength,double *i_totweight,
			 int *i_replength,double *i_haulweight,int *i_tottype)
{
  int        f,h,err,N,ind,n_CC,n;
  double     beta0,beta1,ssq,l,w,r;
  double    *x_CC,*y_CC,*x,*y;



  /* Sufficient Statistics */
  i_D_lin->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  i_D_lin_CC->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin_CC->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin_CC->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    N = max(N,i_nFishBoat[h]);

  x_CC = CALLOC(N,double);               // Free ok
  y_CC = CALLOC(N,double);               // Free ok
  x = CALLOC(N,double);               // Free ok
  y = CALLOC(N,double);               // Free ok

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      i_D_lin->glm->suff[h][0][0] = G_ZERO;
      i_D_lin->glm->suff[h][0][1] = G_ZERO;
      i_D_lin->glm->suff[h][1][1] = G_ZERO;
      n_CC = 0;
      n = 0;
      for(f=0;f<i_nFishBoat[h];f++)
	{
	  l = i_totlength[ind];
	  w = i_totweight[ind];
	  if(l > -9000 && w > -9000)
	    {
	      if(i_tottype[ind]==1) //certain coastal cod
		{
		  r = (double) i_replength[ind];
		  x_CC[n_CC] = l;
		  y_CC[n_CC] = w;
		  i_D_lin_CC->glm->suff[h][0][0] += r;
		  i_D_lin_CC->glm->suff[h][0][1] += r*l;
		  i_D_lin_CC->glm->suff[h][1][1] += r*l*l;
		  n_CC += 1;
		}
	      else if(i_tottype[ind]==5) //certain skrei
		{
		  r = (double) i_replength[ind];
		  x[n] = l;
		  y[n] = w;
		  i_D_lin->glm->suff[h][0][0] += r;
		  i_D_lin->glm->suff[h][0][1] += r*l;
		  i_D_lin->glm->suff[h][1][1] += r*l*l;
		  n += 1;
		}
	    }
	  ind++;
	}

      if(n_CC > 0)
	{
	  err = lm_fit(n_CC,x_CC,y_CC,&beta0,&beta1,&ssq); 
	  //Note: does not use replength! Ok for starting values?
          i_D_lin_CC->glm->beta_hat[h][0][0] = beta0;
          i_D_lin_CC->glm->beta_hat[h][0][1] = beta1;
          i_D_lin_CC->glm->ssq[h] = ssq;
	}
      i_D_lin_CC->glm->suff[h][1][0] = i_D_lin_CC->glm->suff[h][0][1];

      if(n > 0)
	{
	  err = lm_fit(n,x,y,&beta0,&beta1,&ssq); 
	  //Note: does not use replength! Ok for starting values?
          i_D_lin->glm->beta_hat[h][0][0] = beta0;
          i_D_lin->glm->beta_hat[h][0][1] = beta1;
          i_D_lin->glm->ssq[h] = ssq;
	}
      i_D_lin->glm->suff[h][1][0] = i_D_lin->glm->suff[h][0][1];


      if(i_D_lin_CC->glm->nxcov==3)
	{
	  i_D_lin_CC->glm->suff[h][0][2] = i_D_lin_CC->glm->suff[h][0][0] * i_haulweight[h];
	  i_D_lin_CC->glm->suff[h][1][2] = i_D_lin_CC->glm->suff[h][0][1]*i_haulweight[h];
	  i_D_lin_CC->glm->suff[h][2][2] = i_D_lin_CC->glm->suff[h][0][0] * i_haulweight[h]*i_haulweight[h];
	  i_D_lin_CC->glm->suff[h][2][0] = i_D_lin_CC->glm->suff[h][0][2];
	  i_D_lin_CC->glm->suff[h][2][1] = i_D_lin_CC->glm->suff[h][1][2];
          i_D_lin_CC->glm->beta_hat[h][0][2] = G_ZERO;    /* needed in find b-vector in sample_gauss_eff() */
	}
      if(i_D_lin->glm->nxcov==3)
	{
	  i_D_lin->glm->suff[h][0][2] = i_D_lin->glm->suff[h][0][0] * i_haulweight[h];
	  i_D_lin->glm->suff[h][1][2] = i_D_lin->glm->suff[h][0][1]*i_haulweight[h];
	  i_D_lin->glm->suff[h][2][2] = i_D_lin->glm->suff[h][0][0] * i_haulweight[h]*i_haulweight[h];
	  i_D_lin->glm->suff[h][2][0] = i_D_lin->glm->suff[h][0][2];
	  i_D_lin->glm->suff[h][2][1] = i_D_lin->glm->suff[h][1][2];
          i_D_lin->glm->beta_hat[h][0][2] = G_ZERO;    /* needed in find b-vector in sample_gauss_eff() */
	}
    }
  FREE(x_CC);
  FREE(y_CC);
  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_wgl_suff_CC */


/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in makedata_wgl_suff_CC
*/
int re_makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC)
{
  Fmatrix_3d(&i_D_lin->glm->beta_hat[0][0][0],&i_D_lin->glm->beta_hat[0][0],&i_D_lin->glm->beta_hat[0]);
  FREE(i_D_lin->glm->ssq);
  Fmatrix_3d(&i_D_lin->glm->suff[0][0][0],&i_D_lin->glm->suff[0][0],&i_D_lin->glm->suff[0]);

  Fmatrix_3d(&i_D_lin_CC->glm->beta_hat[0][0][0],&i_D_lin_CC->glm->beta_hat[0][0],&i_D_lin_CC->glm->beta_hat[0]);
  FREE(i_D_lin_CC->glm->ssq);
  Fmatrix_3d(&i_D_lin_CC->glm->suff[0][0][0],&i_D_lin_CC->glm->suff[0][0],&i_D_lin_CC->glm->suff[0]);

  return(0);
}		/* end of re_makedata_wgl_suff_CC */


/*!
  \author Geir Storvik
  \brief Make a struct of type Data_cov (see caa.h) describing the covariate structure. 

  This structure will typically be a part of a Data_glm struct.

  Memory allocated by this routine is reallocated by re_make_c_cov
*/
static int make_c_cov(int i_n_cov,int i_ispat,int *i_fix,int *i_x_cov,
                      int i_nHaul,Data_cov **o_xcov)
{
  int       f,h,ind;
  Data_cov *xcov;

  xcov = CALLOC(1,Data_cov);    // Free ok
  /*number of factors */
  xcov->n_cov = i_n_cov;
  xcov->n_fac = CALLOC(xcov->n_cov,int);  // Free ok
  xcov->fix = CALLOC(xcov->n_cov,int);    // Free ok
  /* index for spatial term */
  xcov->ispat = i_ispat-1;   
  /* fixed effect */
  for(f=0;f<i_n_cov;f++)
    xcov->fix[f] = i_fix[f];
  xcov->c_cov = CALLOC(i_nHaul,int *);            // Free ok
  for(h=0;h<i_nHaul;h++)
    xcov->c_cov[h] = CALLOC(xcov->n_cov,int);   // Free ok
  ind = 0;
  for(h=0;h<i_nHaul;h++)
    for(f=0;f<i_n_cov;f++)
      {
	xcov->c_cov[h][f] = i_x_cov[ind];
	ind++;
      }
      
  *o_xcov = xcov;
  return(0);
}		/* end of make_c_cov */


 
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in make_c_cov
*/
static int re_make_c_cov(int i_n_cov,int i_ispat,int *i_fix,int *i_x_cov,
                      int i_nHaul,Data_cov **o_xcov)
{
  int       h;
  Data_cov *xcov;

  xcov = *o_xcov;

  FREE(xcov->n_fac);
  FREE(xcov->fix);

  for(h=0;h<i_nHaul;h++)
    FREE(xcov->c_cov[h]);
  FREE(xcov->c_cov);
      
  FREE(xcov);

  return(0);
}		/* end of re_make_c_cov */

    

/*!
  \author Geir Storvik
  \brief Make spatial structure 

  Memory allocated in this routine is reallocated by re_make_spat_struct
*/
static int make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov)
{
  char  string[150];
  int   ind,k,r,rr,nArea,isp;

  isp = x_xcov->ispat;
  nArea = x_xcov->n_fac[isp];
  x_xcov->num_adj_area = CALLOC(nArea,int);          // Free ok
  x_xcov->adj_area = CALLOC(nArea,int *);            // Free ok

  for(r=0;r<nArea;r++)
    {
      x_xcov->num_adj_area[r] = i_num[r];
      x_xcov->adj_area[r] = CALLOC(x_xcov->num_adj_area[r],int);  // Free ok
    }

  ind = 0;
  for(r=0;r<nArea;r++)
    {
      for(rr=0;rr<x_xcov->num_adj_area[r];rr++)
	{
       
	  x_xcov->adj_area[r][rr] = i_adj_area[ind]-1;
          ind++;
	  /* convert neighbor areas  */
          k = 0;
          while(k < (nArea-1)  && x_xcov->adj_area[r][rr]!=x_xcov->conv[isp][k])
	    k++;
          if(x_xcov->adj_area[r][rr]==x_xcov->conv[isp][k])
	    x_xcov->adj_area[r][rr] = k;
          else
	    {
              printf("make_spat_struct:Missing neighbor area %d %d %d\n",
		     r,rr,(int) i_adj_area[ind]);
              sprintf(string,"make_spat_struct:Missing neighbor area %d %d %d\n",
		      r,rr,(int) i_adj_area[ind]);
              write_warning(string);
	      return(1);
	    }
	}
    }

  return(0);
}		/* end of make_spat_struct */

    
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in make_spat_struct
*/
static int re_make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov)
{
  int   r,nArea,isp;

  isp = x_xcov->ispat;
  nArea = x_xcov->n_fac[isp];

  FREE(x_xcov->num_adj_area);

  for(r=0;r<nArea;r++)
    FREE(x_xcov->adj_area[r]);
  FREE(x_xcov->adj_area);


  return(0);
}		/* end of re_make_spat_struct */



/*!
  \author Geir Storvik
  \brief make a struct of type Data_l containing length only and age-stratified-by-length data.

  Memory allocated by this routine is reallocated by re_makedata_only_length
*/
int makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
                         Data_lin *i_D_lga,Data_l **o_D_l)
{
  int         a,i,ind,l;
  Data_l     *D_l;

  D_l = CALLOC(1,Data_l);     // Free ok

  D_l->nLengths = (int) i_nLengths;

  if(D_l->nLengths>0)
    {
      D_l->count = CALLOC(D_l->nLengths,int);   // Free ok
      for(l=0;l<D_l->nLengths;l++)
	D_l->count[l] = (int) i_lengthCount[l];

      D_l->length = CALLOC(D_l->nLengths,double); // Free ok
      for(l=0;l<D_l->nLengths;l++)
	D_l->length[l] = i_length[l];

      D_l->journey = CALLOC(D_l->nLengths,int);   // Free ok
      for(l=0;l<D_l->nLengths;l++)
	{
	  D_l->journey[l] = (int) i_journey[l];
	  D_l->journey[l]--;    /* Change for hauls starting at zero */
	}
    }
  /* Aged fish */
  D_l->nAgeLengths = i_lga_nAgeLengths;
  if(D_l->nAgeLengths>0)
    {
      D_l->ageLength = i_lga_ageLength;
      D_l->ageLengthCount = 
	Mmatrix_2d(0,i_lga_nAgeLengths-1,0,i_nAges-1,sizeof(int),1); // Free ok
      ind = 0;
      for(i=0;i<i_lga_nAgeLengths;i++)
	for(a=0;a<i_nAges;a++)
	  {
	    D_l->ageLengthCount[i][a] = i_lga_ageLengthCount[ind];
	    ind++;
	  }
      D_l->ageJourney = CALLOC(i_lga_nAgeLengths,int);   // Free ok
      for(i=0;i<i_lga_nAgeLengths;i++)
	{
	  D_l->ageJourney[i] = i_lga_ageJourney[i];
          D_l->ageJourney[i]--;
	}
    }
  *o_D_l = D_l;

  return(0);
}		/* end of makedata_only_length */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in makedata_only_length
*/
int re_makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,
			    int *i_journey,
			    int i_lga_nAgeLengths,double *i_lga_ageLength,
			    int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
			    Data_lin *i_D_lga,Data_l **o_D_l)
{
  Data_l     *D_l;

  D_l = *o_D_l;

  if(D_l->nLengths>0)
    {
      FREE(D_l->count);
      FREE(D_l->length);
      FREE(D_l->journey);
    }
  /* Aged fish */
  if(D_l->nAgeLengths>0)
    {
      Fmatrix_2d(&D_l->ageLengthCount[0][0],&D_l->ageLengthCount[0]);
      FREE(D_l->ageJourney);
    }

  FREE(D_l);

  return(0);
}		/* end of re_makedata_only_length */


/*!
  \author Geir Storvik
  \brief Make a struct of type Data_totcatch containing  total catch (to be used for prediction).

  Memory allocated in this routine is reallocated by re_makedata_totcatch
*/
int makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,
                      int *i_n_cov,int i_nCell,int i_nFactors,int *i_inc_haul,
		      int *i_fac_age_int,int *i_fac_age_hsz,
                      int *i_fac_lga_int,int *i_fac_lga_slp,int *i_fac_lga_hsz,
                      int *i_fac_wgl_int,int *i_fac_wgl_slp,int *i_fac_wgl_hsz,
                      int *i_factors,double *i_catch,
                      Data_totcatch **o_D_totcatch)
{
  int        err,c,i,ind,ncov;
  Data_totcatch  *D_totcatch;
 
  D_totcatch = CALLOC(1,Data_totcatch);        // Free ok
  
  /* Number of cells and number of factors */
  D_totcatch->nCell = (int) i_nCell;
  D_totcatch->nFactors = (int) i_nFactors;

  ncov = 1;
  if(i_n_cov[1]>0)
    ncov = 2;
  D_totcatch->fac_age = Mmatrix_2d(0,ncov-1,0,D_totcatch->nFactors-1,sizeof(int),1);    // Free ok

  /* Factors corresponding to age int model */
  for(i=0;i<i_n_cov[0];i++)
    D_totcatch->fac_age[0][i] = (int) i_fac_age_int[i]-1; /* Assume start on zero */
  /* Factors corresponding to age hsz model */
  for(i=0;i<i_n_cov[1];i++)
    D_totcatch->fac_age[1][i] = (int) i_fac_age_hsz[i]-1;/* Assume start on zero */

  ncov = 2;
  if(i_n_cov[4]>0)
    ncov = 3;
  D_totcatch->fac_lga = Mmatrix_2d(0,ncov-1,0,D_totcatch->nFactors,sizeof(int),1);  // Free ok

  /* Factors corresponding to lga intercept model */
  for(i=0;i<i_n_cov[2];i++)
    D_totcatch->fac_lga[0][i] = (int) i_fac_lga_int[i]-1;  /* Assume start on zero */

  /* Factors corresponding to lga slope model */
  for(i=0;i<i_n_cov[3];i++)
    D_totcatch->fac_lga[1][i] = (int) i_fac_lga_slp[i]-1; /* Assume start on zero */

  /* Factors corresponding to lga hsz model */
  for(i=0;i<i_n_cov[4];i++)
    D_totcatch->fac_lga[2][i] = (int) i_fac_lga_hsz[i]-1; /* Assume start on zero */

  ncov = 2;
  if(i_n_cov[7]>0)
    ncov = 3;
  D_totcatch->fac_wgl = Mmatrix_2d(0,ncov-1,0,D_totcatch->nFactors,sizeof(int),1);  // Free ok

  /* Factors corresponding to wgl intercept model */
  for(i=0;i<i_n_cov[5];i++) 
    D_totcatch->fac_wgl[0][i] = (int) i_fac_wgl_int[i]-1; /* Assume start on zero */

  /* Factors corresponding to wgl slope model */
  for(i=0;i<i_n_cov[6];i++)
    D_totcatch->fac_wgl[1][i] = i_fac_wgl_slp[i]-1; /* Assume start on zero */

  /* Factors corresponding to wgl hsz model */
  for(i=0;i<i_n_cov[7];i++)
    D_totcatch->fac_wgl[2][i] = i_fac_wgl_hsz[i]-1; /* Assume start on zero */

  /* factors */
  D_totcatch->factors = Mmatrix_2d(0,D_totcatch->nCell-1,         // Free ok
				    0,D_totcatch->nFactors,sizeof(int),1);
  D_totcatch->catch = CALLOC(D_totcatch->nCell,double);           // Free ok
  ind = 0;
  for(c=0;c<D_totcatch->nCell;c++)
    {
      for(i=0;i<D_totcatch->nFactors;i++)
	{
	  D_totcatch->factors[c][i] = (int) i_factors[ind]-1;
          ind++;
	}
      /* New hauls */
      D_totcatch->factors[c][D_totcatch->nFactors] = -1;
    }

  /* total weight */
  for(c=0;c<D_totcatch->nCell;c++)
      D_totcatch->catch[c] = i_catch[c];

  /* Convert factors for age model */
  D_totcatch->age_xcov = CALLOC(i_D_age->glm->nxcov,Data_cov *);   // Free ok
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      D_totcatch->age_xcov[i] = CALLOC(i_D_age->glm->nxcov,Data_cov);   // Free ok
      err = convert_tot_cov(D_totcatch,D_totcatch->fac_age[0],i_D_age->glm->xcov[i],D_totcatch->age_xcov[i],0);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }

  /* Convert factors for lga model */
  D_totcatch->lga_xcov = CALLOC(i_D_lga->glm->nxcov,Data_cov *);    // Free ok
  for(i=0;i<i_D_lga->glm->nxcov;i++)  
    {
      D_totcatch->lga_xcov[i] = CALLOC(i_D_lga->glm->nxcov,Data_cov);    // Free ok
      if(i==0)
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],1);
      else
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],0);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }

  /* Convert factors for wgl model */
  D_totcatch->wgl_xcov = CALLOC(i_D_wgl->glm->nxcov,Data_cov *);     // Free ok
  for(i=0;i<i_D_wgl->glm->nxcov;i++) 
    {
      D_totcatch->wgl_xcov[i] = CALLOC(i_D_wgl->glm->nxcov,Data_cov);     // Free ok
      if(i==0)
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],1);
      else
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],0);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }

  *o_D_totcatch = D_totcatch;

  return(0);
}		/* end of makedata_totcatch */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated by makedata_totcatch
*/
int re_makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,
			 int i_nCell,int i_nFactors,
			 int *i_fac_age_int,int *i_fac_age_hsz,
			 int *i_fac_lga_int,int *i_fac_lga_slp,int *i_fac_lga_hsz,
			 int *i_fac_wgl_int,int *i_fac_wgl_slp,int *i_fac_wgl_hsz,
			 int *i_factors,double *i_catch,
			 Data_totcatch **o_D_totcatch)
{
  int        err,i;
  Data_totcatch  *D_totcatch;
 
  D_totcatch = *o_D_totcatch;

  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_age[0],i_D_age->glm->xcov[i],D_totcatch->age_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->age_xcov[i]);
    }
  FREE(D_totcatch->age_xcov);
  Fmatrix_2d(&D_totcatch->fac_age[0][0],&D_totcatch->fac_age[0]);

  for(i=0;i<i_D_lga->glm->nxcov;i++)  
    {
      if(i==0)
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],1);
      else
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->lga_xcov[i]);
    }
  FREE(D_totcatch->lga_xcov);         
  Fmatrix_2d(&D_totcatch->fac_lga[0][0],&D_totcatch->fac_lga[0]);


  for(i=0;i<i_D_wgl->glm->nxcov;i++) 
    {
      if(i==0)
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],1);
      else
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->wgl_xcov[i]);
    }
  FREE(D_totcatch->wgl_xcov);
  Fmatrix_2d(&D_totcatch->fac_wgl[0][0],&D_totcatch->fac_wgl[0]);


  Fmatrix_2d(&D_totcatch->factors[0][0],&D_totcatch->factors[0]);
  FREE(D_totcatch->catch);

  FREE(D_totcatch);
  

  return(0);
}		/* end of re_makedata_totcatch */

    
/*!
  \author Geir Storvik
  \brief Converts covariates for total catch in the approperiate order.

  i_neff is the number of effects not being converted 
  (in order to not convert haul effects for linear models)
*/
static int convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff)
{
  int   c,j,k;

  /* Copy from i_xcov */
  o_xcov->n_cov = i_xcov->n_cov;
  o_xcov->n_fac = i_xcov->n_fac;
  o_xcov->fix = i_xcov->fix;
  o_xcov->ispat = i_xcov->ispat;
  o_xcov->icell = i_xcov->icell;
  o_xcov->num_adj_area = i_xcov->num_adj_area;
  o_xcov->adj_area = i_xcov->adj_area;

  /* Convert factors */
  o_xcov->c_cov = Mmatrix_2d(0,i_D_totcatch->nCell-1,     // Free ok
			     0,i_xcov->n_cov-1,sizeof(int),1);

  for(j=0;j<(i_xcov->n_cov-i_neff);j++)
    {
      if(j!=i_xcov->icell)
	{
	  for(c=0;c<i_D_totcatch->nCell;c++)
	    {
	      k = 0;
	      while(k < (i_xcov->n_fac[j]-1)  && 
		    i_D_totcatch->factors[c][i_fac[j]]!=i_xcov->conv[j][k])
		k++;
	      if(i_D_totcatch->factors[c][i_fac[j]]==i_xcov->conv[j][k])
		o_xcov->c_cov[c][j] = k;
	      else
		o_xcov->c_cov[c][j] = -1;
	    }
	}
      else
	{
	  // Assume cell effects are numbered 1,2,....
	  for(c=0;c<i_D_totcatch->nCell;c++)
	    o_xcov->c_cov[c][j] = i_D_totcatch->factors[c][i_fac[j]];
	}
    }
  return(0);
}		/* end of convert_tot_cov */

    
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated bo convert_tot_cov
*/
static int re_convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff)
{
  Fmatrix_2d(&o_xcov->c_cov[0][0],&o_xcov->c_cov[0]);

  return(0);
}		/* end of re_convert_tot_cov */


/*!
  \author Geir Storvik
  \brief Convert covariates to be 0,1,...,
*/
static int convert_cov(int i_nHaul,int *i_nFac,Data_cov *i_xcov)
{
  char    string[150];
  int     h,j,n, *c_ord, *c_ind;
  static int (*cmp)();

  cmp = compare;
  c_ord = CALLOC(i_nHaul,int);    // Free ok
  c_ind = CALLOC(i_nHaul,int);    // Free ok
  i_xcov->conv = CALLOC(i_xcov->n_cov,int *);   // Free ok

  for(j=0;j<i_xcov->n_cov;j++)
    {
      n = 0;
      for(h=0;h<i_nHaul;h++)
	{
	  n = MAX(n,i_xcov->c_cov[h][j]);
          i_xcov->c_cov[h][j]--;
	}
      if(j==i_xcov->ispat)
	n = i_nFac[j];
      i_xcov->n_fac[j] = n;
      if(n!= i_nFac[j])
	{
	  write_warning("convert_cov:Something is wrong\n");
	  sprintf(string,"j=%d,n=%d,nFac=%d\n",j,n,(int) i_nFac[j]);
          write_warning(string);
          for(h=0;h<i_nHaul;h++)
	    {
              sprintf(string,"c_cov[%d][%d]=%d\n",h,j,i_xcov->c_cov[h][j]);
              write_warning(string);
	    }
	  return(1);
	}
      i_xcov->conv[j] = CALLOC(n,int);   // Free ok
      for(n=0;n<i_xcov->n_fac[j];n++)
	i_xcov->conv[j][n] = n;
    }
  FREE(c_ord);
  FREE(c_ind);

  return(0);
}		/* end of convert_cov */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated by convert_cov
*/
static int re_convert_cov(int i_nHaul,int *i_nFac,Data_cov *i_xcov)
{
  int j;

  for(j=0;j<i_xcov->n_cov;j++)
    FREE(i_xcov->conv[j]);
  FREE(i_xcov->conv);

  return(0);
}		/* end of re_convert_cov */


/*!
  \author Geir Storvik
  \brief  Compare two integers givin -1 if x<y and 1 otherwise

   Routine used for qsort.
*/
int compare(int *i_x,int *i_y)
{
 int r=0;
 if(*i_x < *i_y)
   r= -1;
 else if(*i_x > *i_y)
   r=1;
 return(r);
}		/* end of compare */


/*!
  \author Geir Storvik
  \brief  Compare two integers givin -1 if x<y and 1 otherwise

   Routine used for qsort.
*/
int compd(double *i_x,double *i_y)
{
 int r=0;
 if(*i_x < *i_y)
   r= -1;
 else if(*i_x > *i_y)
   r=1;
 return(r);
}		/* end of compare */



/*!
  \author Geir Storvik
  \brief Define node number of effects in graph structure

  \param i_xcov Struct containing covariates for the model
  \param i_nxcov The number of main terms in model 
         (currently 1 for age model, 2 for lga and wgl model)
  \param i_c_in_gr A matrix indicating if factor j of main term i
         is to be included in the graph
  \param o_node An array giving the node index in the graph
      - First index is for the main term (intercept or slope)
      - Second index is for the covariates (constant term, year, seas,...)
      - Third index is for the factor-level inside the covariate

  The node structure is as follows:
  - list For each main terms
     - sub Node for constant term
     - sub Node for year (if present in model)
     - sub Node for seas (if present in model)
     - sub Node for gear (if present in model)
     - sub Node for area (if present in model)
     - sub Node for cell (if present in model)
  For the age model there is only one main term. On the other hand the above
  structure is inside each age-category with all effects for the second age group
  comming after all effects for the first age group and so on.

  For the lga and wgl models there are two main terms, intercept and slope.
*/
int find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node)
{
  int       i,j,k,ind,n_cov;
  int    ***node;
  Data_cov *xcov;

  node = CALLOC(i_nxcov,int **);     // Free ok
  ind = 0;
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      n_cov = xcov->n_cov;
      node[i] = CALLOC(n_cov,int *);     // Free ok

      for(j=0;j<xcov->n_cov;j++)
	{
          if(i_c_in_gr[i][j])
	    {
  	      node[i][j] = CALLOC(xcov->n_fac[j],int);   // Free ok
	      for(k=0;k<xcov->n_fac[j];k++)
		{
		  node[i][j][k] = ind;
		  ind++;
		}
	    }
	}
    }

  *o_node = node;

  return(0);
}		/* end of find_node_effect */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated in find_node_effect
*/
int re_find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node)
{
  int       i,j;
  int    ***node;
  Data_cov *xcov;

  node = *o_node;

  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(i_c_in_gr[i][j])
	    {
  	      FREE(node[i][j]);
	    }
	}
      FREE(node[i]);
    }
  FREE(node);


  return(0);
}		/* end of re_find_node_effect */


/*!
  \author Geir Storvik
  \brief Allocates space for Eff_str
*/
int alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par)
{
  int      a,i,j;
  Eff_str *par;

  par = CALLOC(1,Eff_str);         // Free ok
  par->eff = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(double **),1); //Free ok
  par->ssq = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(double *),1); //Free ok
  par->n_ssq = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(int *),1); //Free ok
  for(a=0;a<i_ncat;a++)
  for(i=0;i<i_nxcov;i++)
    {
      par->eff[a][i] = CALLOC(i_xcov[i]->n_cov,double *);  // Free ok
      for(j=0;j<i_xcov[i]->n_cov;j++)
	par->eff[a][i][j] = CALLOC(i_xcov[i]->n_fac[j],double); // Free ok
      par->ssq[a][i] = CALLOC(i_xcov[i]->n_cov,double);  // Free ok
      par->n_ssq[a][i] = CALLOC(i_xcov[i]->n_cov,int);  // Free ok
    }

  par->ar = CALLOC(i_nxcov,double);    // Free ok
  par->prior_mean = CALLOC(i_nxcov,double **); // Free ok
  par->prior_prec = CALLOC(i_nxcov,double **); 
  par->prior_ar = CALLOC(i_nxcov,double *); 
  par->tau = CALLOC(i_nxcov,double *); // Free ok
  for(i=0;i<i_nxcov;i++)
    {
      par->prior_mean[i] = CALLOC(i_xcov[i]->n_cov,double *);  // Free ok
      par->prior_prec[i] = CALLOC(i_xcov[i]->n_cov,double *); 
      par->prior_ar[i] = CALLOC(2,double);  
      for(j=0;j<i_xcov[i]->n_cov;j++)
	{
	  par->prior_mean[i][j] = CALLOC(i_xcov[i]->n_fac[j],double); // Free ok
	  par->prior_prec[i][j] = CALLOC(2,double); 
	}
      par->tau[i] = CALLOC(i_xcov[i]->n_cov,double);  // Free ok
    }
  par->prior_prec_obs = CALLOC(2,double);

  *o_par = par;
  return(0);
}		/* end of alloc_Eff_str */



/*!
  \author Geir Storvik
  \brief Reallocates memory allocated by alloc_Eff_str
*/
int re_alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par)
{
  int      a,i,j;
  Eff_str *par;

  par = *o_par;

  for(a=0;a<i_ncat;a++)
  for(i=0;i<i_nxcov;i++)
    {
      for(j=0;j<i_xcov[i]->n_cov;j++)
	FREE(par->eff[a][i][j]);
      FREE(par->eff[a][i]);
      FREE(par->ssq[a][i]);
      FREE(par->n_ssq[a][i]);
    }
  Fmatrix_2d(&par->eff[0][0],&par->eff[0]);
  Fmatrix_2d(&par->ssq[0][0],&par->ssq[0]);
  Fmatrix_2d(&par->n_ssq[0][0],&par->n_ssq[0]);

  FREE(par->ar);
  for(i=0;i<i_nxcov;i++)
    {
      for(j=0;j<i_xcov[i]->n_cov;j++)
	{
	  FREE(par->prior_mean[i][j]);
	  FREE(par->prior_prec[i][j]);
	}
      FREE(par->prior_mean[i]);
      FREE(par->prior_prec[i]);
      FREE(par->prior_ar[i]);
      FREE(par->tau[i]);
    }
  FREE(par->prior_mean);
  FREE(par->prior_prec);
  FREE(par->prior_ar);
  FREE(par->prior_prec_obs);
  FREE(par->tau);

  FREE(par);

  return(0);
}		/* end of re_alloc_Eff_str */



/*!
  \author Geir Storvik
  \brief Calculated sum of effects for a given haul
*/
double calc_eff(Data_cov *i_xcov,double **i_eff,int i_h)
{
  int    j,k;
  double mu;

  mu = G_ZERO;
  for(j=0;j<i_xcov->n_cov;j++)
    {
      k = i_xcov->c_cov[i_h][j];
      mu += i_eff[j][k];
    }
  return(mu);
}


/*!
  \author Geir Storvik
  \brief Calculated sum of effects for a given haul
*/
double calc_eff_no_haul(Data_cov *i_xcov,double **i_eff,int i_h)
{
  int    j,k;
  double mu;

  mu = G_ZERO;
  for(j=0;j<(i_xcov->n_cov-1);j++)
    {
      k = i_xcov->c_cov[i_h][j];
      mu += i_eff[j][k];
    }
  return(mu);
}


/*!
  \author Geir Storvik
  \brief Calculated log-density of univariate normal
*/
double ldnorm(double x,double mu,double sigma,double logsigma)
{
  double d,res;

  res = (x-mu)/sigma;

  d = -G_HALF*res*res-logsigma;

  return(d);
}		/* end of ldnorm */


/*!
  \author Geir Storvik
  \brief Updates average of simulations for age-parameters
*/
int update_average_age(int i_n,Age_struct *i_age,Data_age *i_D_age,
                       Age_struct *x_age_mean)
{
  int    a,h,err;
  
  err =  update_average_par(i_n,i_D_age->glm->ncat,i_age->par,x_age_mean->par,
                            i_D_age->glm->nxcov,i_D_age->glm->xcov);
  if(err)
    {
      write_warning("update_average_age:Error calling update_average_par\n");
      return(err);
    }

  /* Update haul effects */
  for(h=0;h<i_D_age->glm->nHaul;h++)
  for(a=0;a<i_D_age->glm->ncat;a++)
    update_mean(&(x_age_mean->alpha[h][a]),i_age->alpha[h][a],i_n);

  /* Update haul precision */
  update_mean(&(x_age_mean->par->tau_obs),i_age->par->tau_obs,i_n);


  return(err);
}		/* end of update_average_age */


/*!
  \author Hanne Rognebakke
  \brief Updates average of simulations for parameters in g_a function
*/
int update_average_g_a(int i_n,Data_g_a *i_D_g_a,Data_g_a *x_g_a_mean)
{
  int    a;
  /* Update g_a */
  for(a=0;a<i_D_g_a->ncat;a++)
    update_mean(&(x_g_a_mean->g_a[a]),i_D_g_a->g_a[a],i_n);
  for(a=0;a<i_D_g_a->g_a_npar;a++)
    update_mean(&(x_g_a_mean->g_a_par[a]),i_D_g_a->g_a_par[a],i_n);

  return(0);
}		/* end of update_average_g_a */


/*!
  \author Geir Storvik
  \brief Updates average of simulations for parameters in lga or wgl model
*/
int update_average_lin(int i_n,LW_struct *i_lin,Data_lin *i_D_lin,
                       LW_struct *x_lin_mean)
{
  int    err;
  
  err =  update_average_par(i_n,i_D_lin->glm->ncat,i_lin->par,x_lin_mean->par,
                            i_D_lin->glm->nxcov,i_D_lin->glm->xcov);
  if(err)
    {
      write_warning("update_average_lin:Error calling update_average_par\n");
      return(err);
    }

  /* Update fish precision */
  update_mean(&(x_lin_mean->par->tau_obs),i_lin->par->tau_obs,i_n);

  return(0);
}		/* end of update_average_lin */


/*!
  \author Geir Storvik
  \brief  Updates average for effects, precisions and ar-coefficients
*/
static int update_average_par(int i_n,int i_ncat,Eff_str *i_par,Eff_str *x_par_mean,
                              int i_nxcov,Data_cov **i_xcov)
{
  int    a,i,j,k;
  Data_cov *xcov;
  
  /* Update linear effects */
  for(a=0;a<i_ncat;a++)
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
      for(k=0;k<xcov->n_fac[j];k++)
	update_mean(&(x_par_mean->eff[a][i][j][k]),i_par->eff[a][i][j][k],i_n);
    }

  /* Update precisions */
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	update_mean(&(x_par_mean->tau[i][j]),i_par->tau[i][j],i_n);
    }
  /* Update ar-coef */
  for(i=0;i<i_nxcov;i++)
    update_mean(&(x_par_mean->ar[i]),i_par->ar[i],i_n);

  return(0);
}		/* end of update_average_par */


/*!
  \author Geir Storvik
  \brief Update the mean of a univariate quantity
*/
int update_mean(double *i_mean,double i_x,int i_n)
{
  double n,n1;
  n = (double) i_n;
  n1 = (double) (i_n+1);

  *i_mean = (n*(*i_mean)+i_x)/n1;

  return(0);
}		/* end of update_mean */


/*!
  \author Geir Storvik
  \brief Writes the current MCMC samples of a model to the mcmc-vector

  For a specific model (age, lga or wgl) all parameters involved are written as
  a int string just after the previous iteration. 

  The format is
  - Effects: categories: main terms:covariates:factors
  - AR-coefficients: main terms
  - Precisions: Random effects and then observation precision
  - Likelihood

  For lga and wgl, there is just one category. For age there is only one main term, for
  lga and wgl there are two, intercept and slope. For age, observation precision correspond to
  haul precision.
*/
int write_it(int i_it,Data_glm *i_glm,Eff_str *i_par)
{
  int  a,i,j,k,n,ind,ncov;
  Data_cov *xcov;

  #ifdef WRITE_HAUL_EFF
  FILE   *unit;
  //unit = fopen("haul_lin.dat","w");
  unit = stderr;
  #endif

  ind = i_it*i_par->num_var;
  /* Write linear effects except haul effect */
  n = 0;
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Effects\n");
  #endif
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      for(k=0;k<xcov->n_fac[j];k++)
		{
		  i_par->mcmc[ind+n] = i_par->eff[a][i][j][k];
		  n++;
                  #ifdef WRITE_HAUL_EFF
		  fprintf(unit,"%d %d %d %d %lf\n",i,a,j,k,i_par->eff[a][i][j][k]);
                  #endif
		}
	    }
	}
    }
  
  /* Write ar-coef */
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Ar\n");
  #endif
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  i_par->mcmc[ind+n] = i_par->ar[i];
	  n++;
          #ifdef WRITE_HAUL_EFF
          fprintf(unit,"%d %lf\n",i_par->ar[i]);
          #endif
	}
    }

  /* Write precisions for random effects*/
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Precisions\n");
  #endif
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      //printf("i=%d,fix[0]=%d\n",i,xcov->fix[0]);
      ncov = xcov->n_cov;
      for(j=0;j<ncov;j++)
	{
          if(!xcov->fix[j])
	    {
	     i_par->mcmc[ind+n] = i_par->tau[i][j];
             n++;
             #ifdef WRITE_HAUL_EFF
             fprintf(unit,"%d %d %lf\n",i,j,i_par->tau[i][j]);
             #endif
	    }
	}
    }

  /* Write observation precision */
  i_par->mcmc[ind+n] = i_par->tau_obs;
  n++;
#ifdef WRITE_HAUL_EFF
  fprintf(unit,"Observation precisions\n");
  fprintf(unit,"%lf \n",i_par->tau_obs);
#endif
  //printf("Tau_obs: n=%d\n",n);

  /* Write loglikelihood */
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Loglikelihood\n");
  fprintf(unit,"%lf\n",i_par->loglik);
  #endif
  i_par->mcmc[ind+n] = i_par->loglik;
  n++;

  if(n!=i_par->num_var)
    {
      printf("write_it:n=%d != num_var=%d\n",n,i_par->num_var);
      write_warning("write_it:Something is wrong\n");
      return(1);
    }

  #ifdef WRITE_HAUL_EFF
  //fclose(unit);
  #endif

  return(0);
}		/* end of write_it */


/*!
  \author Geir Storvik
  \brief  Writes the MCMC samples for the g_a relation to the mcmc-vector

  All parameters are stored sequentially
*/
int write_it_g_a(int i_it,Data_g_a *i_g_a)
{
  int  i,ind;
  
  ind = i_it*i_g_a->g_a_npar;

  for(i=0;i<i_g_a->g_a_npar;i++)
    i_g_a->g_a_mcmc[ind+i] = i_g_a->g_a_par[i];

  return(0);
}		/* end of write_it */


/*!
  \author Geir Storvik
  \brief Picks out parameters for iteration it from mcmc-vectors
*/
int read_it(int i_it,Data_glm *i_glm,Eff_str *i_par)
{
  int  a,i,j,k,n,ind,ncov;
  Data_cov *xcov;
  
  ind = i_it*i_par->num_var;
  /* Read linear effects except haul effect */
  n = 0;
  for(i=0;i<i_glm->nxcov;i++)
    {
      //printf("nxcov:i=%d\n",i);
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  //printf("ncat:a=%d\n",a);
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      //printf("n_cov:j=%d\n",j);
	      for(k=0;k<xcov->n_fac[j];k++)
		{
		  i_par->eff[a][i][j][k] = i_par->mcmc[ind+n];
		  //printf("n_fac:k=%d: eff[%d][%d][%d][%d]=%f\n",k,a,i,j,k,i_par->eff[a][i][j][k]);
		  n++;
		}
	    }
	}
    }

  /* Read ar-coef */
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  i_par->ar[i] = i_par->mcmc[ind+n];
	  //printf("ar[%d]=%f\n",i,i_par->ar[i]);
	  n++;
	}
    }

  /* Read precisions for random effects*/
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      ncov = xcov->n_cov;
      for(j=0;j<ncov;j++)
	{
          if(!xcov->fix[j])
	    {
	     i_par->tau[i][j] = i_par->mcmc[ind+n];
	     //printf("tau[%d][%d]=%f\n",i,j,i_par->tau[i][j]);
             n++;
	    }
	}
    }

  /* Read observation precision */
  i_par->tau_obs = i_par->mcmc[ind+n];
  //printf("tau_obs=%f\n",i_par->tau_obs);
  n++;

  /* Read loglikelihood */
  i_par->loglik = i_par->mcmc[ind+n];
  //printf("loglik=%f\n",i_par->loglik);
  n++;
  //printf("read_it:n=%d i_par->num_var=%d\n",n,i_par->num_var);

  if(n!=i_par->num_var)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"read_it:n=%d != num_par=%d\n",n,i_par->num_var);
      #endif
      printf("read_it:n=%d i_par->num_var=%d\n",n,i_par->num_var);
      write_warning("read_it:Something is wrong\n");
      return(1);
    }

  return(0);
}            /* End of read_it */




/*!
  \author Geir Storvik
  \brief Pick out parameters on the non-linear lga relation for iteration it
*/
int read_it_g_a(int i_it,Data_g_a *i_D_g_a)
{
  int  err,i,ind,a,A;
  
  ind = i_it*i_D_g_a->g_a_npar;

  for(i=0;i<i_D_g_a->g_a_npar;i++)
    i_D_g_a->g_a_par[i] = i_D_g_a->g_a_mcmc[ind+i];

  if(i_D_g_a->g_a_model==0)
    {
      A = i_D_g_a->ncat-1;
      for(a=0;a<i_D_g_a->ncat;a++)
	i_D_g_a->g_a[a] = 
	  (log(i_D_g_a->a_vec[a])-log(i_D_g_a->a_vec[0]))/
	  (log(i_D_g_a->a_vec[A])-log(i_D_g_a->a_vec[0]));
    }
  else if(i_D_g_a->g_a_model==1)
    {
      err = calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,i_D_g_a->g_a);
    }
  return(0);
}		/* end of read_it_g_a */


/*!
  \author Geir Storvik
  \brief  Write mean of simulations to mcmc-vector. 

  Same format as write_it except that there is just one value for each parameter in this case.
*/
int write_it_mean(int i_it,Data_glm *i_glm,Eff_str *i_par,Eff_str *i_par_mean)
{
  int  a,i,j,k,n,ind;
  Data_cov *xcov;
  
  ind = i_it*i_par->num_var;
  /* Write fixed linear effects */
  n = 0;
  for(a=0;a<i_glm->ncat;a++)
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
	  for(k=0;k<xcov->n_fac[j];k++)
	    {
		i_par->mcmc[ind+n] = i_par_mean->eff[a][i][j][k];
                n++;
	    }
	}
    }

  /* Write ar-coef */
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  i_par->mcmc[ind+n] = i_par->ar[i];
	  n++;
	}
    }

  /* Write precisions for random effects*/
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	     i_par->mcmc[ind+n] = i_par_mean->tau[i][j];
             n++;
	    }
	}
    }

  /* Write observation precision */
  i_par->mcmc[ind+n] = i_par_mean->tau_obs;
  n++;

  /* Write loglikelihood */
  i_par->mcmc[ind+n] = i_par_mean->loglik;
  n++;

  if(n!=i_par->num_var)
    {
      write_warning("write_it_mean:Something is wrong\n");
      return(1);
    }

  return(0);
}		/* end of write_it_mean */



/*!
  \author Geir Storvik
  \brief Writes the current MCMC samples of catch at age to the mcmc-vector.

  The format is
  - Total catch: age-categories:length-categories
  - Mean lga and mean wgl for each age-category
*/
int write_it_totcatch(int i_it,int i_nCell,int i_ncat,int i_nlen,
			     TC_struct *i_totcatch,
			     double *i_mean_l,double *i_mean_w)
{
  int  a,l,ind,n;

  ind = i_it*i_ncat*i_nlen;  
  n = 0;
  for(a=0;a<i_ncat;a++)
  for(l=0;l<i_nlen;l++)
    {
      i_totcatch->mcmc[ind+n] = i_totcatch->catch_at_age[a][l];
      n++;
    }

  for(a=0;a<i_ncat;a++)
    {
      i_totcatch->mean_l[i_it*i_ncat+a] = i_mean_l[a];
      i_totcatch->mean_w[i_it*i_ncat+a] = i_mean_w[a];
    }
  return(0);
}		/* end of write_it_totcatch */


/*!
  \author Geir Storvik/Håvard Rue
  \brief propose a new value, scale, on the interval [1/f, f].

  Density of proposal is \f$\propto 1+1/x\f$.  This choice makes \f$q(x,x')/q(x',x) = 1\f$,
  in the acceptance ratio, when \f$x' = scale*x\f$.
*/
double scale_proposal(double x, double f, double *la)
{
    double len = f - 1/f;
    if (la) *la = 0.0;
    if (f == 1.0) return x;
    if ((*GMRFLib_uniform)() < len/(len+2*log(f)))
        return (1/f + len*(*GMRFLib_uniform)())*x;
    else
        return pow(f, 2.0*(*GMRFLib_uniform)()-1.0)*x;
}



/*!
  \author Geir Storvik
  \brief Writes a warning message to screen
*/
void write_warning(char *i_text)
{
  #ifdef LOG_FILE
  fprintf(g_caa_log,i_text);
  #endif
  printf(i_text);

  return;
}		/* end of write_warning */

void write_output(char *filename, char *i_text)
{
  FILE *fp;
  fp = fopen(filename,"a");
  fprintf(fp,i_text);
  fclose(fp);

  return;
}		/* end of write_output */

/*!
  \author Geir Storvik
  \brief Generate an observation from the multinomial distribution 
  \param n Number of events that will be classified into one of
           the categories 1...ncat
  \param p Vector of probabilities.  p(i) is the probability that
           an event will be classified into category i.  Thus, p(i)
           must be [0,1]. Only the first ncat-1 p(i) must be defined
           since P(NCAT) is 1.0 minus the sum of the first
  \param ncat Number of categories.  Length of p and ix.
  \param ix Observation from multinomial distribution.  All ix(i)
            will be nonnegative and their sum will be n.

   Method: Algorithm from page 559 of Devroye, Luc, Non-Uniform Random Variate Generation.  
   Springer-Verlag, New York, 1986.

   This routine is a slight modification of the genmul routine of the ranlib library
*/
void my_genmul(long n,double *p,long ncat,long *ix)
{
static double prob,ptot,sum;
static long i,icat,ntot;
    ptot = 0.0F;
    for(i=0; i<ncat; i++) {
        ptot += *(p+i);
    }
/*
     Initialize variables
*/
    ntot = n;
    sum = 1.0F;
    for(i=0; i<ncat; i++) ix[i] = 0;
/*
     Generate the observation
*/
    for(icat=0; icat<(ncat-1); icat++) {
        prob = *(p+icat)/sum;
        *(ix+icat) = ignbin(ntot,prob);
        ntot -= *(ix+icat);
	if(ntot <= 0) return;
        sum -= *(p+icat);
    }
    *(ix+ncat-1) = ntot;
/*
     Finished
*/
    return;
}

