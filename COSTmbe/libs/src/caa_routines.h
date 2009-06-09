int  write_input_model1_new(SEXP i_mcmc_par,SEXP i_constr,SEXP i_seed,SEXP i_num_par,SEXP i_nBoats, 
			    SEXP i_common_par,SEXP i_dataList,SEXP i_ageList,SEXP i_lgaList,
			    SEXP i_priorList);
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
			int *i_hsz_quad,
			int *i_lga_cens_model,double *i_lga_cens,
			int *i_num_adj_area,int *i_adj_area,
			int *i_num_par,
			double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
			double *i_pri_age_prec_par,double *i_pri_age_ar,
			double *i_pri_lga_eff_mean,double *i_pri_lga_eff_prec,
			double *i_pri_lga_prec_par,double *i_pri_lga_ar);
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
			int *i_inc_haul,
			int *i_tot_factors,double *i_tot_catch,double *i_par_haulsize,
			int *i_N_l_int,double *i_l_int,int *i_nMC);
int write_input_predict_new(SEXP i_mcmc_samp,SEXP i_common_par,
			    SEXP i_data_age,SEXP i_data_lga,SEXP i_data_wgl,
			    SEXP i_data_catch,
			    SEXP i_par_haulsize,SEXP i_dist_cell,
			    SEXP i_N_l_int,SEXP i_l_int,SEXP i_nMC);
int write_input_model2(int *i_mcmc_par,
		       int *i_constr,int *i_seed,
		       int *i_wgl_nBoats,double *i_totlength,double *i_totweight,double *i_haulweight,
		       int *i_replength,int *i_nFishBoat,int *i_n_cov,int *i_ispat,
		       int *i_wgl_int_nFac,int *i_wgl_int_fix,int *i_wgl_int_c_cov,
		       int *i_wgl_slp_nFac,int *i_wgl_slp_fix,int *i_wgl_slp_c_cov,
		       int *i_wgl_hsz_nFac,int *i_wgl_hsz_fix,int *i_wgl_hsz_c_cov,
		       int *i_num_adj_area,int *i_adj_area,
		       int *i_num_par);
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
			  int *i_num_adj_area,int *i_adj_area);
int makedata_orig(SEXP i_dataList, Data_orig **o_D_orig);
int re_makedata_orig(Data_orig **o_D_orig);
int makedata_CC(SEXP i_ageList, Data_CC **o_D_CC);
int re_makedata_CC(Data_CC **o_D_CC);
int makedata_age1(int i_nBoats,int i_nAges,int *i_a_vec,
		  int i_int_n_cov,int *i_int_nFac,int i_int_ispat,
		  int *i_int_fix,int *i_int_x_cov,
		  int *i_int_num,int *i_int_adj_area,
		  int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,
		  int *i_hsz_fix,int *i_hsz_x_cov,
		  int *i_hsz_num,int *i_hsz_adj_area,
		  Data_age **o_D_age);
int re_makedata_age1(int i_nBoats,int i_nAges,int *i_a_vec,
		     int i_int_n_cov,int *i_int_nFac,int i_int_ispat,
		     int *i_int_fix,int *i_int_x_cov,
		     int *i_int_num,int *i_int_adj_area,
		     int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,
		     int *i_hsz_fix,int *i_hsz_x_cov,
		     int *i_hsz_num,int *i_hsz_adj_area,
		     Data_age **o_D_age);
int makedata_age2(Data_orig *i_D_orig,Data_age *i_D_age);
int re_makedata_age2(Data_age *i_D_age);
int makedata_g_a(int lga_g_a_ncat,int lga_g_a_nSeason,double *i_avec,int *i_a2Age_vec,
		 int i_g_a_model,Data_age *i_D_age,Data_g_a **o_g_a);
int re_makedata_g_a(Data_g_a **o_g_a);
int makedata_lin1(int i_nBoats,
		  int i_int_n_cov,int *i_int_nFac,int i_int_ispat,int *i_int_fix,int *i_int_x_cov,
		  int *i_int_num,int *i_int_adj_area,
		  int i_slp_n_cov,int *i_slp_nFac,int i_slp_ispat,int *i_slp_fix,int *i_slp_x_cov,
		  int *i_slp_num,int *i_slp_adj_area,
		  int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,int *i_hsz_fix,int *i_hsz_x_cov,
		  int *i_hsz_num,int *i_hsz_adj_area,
		  Data_lin **o_D_lin);
int re_makedata_lin1(int i_nBoats,
		     int i_int_n_cov,int *i_int_nFac,int i_int_ispat,int *i_int_fix,int *i_int_x_cov,
		     int *i_int_num,int *i_int_adj_area,
		     int i_slp_n_cov,int *i_slp_nFac,int i_slp_ispat,int *i_slp_fix,int *i_slp_x_cov,
		     int *i_slp_num,int *i_slp_adj_area,
		     int i_hsz_n_cov,int *i_hsz_nFac,int i_hsz_ispat,int *i_hsz_fix,int *i_hsz_x_cov,
		     int *i_hsz_num,int *i_hsz_adj_area,
		     Data_lin **o_D_lin);
int makedata_lga_suff(Data_lin *i_D_lga,Data_orig *i_D_orig,Data_g_a *i_D_g_a);
int re_makedata_lga_suff(Data_lin *i_D_lga);
int makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC,Data_orig *i_D_orig,
			 Data_g_a *i_D_g_a,Data_g_a *i_D_g_a_CC);
int re_makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC);
int makedata_wgl_suff(Data_lin *i_D_lin,
		      int *i_nFishBoat,double *i_totlength,double *i_totweight,
		      int *i_replength,double *i_haulweight);
int re_makedata_wgl_suff(Data_lin *i_D_lin,
			 int *i_nFishBoat,double *i_totlength,double *i_totweight,
			 int *i_replength);
int makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC,
			 int *i_nFishBoat,double *i_totlength,double *i_totweight,
			 int *i_replength,double *i_haulweight,int *i_tottype);
int re_makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC);
int makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
                         Data_lin *i_D_lga,Data_l **o_D_l);
int re_makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
			    Data_lin *i_D_lga,Data_l **o_D_l);
int makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,
                      int *i_n_cov,int i_nCell,int i_nFactors,int *i_inc_haul,
		      int *i_fac_age_int,int *i_fac_age_hsz,
                      int *i_fac_lga_int,int *i_fac_lga_slp,int *i_fac_lga_hsz,
                      int *i_fac_wgl_int,int *i_fac_wgl_slp,int *i_fac_wgl_hsz,
                      int *i_factors,double *i_catch,
                      Data_totcatch **o_D_totcatch);
int re_makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,
			 int i_nCell,int i_nFactors,
			 int *i_fac_age_int,int *i_fac_age_hsz,
			 int *i_fac_lga_int,int *i_fac_lga_slp,int *i_fac_lga_hsz,
			 int *i_fac_wgl_int,int *i_fac_wgl_slp,int *i_fac_wgl_hsz,
			 int *i_factors,double *i_catch,
			 Data_totcatch **o_D_totcatch);
int compd(double *i_x,double *i_y);
int find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node);
int re_find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node);
int alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par);
int re_alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par);
double calc_eff(Data_cov *i_xcov,double **i_eff,int i_h);
double calc_eff_no_haul(Data_cov *i_xcov,double **i_eff,int i_h);
double calc_eff_fix(Data_cov *i_xcov,double **i_eff,int i_h);
double calc_eff_ran(Data_cov *i_xcov,double **i_eff,int i_h);
double ldnorm(double x,double mu,double sigma,double logsigma);
int write_it(int i_it,Data_glm *i_glm,Eff_str *i_par);
int write_it_g_a(int i_it,Data_g_a *i_g_a);
int read_it(int i_it,Data_glm *i_glm,Eff_str *i_par);
int read_it_g_a(int i_it,Data_g_a *i_D_g_a);
int write_it_mean(int i_it,Data_glm *i_glm,Eff_str *i_par,Eff_str *i_par_mean);
int write_it_totcatch(int i_it,int i_nCell,int i_ncat,int i_nlen,TC_struct *i_totcatch,
		      double *i_mean_l,double *i_mean_w);
int update_average_age(int i_n,Age_struct *i_age,Data_age *i_D_age,
                       Age_struct *x_age_mean);
int update_average_g_a(int i_n,Data_g_a *i_D_g_a,Data_g_a *x_g_a_mean);
int update_average_lin(int i_n,LW_struct *i_lin,Data_lin *i_D_lin,
                       LW_struct *x_lin_mean);
int update_mean(double *i_mean,double i_x,int i_n);
double scale_proposal(double x, double f, double *la);
void write_error(char *i_text);
void write_warning(char *i_text);
void my_genmul(long n,double *p,long ncat,long *ix);
void write_output(char *filename, char *i_text);
