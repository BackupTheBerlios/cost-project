int init_evaluate(int i_nHaul,int i_ncat);
int re_init_evaluate(int i_ncat);
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
		   int *o_err);
int calc_lik_age(Age_struct *i_age,Data_age *i_D_age, double *o_loglik);
int Bayes_CV_model1(int i_it,Age_struct *i_age,Data_age *i_D_age,
		    LW_struct *i_lga,Data_lin *i_D_lga,
		    double *o_mean_inv_lik_mod1);
int Bayes_CV_age(int i_it,Age_struct *i_age,Data_age *i_D_age,
                 double *o_mean_inv_loglik_age);
double calc_lik_age_h(int i_h,Age_struct *i_age,Data_age *i_D_age);
int calc_lik_lin(LW_struct *i_lin,Data_lin *i_D_lin,double *o_loglik);
int Bayes_CV_lin(int i_it,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *o_mean_inv_loglik_lin);
double calc_lik_lin_h(int i_h,LW_struct *i_lin,Data_lin *i_D_lin,
                 double *w_res);
int calc_resid_lga(int *i_totage, double *i_totlength,int *i_nFishBoat,
		   int *i_replength,
                   Age_struct *i_age, Data_age *i_D_age,LW_struct *i_lin,Data_lin *i_D_lin,
		   Data_g_a *i_D_g_a,double *o_resid);
int calc_resid_wgl(double *i_totlength, double *i_totweight,int *i_nFishBoat,
                   LW_struct *i_weight,Data_lin *i_D_wgl,
		   double *o_resid);
int calc_KS_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_p);
int calc_D_Robins_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_p);
int calc_entropy_age(Age_struct *i_age,Data_age *i_D_age,double *o_d,double *o_d2,
		      double *o_p);

