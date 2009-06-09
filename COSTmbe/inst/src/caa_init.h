int init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_coastal_cod,
	     double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
	     double *i_pri_age_prec_par,double *i_pri_age_ar,
             Age_struct **o_age,Age_struct **o_age_mean);
int re_init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,
		Age_struct **o_age,Age_struct **o_age_mean);
int alloc_age(Data_age *i_D_age,Age_struct **o_age);
int re_alloc_age(Data_age *i_D_age,Age_struct **o_age);
int init_lin(Data_lin *i_D_lin,
	     double *i_pri_lin_eff_mean,double *i_pri_lin_eff_prec,
	     double *i_pri_lin_prec_par,double *i_pri_lin_ar,
	     LW_struct **o_lin,LW_struct **o_lin_mean);
int re_init_lin(Data_lin *i_D_lin,LW_struct **o_lin,LW_struct **o_lin_mean);
int alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin);
int re_alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin);
int init_lga_par(LW_struct *i_length,Data_lin *i_D_lga,
		 double *i_lga_fixed_int,double *i_lga_fixed_slp,double *i_lga_fixed_tau);
int init_glm_sim(Data_glm *i_glm,Data_glm **o_glm_sim);
int re_init_glm_sim(Data_glm *i_glm,Data_glm **o_glm_sim);
