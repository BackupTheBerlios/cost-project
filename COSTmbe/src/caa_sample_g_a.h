int sample_g_a_initialize(int i_ncat,int i_g_a_model);
int sample_g_a_re_initialize();
int suff_g_a_init(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		  Data_g_a *i_D_g_a,int i_start_h,double **o_suff);
int suff_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
	     Data_g_a *i_D_g_a,int i_start_h,double **o_suff);
int calc_g_a(int i_ncat,double *i_avec,double *i_par,double *o_g);
int sample_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
               Data_g_a *i_D_g_a,double **i_suff,int i_it);
int sample_g_a_S_R(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		   Data_g_a *i_D_g_a,double **i_suff);
int fit_g_a_init(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
		 Data_g_a *i_D_g_a,double **i_suff);
