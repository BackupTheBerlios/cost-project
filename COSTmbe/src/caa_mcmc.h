int MCMC_it_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h);
int MCMC_it_g_a(Data_age *i_D_age,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		LW_struct *i_length,int i_start_h,int i_it);
int MCMC_it_wgl(Data_lin *i_D_wgl,LW_struct *i_weight,int start_h);
