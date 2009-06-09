int sample_multi_initialize(int i_ncat);
int sample_multi_re_initialize();
int age_haul_modes(int i_start_h,int i_stop_h,Age_struct *i_age,Data_age *i_D_age);
int make_graph_age_fix(Data_glm *i_glm,Age_struct *i_age,int **i_in_gr,int i_start_h);
int sample_ages_len_only_init_new(Data_orig *i_D_orig,Age_struct *i_age,
				  Data_age *i_D_age,LW_struct *i_length,
				  Data_lin *i_D_lga,Data_g_a *i_D_g_a,int saveSim);
int sample_ages_len_only_new(Data_orig *i_D_orig,Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
			     int saveSim,int i_it);
int sample_ages_len_only_init_new_CC(Data_orig *i_D_orig,Data_CC *i_D_CC,
				     Age_struct *i_age,Data_age *i_D_age,
				     Data_lin *i_D_lga,Data_g_a *i_D_g_a,
				     Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,int printSim);
int sample_ages_len_only_new_CC(Data_orig *i_D_orig,Data_CC *i_D_CC,
				Age_struct *i_age,Data_age *i_D_age,
				LW_struct *i_length,LW_struct *i_length_CC,
				Data_lin *i_D_lga,Data_g_a *i_D_g_a,
				Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,
				int printSim);
int make_suff_age(int i_ncat,Age_struct *i_age,Data_age *i_D_age,double *i_haulweight,
		  int i_start_h);
int make_suff_lga(Data_lin *i_D_lga,Data_g_a *i_g_a,int i_start_h);
int sample_g_a_vonBert(Age_struct *i_age,LW_struct *i_length,Data_age *i_D_age,
		       Data_lin *i_D_lga,double **i_suff);
int make_graph_age_ran(Age_struct *i_age,Data_age *i_D_age);
int re_make_graph_age_ran(Age_struct *i_age,Data_age *i_D_age);
int age_haul_init2mean(int i_start_h,int i_stop_h,Age_struct *i_age,Data_age *i_D_age);
int age_stratified(int i_nHaul_noAges,Age_struct *i_age,Data_age *i_D_age,
                   Data_lin *i_D_lga,Data_l *i_D_l);
int sample_age_alpha(Age_struct *i_age,Data_age *i_D_age,
		    int i_start_h,int i_acc,int i_it,int *acc_h);
int sample_age_ran(Age_struct *i_age,Data_age *i_D_age,int i_start_h);
int sample_precision_age_haul(int i_start_h,Eff_str *i_par,
			      double **i_alpha,Data_glm *i_glm);
int sample_discard_init_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,
			     Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_orig *i_D_orig,
			     Data_COST *i_D_COST,int printDisc);
int sample_discard_init2_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,
			      Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_orig *i_D_orig,
			      Data_COST *i_D_COST,int printDisc);
int sample_discard_COST(Age_struct *i_age,Data_age *i_D_age,LW_struct *i_length,
			Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_orig *i_D_orig,
			Data_COST *i_D_COST,int printDisc, int i_it);
int sample_cens_par(Data_lin *i_D_lga, Data_orig *i_D_orig, Data_COST *i_D_COST, int i_it);
int sample_ages_init_COST(Data_orig *i_D_orig,Age_struct *i_age,
			  Data_age *i_D_age,LW_struct *i_length,
			  Data_lin *i_D_lga,Data_g_a *i_D_g_a,Data_COST *i_D_COST,int saveSim);
int sample_ages_init2_COST(Data_orig *i_D_orig,Age_struct *i_age,
			   Data_age *i_D_age,LW_struct *i_length,
			   Data_lin *i_D_lga,Data_g_a *i_D_g_a,
			   Data_COST *i_D_COST,int saveSim);
int sample_ages_COST(Data_orig *i_D_orig,Age_struct *i_age,Data_age *i_D_age,
		     LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		     int saveSim, int i_it);
double cens_function(double lstart,double lmid,double lend,double k,double m,double r);

/*
double haul_loglik(double x, void *i_d,Age_struct *i_age);
*/
