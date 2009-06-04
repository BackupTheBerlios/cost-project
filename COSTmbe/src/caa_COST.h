int makedata_COST(SEXP i_COSTList, Data_orig **o_D_orig, Data_COST **o_D_COST);
int re_makedata_COST(Data_orig **i_D_orig, Data_COST **i_D_COST);
int init_cens_COST(Data_orig *i_D_orig, Data_COST *i_D_COST,
		   double *i_cens_mu, double *i_cens_tau, double *i_cens_pri);
int write_input_model1_COST(Data_orig *i_D_orig, Data_COST *i_D_COST,
			    SEXP i_ageList,SEXP i_lgaList,SEXP i_priorList);
int write_it_COST(int i_it,Data_COST *i_D_COST);
int read_it_COST(int i_it, Data_COST *i_D_COST);
int makedata_COST_predict(int i_nHaul,Data_COST **o_D_COST);
int sample_lambda_init_COST(Data_age *i_D_age,Data_COST *i_D_COST);
int sample_lambda_prior_COST(Data_COST *i_D_COST);
int resample_data_COST(Data_orig *i_D_orig, Data_COST *i_D_COST);
