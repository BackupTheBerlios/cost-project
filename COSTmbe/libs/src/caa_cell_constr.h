int make_cell_constr_age(Data_age *i_D_age,int *i_icell,
			 double *i_age_int_Sigma_cell,double *i_age_int_constr_cell,
			 int i_age_int_nconstr_cell,
			 double *i_age_hsz_Sigma_cell,double *i_age_hsz_constr_cell,
			 int i_age_hsz_nconstr_cell);
int make_cell_constr_lin(Data_lin *i_D_lin,int *i_icell,
			 double *i_lin_int_Sigma_cell,double *i_lin_int_constr_cell,
			 int i_lin_int_nconstr_cell,
			 double *i_lin_slp_Sigma_cell,double *i_lin_slp_constr_cell,
			 int i_lin_slp_nconstr_cell,
			 double *i_lin_hsz_Sigma_cell,double *i_lin_hsz_constr_cell,
			 int i_lin_hsz_nconstr_cell);
int makedata_cell_dist(SEXP i_dist_cell,int *i_icell,
		       Age_struct *i_age,Data_age *i_D_age,
		       LW_struct *i_length,Data_lin *i_D_lga,
		       LW_struct *i_weight,Data_lin *i_D_wgl);
int re_makedata_cell_dist(SEXP i_dist_cell,int *i_icell,
			  Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl);
int simulate_cell_effects(Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl);
int makedata_cell_dist_CC(SEXP i_dist_cell,int *i_icell,
			  Age_struct *i_age,Data_age *i_D_age,
			  LW_struct *i_length,Data_lin *i_D_lga,
			  LW_struct *i_weight,Data_lin *i_D_wgl,
			  LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			  LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC);
int re_makedata_cell_dist_CC(SEXP i_dist_cell,int *i_icell,
			     Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC);
int simulate_cell_effects_CC(Age_struct *i_age,Data_age *i_D_age,
			     LW_struct *i_length,Data_lin *i_D_lga,
			     LW_struct *i_weight,Data_lin *i_D_wgl,
			     LW_struct *i_length_CC,Data_lin *i_D_lga_CC,
			     LW_struct *i_weight_CC,Data_lin *i_D_wgl_CC);
