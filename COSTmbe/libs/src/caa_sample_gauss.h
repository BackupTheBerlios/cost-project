
int make_graph_gauss(Data_glm *i_glm,Graph_str *i_gr,Eff_str *i_par,
                     int i_start_h,int i_constr,int i_hsz_quad);
int remake_graph_gauss(Graph_str *i_gr);
int make_Q_gauss(Graph_str *i_gr,Data_glm *i_glm,Eff_str *i_par,
                 int i_start_h,double **o_Q,int *o_no_data);
int make_constr(Graph_str *i_gr,Data_glm *i_glm,int i_hsz_quad);
int re_make_constr(Graph_str *i_gr);
int make_constr2(Graph_str *i_gr,Data_glm *i_glm);
int sample_gauss_eff(Graph_str *i_gr,Eff_str *i_par,Data_glm *i_glm,int i_start_h);
int re_make_graph_gauss(Graph_str *i_gr);
int sample_precision(int i_start_h,Eff_str *i_par,Data_glm *i_glm);
int sample_precision_lin(int i_start_h,Eff_str *i_par,Data_glm *i_glm);
