/*!
  \file caa.h
  \brief Main header file for the caa system
  \author Geir Storvik
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "GMRFLib.h"
#include "utl.h"
#include "ranlib.h"
#include "numrec.h"
#include "caa_substruct.h"

/* Universal constants */

/*
#define  REAL double
*/
#define  G_PI                3.14159265  /* ? */
#define  G_E                 2.718281    /* ? */
#define  G_ZERO              0.00000000
#define  G_ONE_THIRD         0.33333333 
#define  G_ONE_FOURTH        0.25000000 
#define  G_HALF              0.50000000
#define  G_ONE               1.00000000
#define  G_TWO               2.00000000
#define  G_THREE             3.00000000
#define  G_FOUR              4.00000000
#define  TRUE              1
#define  FALSE             0
#define  CALLOC(n,type) (type *)calloc((size_t)(n),sizeof(type))
#define  FREE(ptr) {if(ptr)free((void *)(ptr));ptr=NULL;}

//Priors or initial values for hyper parameters
#define P_TAU_CON            0.0001
#define P_TAU_FIX            0.001
#define P_TAU_RAN            100.0
#define P_TAU_OBS            1.0
#define P_AR                 0.5
#define P_GAMMA_A            0.01 //parameters in Gamma priors for precisions
#define P_GAMMA_B            0.01 //parameters in Gamma priors for precisions


//#define  LOG_FILE          1 //Writes output to log-file
//#define  DEBUG             1 //Different output for each iteration
//#define  DEBUG_TEST        1 //Testing under programming
//#define  DEBUG_INPUT       1 //Writes input from all main routines to files
//#define  DEBUG_GAUSS       0 //For debugging simulation of linear parts of model
//#define  DEBUG_GAUSS_FILE  0 //For debugging simulation of linear parts of model
//#define  DEBUG_GAUSS_CONSTR  0 //For debugging simulation of linear parts of model
//#define  DEBUG_MULTI       0 //For debugging simulation of nonlinear parts of age-model
//#define  DEBUG_MULTI_TEST  0
//#define  DEBUG_MULTI_H     0
//#define  DEBUG_MULTI_INIT  0
//#define  DEBUG_MULTI_FILE  0
//#define  DEBUG_MULTI_FILE2 0
//#define  DEBUG_MULTI_ALPHA 1
//#define DEBUG_LOGFILE      1
//#define DEBUG_PREDICT      1
//#define DEBUG_G_A          0
//#define DEBUG_G_A2         0
//#define DEBUG_G_A_FILE     0
//#define DEBUG_EVALUATE     0
//#define DEBUG_EVALUATE_FILE 1
//#define WRITE_HAUL_EFF     0
//#define  DEBUG_CELL      0

FILE       *g_caa_log;    /*!< file-unit for log-file */
FILE       *g_caa_alpha;  /*!< file-unit for writing alpha's to file, for debugging */
/*
#define MEMCHK(ptr) GMRFLib_ASSERT(ptr,GMRFLib_EMEMORY)
*/
/* Application specific constants: */

/* ran_unif01() */
/*
#define MULTIPLIER   69069
#define SHIFT            1
#define MODULUS      256*256*256*128
#define INVMOD       ( (REAL) 1 / ((REAL) MODULUS )) / ( (REAL) 2 )
*/

/* Macro definitions */
#define  MAX(A,B)  ((A) > (B) ? (A) : (B))
#define  MIN(A,B)  ((A) < (B) ? (A) : (B))
typedef  long int BOOL;

#define SWAP(a,b) tempr=a;a=b;b=tempr
#define  MAX_FILENAME    120


/*!
  \struct Age_struct caa.h
  \brief Define structure for age parameters

  A struct containing the model-description and the parameters in
  the age model. Typically names x_age where x_ is s_ if a global 
  variable, i_ if an input variable and o_ if an output variable.
*/
typedef struct
{
  Eff_str   *par;      /*!< \brief Parameters to be simulated  */
  /*! \brief  Parameters to be simulated  */
  Eff_str   *par_new;         
  /*! \brief alpha[a][h] in age model */
  double   **alpha;
  /*! \brief mu in Poison approximation of multinomial model */
  double    *mu;
  /*! \brief Graph struct for effect */ 
  Graph_str *gr_str_f;
  /*! \brief Graph struct for random effects, not in use now */
  Graph2_str *gr_str_r;
  /*! \brief Transition matrix for errors in age-readings */
  double   **A2A;
  /*! \brief Number of neighbors in A2A */
  int       *A_Nneigh;
  /*! \brief Neighbors in A2A */
  int      **A_neigh;
  /*! \brief =1 if quadradic shape on constant term of haulsize effect */
  int       hsz_quad;
} Age_struct;

/*! 
  \struct LW_struct caa.h
  \brief Define structure for linear models lga and wgl

  A struct containing the model-description and the parameters in
  the  linear models lga or wgl. Typically named x_length and x_weight.
*/
typedef struct
{
  Eff_str   *par;           /*! \brief Parameters to be simulated */
  Graph_str *gr_str;        /*! \brief Graph struct for effects */ 
  int        cens_model;    /*! \brief = 1 if sampling discards, = 0 otherwise */
  double     cens_k;        /*! \brief k-parameter in censoring */
  double     cens_m;        /*! \brief m-parameter in censoring */
  double     cens_r;        /*! \brief r-parameter in censoring */
  double     cens_Nlim;     /*! \brief limit for number of sampled discards */
  int        fixed_model;   /*! \brief = 0 for simulated parameters, = 1 for fixed parameters */
} LW_struct;

/*!
  \struct TC_struct caa.h
  \brief Define structure for catch-at-age

  A struct containing catch-by-age for the current simulation as well as
  vectors containing simulations from all iterations
*/
typedef struct
{
  /*! \brief catch_at_age for current simulation */
  double     **catch_at_age;
  /*! \brief catch at age for all iterations */
  double      *mcmc;          
  /*! \brief Mean in length given age for current simulation */
  double      *mean_l;
  /*! \brief Mean in weight given age for current simulation */
  double      *mean_w;
} TC_struct;

/*! \struct Data_orig
    \brief Contain the original age and length data of amigo type. 

    Needed when a non-linear lga model is used and sufficient statistics are to be 
    recalculated when non-linear model changes.
    This struct could perhaps be merged together with Data_l.
*/
typedef struct
{
  int         *nFishBoat;   /*!< \brief Number of fish pr boat */
  int         *totage;      /*!< \brief ages of fish */     
  int         *start_noAge; /*!< \brief First position of missing age for haul */
  int         *start_Age;   /*!< \brief First position of observation for haul */
  int         *num_noAge;   /*!< \brief Number of noAge fish in haul */
  int          n_int_len;   /*!< \brief Number of intervals for lengths */
  double      *int_len_lim; /*!< \brief lower limits on intervals for lengths */
  double      *int_len;     /*!< \brief length value for intervals */
  double      *totlength;   /*!< \brief lengths of fish */ 
  double      *haulweight;  /*!< \brief haul sizes */
  int         *replength;   /*!< \brief repetitions fish with equal lengths */        
  int         *season;      /*!< \brief seasons or months */
  int          coastal_cod; /*!< \brief =1 if including coastal cod, =0 otherwise */
  int         *tottype;     /*!< \brief if cod, type of cod (coastal or skrei) */
  int         *discard;     /*!< \brief repetitions with discarded */ 
  int         *n_discard;   /*!< \brief number of discards in each haul */
  int         *landed;      /*!< \brief repetitions with landed */
  int         *n_landed;    /*!< \brief number of landed in each haul */
} Data_orig;

/*!
  \struct Data_CC caa.h
  \brief Contain data for coastal cod
*/
typedef struct
{
  double      *ptype1_CC;   /*!< \brief p(type==1|Coastal cod) */
  double      *ptype1_S;    /*!< \brief p(type==1|Skrei) */
  double      *ptype2_CC;
  double      *ptype2_S;
  double      *ptype4_CC;
  double      *ptype4_S;
  double      *ptype5_CC;
  double      *ptype5_S;
} Data_CC;

/*! 
  \struct Data_age caa.h
  \brief Contain data and covariates for age model
*/
typedef struct
{
  Data_glm  *glm;          /*!< \brief covariates */
  int       *n_h;          /*!< \brief Number of fish in haul */
  int      **Ages;         /*!< \brief Number in each age group for each haul */
  int      **Ages_fix;     /*!< \brief Numbers fixed by observed ages */
  int      **Ages_disc;    /*! \brief Number of discards in each age group for each haul */
  int      **Ages_land;    /*! \brief Number of landed in each age group for each haul */
  int       *a_vec;        /*!< \brief Vector defining sequence of ages */
  int       *type_age;        
  /*!< \brief Type of age observations
  = 0 if ages are given without errors,
  = 1 if ages are given with errors,
  = 2 if only lengths are given,
  = 3 if only lengths with some ages are given,
  = 4 if only lengths with some ages with errors are given 
  */
} Data_age;

/*! 
  \struct Data_lin caa.h
  \brief Contain data and covariates for linear model (lga and wgl)
*/
typedef struct
{
  
  Data_glm  *glm; /*!< \brief covariates */
  
  int      **Ages;             /*!< \brief Number in each age group for each haul
			            Ages used in lga model, could be smaller than Ages 
			            in Data_age.  */
  int      **Ages_fix;         /*!< \brief Numbers fixed by observed ages */
  double   **Ages_miss_sim;    /*!< \brief Matrix to store simulated missing ages 
                                    and corresponding observed lengths*/
  double    *haulweight;       /*!< \brief weight in hauls */
  double   **sum_by_cat;       /*!< \brief For lga model, sum of length by age */
  double   **sum_by_cat_fix;   /*!< \brief Fixed by the aged observations */
  double   **sqsum_by_cat;     /*!< \brief For lga model, square sum of length by age */
  double   **sqsum_by_cat_fix; /*!< \brief Fixed by the aged observations */
} Data_lin;

/*!
  \struct Data_g_a caa.h
  \brief Contain parameters for g-function in lga model
*/
typedef struct
{
  /*! \brief Function of a giving linear model between log-length and age */
  double    *g_a;
  /*! \brief  = 0 for log-linear model, = 1 for Schute-Richards model */
  int        g_a_model;
  /*! \brief Number of age intervals */
  int        ncat;
  /*! \brief Number of seasons */
  int        nSeason;
  /*! \brief Vector defining sequence of ages */
  double     *a_vec;           
  /*! \brief Vector defining a_vec corresponding to a_vec in Data_age */
  int       *a2Age_vec;           
  /*! \brief Number of parameters describing g_a function */
  double     g_a_npar;
  /*! \brief Parameters describing g_a function */
  double    *g_a_par;
  /*! \brief vector containing all simulations of g_a parameters */
  double    *g_a_mcmc;
  /*! \brief Sufficient statistics */
  double   **suff;             
} Data_g_a;

/*! 
  \struct Data_l caa.h
  \brief Contain data and covariates for length-only and stratified by length data
*/
typedef struct
{
  /*! \brief Number of lengths in data, for non-Amigo boats */
  int     nLengths;      
  /*! \brief Integers corresponding to different hauls */
  int    *journey;       
  /*! \brief Number of fish in each length category */
  int    *count;         
  /*! \brief Length for each category  */
  double *cat;           
  /*! \brief Length for length category */
  double *length;        
  /*! \brief Number of length-categories for which fish are aged */
  int     nAgeLengths;   
  /*! \brief Lengths for which fish are aged */
  double *ageLength;     
  /*! \brief Number of aged fish in each length/age combination */
  int   **ageLengthCount;
  /*! \brief Corresponding journey */
  int    *ageJourney;    
  /*! \brief Number in each age category corresponding to length */
  int   **Ages;          
} Data_l;                     /* Only length data */


/*! 
  \struct Data_totcatch caa.h
  \brief Contain data and covariates for total catch

  To be used for prediction of catch-at-age
*/
typedef struct
{
  /*! \brief Number of cells with total catch */
  int          nCell;    
  /*! \brief Number of length intervals to divide predictions on */
  int          nlint;    
  /*! \brief Limits of length intervals */
  double      *l_int;    
  /*! \brief Number of factors */
  int          nFactors; 
  /*! \brief Factors corresponding to age

  negative means this factor is not included in age model */
  int        **fac_age;  
  /*! \brief Factors corresponding to lga

  Negative means this factor is not included in age model */
  int        **fac_lga;  
  /*! \brief Factors corresponding to wgl

    Negative means this factor is not included in age model */
  int        **fac_wgl;  
  /*! \brief Factors corresponding to catch */
  int        **factors;  
  /*! \brief Catch  */
  double      *catch;    
  /*! \brief Factors for age model */
  Data_cov   **age_xcov; 
  /*! \brief Factors for lga model */
  Data_cov   **lga_xcov; 
  /*! \brief Factors for wgl model */
  Data_cov   **wgl_xcov; 
} Data_totcatch;    /* Data for total catch */




/*! 
  \struct Data_obs caa.h
  \brief Contain the original data for observer data

  To be used in COST project
*/
typedef struct
{
  /*! \brief number of trips */
  int          n_trip;
  /*! \brief number of hauls pr trip */
  int         *num_trip;
  /*! \brief number of measured discarded fish pr haul */
  int         *num_haul_disc;
  /*! \brief observed month */
  int         *season;
  /*! \brief length categories for discard samples */
  double      *l_disc;
  /*! \brief number at length for discards */
  int         *lfreq_disc;
  /*! \brief number of discards in haul */
  double      *haulsize_disc;
  /*! \brief number of discards sampled */
  double      *sampsize_disc;
  /*! \brief numbers of discard age-length data in trip */
  int         *num_alk;
  /*! \brief ages for discard age-length data */
  int         *alk_a;
  /*! \brief lengths for discard age-length data */
  double      *alk_l;
  /*! \brief number at length for discard age-length data */
  int         *alk_lfreq;
  /*! \brief number of size classes pr trip with landings */
  int         *num_trip_land;
  /*! \brief number of measured landed fish pr size class */
  int         *num_size_land;
  /*! \brief length categories for landing samples */
  double      *l_land;
  /*! \brief number at length for landings */
  int         *lfreq_land;
  /*! \brief total weight landed in size class */
  double      *totsize_land;
  /*! \brief weight of landings sampled for lengths in size class */
  double     *sampsize_land;
} Data_obs; /* Observer data */

/*! 
  \struct cens_struct caa.h
  \brief Contain the parameters in the censoring model

  To be used in COST project
*/
typedef struct
{
  int         ncat;          /*! \brief number of categories for censoring parameters, e.g. n_trip */
  double      k;             /*! \brief k-parameter in censoring (same for all hauls) */
  double      m;             /*! \brief m-parameter in censoring (same for all hauls) */
  double     *r;             /*! \brief r-parameter in censoring */
  double     *mu;            /*! \brief expected value for censoring parameters */
  double     *tau;           /*! \brief precision for censoring parameters */
  double      a_prior;       /*! \brief parameter in gamma prior */ 
  double      b_prior;       /*! \brief parameter in gamma prior */ 
  double      mu_prior_mean; /*! \brief prior mean for expected value */
  double      mu_prior_prec; /*! \brief prior precision for expected value */
} cens_struct;  /* Parameters in censoring model */

/*!
  \struct Data_mland caa.h
  \brief Contain the original data for market landings

  To be used in COST project
*/
typedef struct
{
  /*! \brief number of trips */
  int          n_trip;
  /*! \brief number of size classes pr trip*/
  int         *num_trip;
  /*! \brief observed month */
  int         *season;
  /*! \brief number of measured fish pr size class */
  int         *num_size;
  /*! \brief length categories for landing samples */
  double      *l;
  /*! \brief number at length */
  int         *lfreq;
  /*! \brief total weight landed in size class */
  double      *totsize;
  /*! \brief weight of landings sampled for lengths in size class */
  double     *sampsize;
  /*! \brief numbers of age-length data in trip */
  int         *num_alk;
  /*! \brief ages for age-length data */
  int         *alk_a;
  /*! \brief lengths for age-length data */
  double      *alk_l;
  /*! \brief number at length for age-length data */
  int         *alk_lfreq;
  /*! \brief number of length categories for simulated discards */
  int          N_int_disc;
  /*! \brief length categories for simulated discards */
  double      *l_disc;
  /*! \brief number at length for simulated discards */
  int         *lfreq_disc;
  /*! \brief upper limits on intervals for lengths */
  double      *int_len_lim;
  /*! \brief parameter in Poisson distribution for number of fish */
  double      *lambda;
  /*! \brief hyperparameter for lambda */
  double       c;
  /*! \brief hyperparameter for lambda */
  double       d;
} Data_mland; /* Market landing data */

/*! 
  \struct Data_COST caa.h
  \brief Contain the original data with both observer data and market landing data
  and the simulated parameters in censoring model

  To be used in COST project
*/
typedef struct
{
  /*! \brief Observer data */
  Data_obs    *obs;
  /*! \brief Market landing data */
  Data_mland  *mland;
  /*! \brief parameters in censoring model */
  cens_struct *cens;
  /*! \brief Number of COST specific variables to be simulated */
  int        num_var;
  /*! \brief MCMC samples of COST specific variables from all iterations */
  double      *mcmc; 
} Data_COST;  /* Original observer and market landing data */




#include "caa_main.h"
#include "caa_mcmc.h"
#include "caa_init.h"
#include "caa_routines.h"
#include "caa_cell_constr.h"
#include "caa_sample_gauss.h"
#include "caa_sample_multi.h"
#include "caa_sample_g_a.h"
#include "caa_evaluate.h"
#include "caa_predict.h"
#include "caa_chol.h"
#include "caa_lqp.h"
#include "caa_util.h"
#include "caa_utl.h"
#include "caa_summaries.h"
#include "caa_evaluate.h"
#include "caa_Rextensions.h"
#include "caa_COST.h"

