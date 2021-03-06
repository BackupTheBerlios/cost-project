/* geo.h
 * 
 * Copyright (C) 2001 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The author's contact information:
 *
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 * RCSId: $Id: geo.h,v 1.1 2009/06/09 10:30:47 mmerzere Exp $
 *
 */

/* This code was initially written by Hanne T Wist, later modified by H. Rue */

/*!
  \file geo.h
  \brief Typedefs and defines for \ref geo.c
*/

#ifndef __GMRFLib_GEO_H__
#define __GMRFLib_GEO_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/*!
  \brief Use Matern CF
*/
#define GMRFLib_CORTP_MATERN    0

/*!
  \brief Use Exponential CF
*/
#define GMRFLib_CORTP_EXP       1

/*!
  \brief Use Gaussian CF
*/
#define GMRFLib_CORTP_GAUSS     2

/*!
  \struct GMRFLib_geo_problem_tp geo.h
  \brief Specification of the geo sampling problem.
*/
typedef struct
{
    /*!
      \brief The graph
    */
    GMRFLib_graph_tp *graph;

    /*!
      \brief The Q-function, computing Q(i,j).
    */
    GMRFLib_Qfunc_tp *Qfunc;

    /*!
      \brief The arguments to the Q-function.
    */
    char *Qfunc_arg;
}
GMRFLib_geo_problem_tp;
 

/*!
  \struct GMRFLib_geoQfunc_arg_tp geo.h
*/
typedef struct
{
    int       name;
    int       neigh;
    double    param;
    double    range;
    int       nrow;
    int       ncol;
    double    *prec;
    int       cyclic_flag; /* if != 0, the grid is made cyclic */
    double    *coef;
    map_ii    *hash_table;
}
GMRFLib_geoQfunc_arg_tp;


typedef struct
{
    map_stri coef2;
    map_stri coef3;
}
GMRFLib_Global_geo_tp;


char *GMRFLib_geo_translate_cortp(int name);
char *GMRFLib_geo_translate_neigh(int neigh);
double *GMRFLib_get_geo_coefs(int name, int neigh, double param, double range);
double *GMRFLib_get_geo_coefs2(int name, double param, double range);
double *GMRFLib_get_geo_coefs3(int name, double param, double range);
double GMRFLib_geoQfunc(int node, int nnode, char *arg);
int GMRFLib_free_geo_problem(GMRFLib_geo_problem_tp *geo_problem);
int GMRFLib_init_geo_problem(GMRFLib_geo_problem_tp **geo_problem, int name, int neigh,
			     double param, double range, int nrow, int ncol, double *prec, int cyclic_flag); 
int GMRFLib_revise_geo_problem(GMRFLib_geo_problem_tp *geo_problem,
			       int name, double param, double range, double *prec);
int GMRFLib_is_geo_coefs(int name, int neigh, double param, double range);
int GMRFLib_print_geo_coefs(FILE *fp);

__END_DECLS

__END_DECLS
#endif
