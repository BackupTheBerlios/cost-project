/*!
  \file caa_Rextensions.c
  \brief Routines for writing R extensions
*/
#include <R.h>
#include <Rdefines.h>
#include "caa.h"

/*!
  \brief Get element in R object with name i_str

  \author Hanne Rognebakke
*/
SEXP getListElement(SEXP i_list, char *i_str)
{
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(i_list, R_NamesSymbol);
  int i, inList = 0;
  
  for (i = 0; i < length(i_list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), i_str) == 0) 
      {
	elmt = VECTOR_ELT(i_list, i);
	inList = 1;
	break;
      }
  if(inList == 0)
    printf("WARNING: %s not in list\n",i_str);
  
  return elmt;
}				/* end of getListElement */


