#ifndef _BC_CPP

#include "BC.h"
#include <cmath>

BC::BC() 
{
  std::cout << "La lecture des fonctions relatives au problème" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}


double BC::Initial_condition(const double x, const double y) const
{
  return 0.; 
}


double BC::Source_term(const double x, const double y, const double t) const
{
  return 0.;
  
}

double BC::Neumann_Function(const double x, const double y, const double t) const
{
  return 0.;
 
}

double BC::Dirichlet_Function(const double x, const double y, const double t) const
{
  return 0.;
 
}

double BC::Exact_solution(const double x, const double y, const double t) const
{
  return 0.;
}


#define _BC_CPP
#endif