#ifndef BCDEF
#define BCDEF


#include <fstream>
#include <iostream>
#include <vector>

#include "Problem.h"




class BC
{
private:
    const Problem* Pb_;
    int Nx,Ny;
    double Lx,Ly;
public:
    //constructor
    BC(Problem* problem);
    //destructor
    ~BC(){};
    double Exact_solution(const double x, const double y, const double t) const;


    double Initial_condition(const double x, const double y, const double t, const int cas) const;


    double Source_term(const double x, const double y, const double t) const;
    double Neumann_Function(const double x, const double y, const double t) const;
    double Dirichlet_Function(const double x, const double y, const double t) const;
};


#endif