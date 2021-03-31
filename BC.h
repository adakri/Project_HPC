#ifndef BCDEF
#define BCDEF


#include <fstream>
#include <iostream>
#include <vector>

#include "Problem.h"




class BC
{
private:
    int Nx_,Ny_;
    double Lx_,Ly_;
public:
    //constructor
    BC(int Nx, int Ny, double Lx, double Ly);
    //destructor
    ~BC(){};
    double Exact_solution(const double x, const double y, const double t) const;


    double Source_term(const double x, const double y, const double t, const int cas) const;


    double Initial_condition(const double x, const double y, const double t) const;
    double Neumann_Function(const double x, const double y, const double t) const;
    double Dirichlet_Function(const double x, const double y, const double t) const;
};


#endif