#ifndef PBDEF
#define PBDEF


#include <fstream>
#include <iostream>
#include <vector>

#include "BC.h"
#include "GradConj.h"


class Problem
{
private:
    int Nx_,Ny_,Nt_;
    double Lx_,Ly_,t;
    double deltax_;
    double deltay_;
    double deltat_;
    const double D_=1.;
    std::vector<double> F_;
    std::vector<std::vector<double>> A_;

    BC* BC_functions_;

public:
    //constructor
    Problem(BC*, int , int , int , double , double, double);
    //destructor
    ~Problem(){};
    std::vector<std::vector<double>> Construct_Matrix();
    std::vector<double> Construct_F(int, double, double);
};


#endif