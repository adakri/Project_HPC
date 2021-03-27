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
    int N_,Nt_;
    double Lx_,Ly_,t;
    const double D=1.;
    std::vector<double> A_,F_;
    BC* BC_functions_;

public:
    //constructor
    Problem(BC* BC_functions);
    //destructor
    ~Problem();
    std::vector<double> Construct_Matrix(int, double, double){return A_;};
    std::vector<double> Construct_F(int, double, double){return F_;};
};


#endif