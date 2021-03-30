#ifndef PBDEF
#define PBDEF


#include <fstream>
#include <iostream>
#include <vector>

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


public:
    //constructor
    Problem(int , int , int , double , double, double);
    //destructor
    ~Problem(){};
    std::vector<std::vector<double>> Construct_Matrix();
    std::vector<double> Construct_F(int, double, double);
    double get_Lx() const {return Lx_;} ;
    double get_Ly() const {return Ly_;} ;
    double get_Nx() const {return Nx_;} ;
    double get_Ny() const {return Ny_;} ;

};


#endif