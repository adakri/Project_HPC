#ifndef BCDEF
#define BCDEF


#include <fstream>
#include <iostream>
#include <vector>



class BC
{
private:

public:
    //constructor
    BC();
    //destructor
    ~BC();
    double Exact_solution(const double x, const double y, const double t) const;
    double Initial_condition(const double x, const double y) const;
    double Source_term(const double x, const double y, const double t) const;
    double Velocity_x(const double x, const double y, const double t) const;
    double Velocity_y(const double x, const double y, const double t) const;
    double Neumann_Function(const double x, const double y, const double t) const;
    double Dirichlet_Function(const double x, const double y, const double t) const;
};


#endif