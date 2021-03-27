#ifndef _PB_CPP


#include "Problem.h"
#include <fstream>
#include <iostream>

// Constructeur
Problem::Problem(BC* BC_functions) : BC_functions_(BC_functions)
{
	std::cout << "Construction de la classe Matrice" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
std::vector<double> Construct_Matrix(int N, double Lx, double Ly)
{

}


#define _PB_CPP
#endif