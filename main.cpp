#include <iostream>
#include <fstream>
#include <chrono>
#include "Problem.h"

using namespace std;


int main(int argc, char** argv)
{


    BC BC_functions();
    //donées du problème
    double Lx=1.,Ly=1.,D=1.;
    int N=500;
    Problem Pb(BC* BC_functions);
    

  
  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();


  
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}