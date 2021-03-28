#include <iostream>
#include <fstream>
#include <chrono>

#include <math.h>
#include "Problem.h"

using namespace std;


void print_matrix(std::vector<double> A)
{
  int n=sqrt(A.size());
  cout<<"la matrice de taille "<<n<<"x"<<n<<endl;
  for (int i = n; i --> 0; )
  {
    for (int j = n; j --> 0; )
    {
      cout<<A[i*n+j]<<" ";
    }
    cout<<endl;
  }
}

void print_vector(std::vector<double> x)
{
  int n=x.size();
  cout<<"le vecteur de taille "<<n<<endl;
  for (int i = 0; i<n; i++)
  {
    cout<<x[i]<<" ";
  }
  cout<<endl;
 
}



int main(int argc, char** argv)
{


    BC BC_functions();
    //donées du problème
    double Lx=1.,Ly=1.,D=1.;
    int N=500;
    Problem Pb(BC* BC_functions);


  
    // Démarrage du chrono
    auto start = chrono::high_resolution_clock::now();

    //chech the algebra
    std::vector<double> A(4),b(2);
    A[0]=1.;A[1]=2.;A[2]=3.;A[3]=4.;
    b[0]=5.;b[1]=6.;
    GradConj Gc(A,b);

    //output test
    print_matrix(A);
    print_vector(b);

  
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}