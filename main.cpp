#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include "Problem.h"

using namespace std;



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


void print_matrix(std::vector<std::vector<double>> A)
{
  cout<<"la diagonale de la matrice:"<<endl;
  cout<<"[ ";
  for(int i=0; i<A[0].size();i++)
  {
    cout<<A[0][i]<<" ";
  }
  cout<<" ]"<<endl;
  cout<<"la sur diagonale de la matrice:"<<endl;
  cout<<"[ ";
  for(int i=0; i<A[1].size();i++)
  {
    cout<<A[1][i]<<" ";
  }
  cout<<" ]"<<endl;
  cout<<"la sur-sur diagonale de la matrice:"<<endl;
  cout<<"[ ";
  for(int i=0; i<A[2].size();i++)
  {
    cout<<A[2][i]<<" ";
  }
  cout<<" ]"<<endl;
  cout<<endl;

}


void print_matrix_verbose(std::vector<std::vector<double>> A)
{
  cout<<"version verbose de la matrice creuse"<<endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if(A.size()!=3)
  {
    cout<<"mauvais format de matrice"<<endl;
  }else{
    //retrouvant Nx et Ny
    int Ny=A[0].size()/(A[0].size()-A[2].size());
    int Nx=A[0].size()/Ny;
    //cout<<Nx<<" "<<Ny<<endl;
    for(int i=0; i<Nx*Ny; i++)
    {
      for(int j=0; j<Nx*Ny; j++)
      {
        if(i==j)
        {
          cout<<A[0][i]<<" ";
        }else if( (j==i+1) || (j==i-1) ){
          cout<<A[1][i]<<" ";
        }else if(j==Nx+i){
          cout<<A[2][i];
        }else{
          cout<<"0 ";
        }
      }
      cout<<endl;
    }
  }
}



int main(int argc, char** argv)
{

  
  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  //chech the algebra
  std::vector<double> A(4),b(2);
  A[0]=1.;A[1]=2.;A[2]=3.;A[3]=4.;
  b[0]=5.;b[1]=6.;
  GradConj Gc(A,b);

  //output test
  //print_matrix(A);
  print_vector(b);

  //Problem test
  BC* BC_functions;
  //donées du problème
  double Lx=20.,Ly=30.,D=1.,deltat=1.;
  int Nx(2),Ny(6),Nt(10);


  Problem P=Problem(BC_functions, Nx ,  Ny,  Nt,  Lx,  Ly, deltat);
  std::vector<std::vector<double>> B(3);
  B=P.Construct_Matrix();

  print_matrix(B);

  print_matrix_verbose(B);

  
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  std::cout << "-------------------------------------------------" << std::endl;
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}