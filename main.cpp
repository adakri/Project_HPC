#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include "Problem.h"
#include "GradConj.h"
using namespace std;


std::vector<double> product1(std::vector<std::vector<double>> A,std::vector<double> x, int Nx, int Ny)
{
	//ajouter la partie sim (seulement riangle sup * vect)
	std::vector<double> y;
  bool a,b,c,d;
	for(int i=0; i<Nx*Ny; i++)
	{
		/*if(i>Nx*(Ny-1))
		{
			y.push_back(A[0][i]*x[(i%Nx)* Nx + i]+A[1][i]*x[(i%Nx)* Nx + i + 1 ]);
		}else if(i==Nx*Ny-1)
		{
			y.push_back(A[0][i]*x[(i%Nx) * Nx + i]);
		}/*else{
			y.push_back(A[0][i]*x[(i%Nx)* Nx + i]+A[1][i]*x[(i%Nx)* Nx + i + 1 ]+A[2][i]*x[(i%Nx)* (Nx+1) + i]);
		}*/
    a=(i>Nx-1);
    b=(i<Nx*(Ny-1));
    c=(i>0);
    d=(i<Nx*Ny-1);
    y.push_back(A[0][i]*x[i]+d*A[1][i]*x[i+1]+b*A[2][i]*x[i+Nx]+c*A[1][i-1]*x[i-1]+a*A[2][i-Nx]*x[i-Nx]);

	}
  return y;
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
//  GradConj Gc(A,b);

  //output test
  //print_matrix(A);
  print_vector(b);

  //Problem test
  BC* BC_functions;
  //donées du problème
  double Lx=5.,Ly=4.,D=1.,deltat=1.;
  int Nx=4,Ny=3,Nt=1;


  Problem P=Problem(BC_functions, Nx ,  Ny,  Nt,  Lx,  Ly, deltat);
  std::vector<std::vector<double>> B(3);
  B=P.Construct_Matrix();

  print_matrix(B);

  //print_matrix_verbose(B);
//  GradConj mc(A,b);
  std::vector<double> y(Nx*Ny);
  std::vector<double> g(Nx*Ny,1.);
  y=product1(B,g,Nx,Ny);
  print_vector(y);
  GradConj mc(B,y);
  double state=0.001;
  mc.Solve(state,g);
  print_vector(g);
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  std::cout << "-------------------------------------------------" << std::endl;
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}
