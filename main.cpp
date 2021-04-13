#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
//#include <mpi.h>

#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"


using namespace std;



#define bloc std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
#define SHOW(a) std::cout << #a << std::endl;

/* std::vector<double> product1(std::vector<std::vector<double>> A,std::vector<double> x, int Nx, int Ny)
{
	//ajouter la partie sim (seulement riangle sup * vect)
	std::vector<double> y;
  bool a,b,c,d;
	for(int i=0; i<Nx*Ny; i++)
	{
    a=(i>Nx-1);
    b=(i<Nx*(Ny-1));
    c=(i>0);
    d=(i<Nx*Ny-1);
    y.push_back(A[0][i]*x[i]+d*A[1][i]*x[i+1]+b*A[2][i]*x[i+Nx]+c*A[1][i-1]*x[i-1]+a*A[2][i-Nx]*x[i-Nx]);

	}
  return y;
} */
void print_vector(std::vector<double> x)
{
  int n=x.size();
  cout<<"le vecteur de taille "<<n<<endl;

  for (int i = 0; i<n; i++)
  {
    cout<<x[i]<<" ";
  }
  cout<<endl;
  cout<<"----------------------------------"<<endl;

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




  //check the algebra and simple output
  std::vector<double> A(4),b(2);
  A[0]=1.;A[1]=2.;A[2]=3.;A[3]=4.;
  b[0]=5.;b[1]=6.;
  print_vector(b);





  //Problem test
  //donées du problème
  double Lx=1.,Ly=1.,D=1.,deltat=0.1,tf=10.;
  int Nx=40,Ny=40,Nt=1;





  BC bc=BC(Nx,Ny,Lx,Ly);
  Problem P=Problem(&bc,Nx ,  Ny,  Nt,  Lx,  Ly, deltat);
  std::vector<std::vector<double>> B(3);
  B=P.Construct_Matrix();

  //test of matrix construction (**** checkede)
  print_matrix(B);


  // operator = is working :(
  std::vector<double> v(4,1.);
  std::vector<double> u(4);
  u=v;
  SHOW(u);
  print_vector(v);
  SHOW(v);
  print_vector(u);


  //test od static methods-----------------------------
  std::vector<double> y(Nx*Ny,2.);
  std::vector<double> g(Nx*Ny,4.);

  GradConj mc(B,g,Nx,Ny);

  y=mc.product(B,g,Nx,Ny);
  SHOW(y);
  print_vector(y);

  y=mc.prod_scal(g,5.);
  SHOW(y);
  print_vector(y);


  double test(0.);
  test=mc.dot_product(g,y);
  SHOW(y);
  cout<<"prod scalaire "<<test<<endl;

  test=mc.norm(g);
  SHOW(y);
   cout<<"norm "<<test<<endl;



  //test of conjugate gradient
  //grad conj B,g
  int state=100; //nb_iter
  mc.Solve(state,y);
  print_vector(y);

  //vérif
  print_vector(mc.product(B,y,Nx,Ny));
  print_vector(g);


  //test Bc functions
  test=bc.Source_term(10,20,25,4);
  cout<<"test functions "<<test<<endl;



  //test de construction de F cas 4, fait en dehors de la classe, non initié avec sol ça marche
  /*
  double tt=0.;
  std::vector<double> l;
  P.Construct_F(4,tt,l);
  print_vector(y);
  */

  //implémentation du cas 4 méthode create second term, class output to print and splot
  bloc
  int cas=5;
  P.Solve_problem(cas,tf); //cas 4
  y=P.get_sol();
  bloc
  print_vector(y);

  //sauvegarde du fichier
  Output io=Output(&P);
  io.Save_sol("sol.dat");
  io.splot_solution("sol.dat");


  bloc
  //test of parallel region
//  MPI_Init( &argc, &argv );








  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  std::cout << "-------------------------------------------------" << std::endl;
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}
