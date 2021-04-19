#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <mpi.h>

#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"
#include "Readfile.h"


using namespace std;



#define bloc std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
#define SHOW(a) std::cout << #a << std::endl;





void print_vector(std::vector<double> x)
{
  int n=x.size();
  cout<<"le vecteur de taille "<<n<<endl;

  for (int i = 0; i<n; i++)
  {
    //cout<<x[i]<<" ";
  }
  cout<<endl;
  cout<<"----------------------------------"<<endl;

}


void print_matrix(std::vector<std::vector<double>> A)
{
  //cout<<"la diagonale de la matrice:"<<endl;
  //cout<<"[ ";
  for(int i=0; i<A[0].size();i++)
  {
    //cout<<A[0][i]<<" ";
  }
  //cout<<" ]"<<endl;
  //cout<<"la sur diagonale de la matrice:"<<endl;
  //cout<<"[ ";
  for(int i=0; i<A[1].size();i++)
  {
    //cout<<A[1][i]<<" ";
  }
  //cout<<" ]"<<endl;
  //cout<<"la sur-sur diagonale de la matrice:"<<endl;
  //cout<<"[ ";
  for(int i=0; i<A[2].size();i++)
  {
    //cout<<A[2][i]<<" ";
  }
  //cout<<" ]"<<endl;
  //cout<<endl;

}


void print_matrix_verbose(std::vector<std::vector<double>> A)
{
  //cout<<"version verbose de la matrice creuse"<<endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if(A.size()!=3)
  {
    //cout<<"mauvais format de matrice"<<endl;
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
          //cout<<A[0][i]<<" ";
        }else if( (j==i+1) || (j==i-1) ){
          //cout<<A[1][i]<<" ";
        }else if(j==Nx+i){
	  // cout<<A[2][i];
        }else{
          //cout<<"0 ";
        }
      }
      //cout<<endl;
    }
  }
}

vector<int> charge (int n ,int Np, int me )
{
  int limite  =n - Np*(n/Np) ;
  vector <int > res (2) ;

 if ( me < limite)
   {
  res [0] = me*(n/Np+1);
  res [1] = res[0]+n/Np;
   }
  else
    {
      res[0]= limite*(n/Np+1)+(me-limite)*(n/Np) ;
      res[1]= res[0]+(n/Np)-1;
    }
  return res ;
}





int main(int argc, char** argv)
{


  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  if (argc < 2)
  {
    cout << "Ooops, forgot to give the data_file name" << endl;
    exit(0);
  }


  //récuperer le nom
  const string data_file_name = argv[1];
  bloc
  cout<<data_file_name<<endl;
  Readfile* Rf = new Readfile(data_file_name);
  Rf->Read_data_file();
  // ------------------------------------------------------------




  //check the algebra and simple output
  std::vector<double> A(4),b(2);
  A[0]=1.;A[1]=2.;A[2]=3.;A[3]=4.;
  b[0]=5.;b[1]=6.;
  //print_vector(b);





  //Problem test
  //donées du problème

  //double Lx=1.,Ly=1.,D=1.,deltat=0.1,tf=10.;
  //int Nx=40,Ny=40,Nt=1;
  double Lx=Rf->Get_Lx(),Ly=Rf->Get_Ly(),D=Rf->Get_D(),deltat=Rf->Get_dt(),tf=Rf->Get_tfinal();
  int Nx=Rf->Get_Nx(),Ny=Rf->Get_Ny(),Nt=4;  //Nt ou delta t à éliminer





  BC bc=BC(Nx,Ny,Lx,Ly);
  Problem P=Problem(&bc,Nx ,  Ny,  Nt,  Lx,  Ly, deltat);
  std::vector<std::vector<double>> B(3);
  B=P.Construct_Matrix();

  //test of matrix construction (**** checkede)
  //print_matrix(B);


  // operator = is working :(
  std::vector<double> v(4,1.);
  std::vector<double> u(4);
  u=v;
  SHOW(u);
  //print_vector(v);
  SHOW(v);
  //print_vector(u);


  //test od static methods-----------------------------
  std::vector<double> y(Nx*Ny,2.);
  std::vector<double> g(Nx*Ny,4.);
  
  GradConj mc(B,g,Nx,Ny);
  
// test produit parallel 
  
 
  //
  //y=mc.product(B,g,Nx,Ny);
  //SHOW(y);
  //print_vector(y);

  //y=mc.prod_scal(g,5.);
  //SHOW(y);
  //print_vector(y);


  double test(0.);
  test=mc.dot_product(g,y);
  //SHOW(y);
  cout<<"prod scalaire "<<test<<endl;

  test=mc.norm(g);
  //SHOW(y);
   cout<<"norm "<<test<<endl;



  //test of conjugate gradient
  //grad conj B,g
  int state=100; //nb_iter
  mc.Solve(state,y);
  //print_vector(y);

  //vérif
  //print_vector(mc.product(B,y,Nx,Ny));
  //print_vector(g);


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
  int cas=Rf->Get_cas();
  P.Solve_problem(cas,tf); //cas 4
  y=P.get_sol();
  bloc
  print_vector(y);

  //sauvegarde du fichier
  Output io=Output(&P);
  io.Save_sol("sol.dat");
  io.splot_solution("sol.dat");



  //test of parallel region


  /**********************************************************/

 bloc
  //test of sum inside parallel region
  std::vector<double> y1(Nx*Ny,2.);
  std::vector<double> g1(Nx*Ny,4.);
  //print_vector(y1);
  //print_vector(g1);



/*   MPI_Init(&argc,&argv);





  y=mc.MPI_sum(y1,g1,-1);
  //print_vector(y);
  //la somme fonctionne

  double temp=mc.MPI_dot_product(y1,g1);
  bloc
  cout<<temp<<endl;
  //prod scal fonctionne
  bloc
  temp=mc.MPI_norm(y1);
  cout<<temp<<endl;

  //reste à paralléliser le produit et organiser la distribution
  //------------------------------------------------------------
  int me,Np,tag,input,begin,end;
  tag=100;
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  mc.product_parallel(B,g,y,me, Np,Nx,Ny);
  
  MPI_Finalize (); */
  // test produit parallel
  cout << "-----------------------------------------produit parallel ---------------" << endl ;
  print_vector(mc.sum(mc.product(B,g,Nx,Ny), y ,-1));




  bloc 
  bloc 
  MPI_Init(&argc,&argv);

  int me,Np,tag,input,begin,end;
  tag=100;
  std::vector<int> v(2);

  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  v=charge(Nx*Ny,Np,me);
  std::vector<std::vector>double>> C(3);
  int size=v[1]-v[0]+1;
  C[0].resize(size);
  C[1].resize(size);
  C[2].resize(size);

  std::vector<double> x(size,0.);
  



<<<<<<< HEAD

  
=======



  MPI_Finalize();
  //print_vector(y);

  bloc 
  bloc 

  Rf->Assembel_sol_file(3);




>>>>>>> f28755e0cc003106ecb3f3a3dc058d5223b65996


  MPI_Finalize ();



  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  std::cout << "-------------------------------------------------" << std::endl;
  cout << "Cela a pris "<< t << " seconds" << endl;

  return 0;
}
