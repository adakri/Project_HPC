#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <mpi.h>
#include <vector>


#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"
#include "Readfile.h"


using namespace std;



#define bloc std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
#define SHOW(a) std::cout << #a << std::endl;
#define debug std::cout<<"debug"<<std::endl;





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


void print_matrix(std::vector<std::vector<double> > A)
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
  cout<<"la sous diagonale de la matrice:"<<endl;
  cout<<"[ ";
  for(int i=0; i<A[3].size();i++)
  {
    cout<<A[3][i]<<" ";
  }
  cout<<" ]"<<endl;
  cout<<"la sous-sous diagonale de la matrice:"<<endl;
  cout<<"[ ";
  for(int i=0; i<A[4].size();i++)
  {
    cout<<A[4][i]<<" ";
  }
  cout<<" ]"<<endl;
  cout<<endl;

}


void print_matrix_verbose(std::vector<std::vector<double> > A)
{
  cout<<"version verbose de la matrice creuse"<<endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if(A.size()!=5)
  {
    cout<<"mauvais format de matrice"<<endl;
  }else{
    //retrouvant Nx et Ny
    int Ny=A[0].size()/(A[0].size()-A[2].size());
    int Nx=A[0].size()/Ny;
    cout<<Nx<<" "<<Ny<<endl;
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

  MPI_Init(&argc,&argv);




  //initialisation du parallélisme
  int me,Np,tag,input,begin,end;
  int Nx=3 , Ny=4 ;
  double betax=-1 , betay=-1 , alpha=5;
  tag=100;
  std::vector<int> v(2);
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  
  v=charge(Nx*Ny,Np,me);


  int size=v[1]-v[0]+1;


  std::vector<std::vector<double> > C(5,vector<double>(size));
  std::vector<double> x(size,1.), w(size,2.), prod(size,0.);

  //on souhaite faire une somme de deux vecteurs, suffit d'appeler les fonctions somme sur la partition du vecteur selon charge, de meme pour norme et produit scalaire etc





  //connaite le rang sur la matrice NxNy x NxNy
  int rang=0 ;
  for (int k=0; k<me;k++)
  {
    rang+= charge(Nx*Ny,Np, k )[1]-charge(Nx*Ny,Np, k )[0]+1;
  }


  //construction de la matrice fonctionne
  for (int i=0;i<size ; i++)
  {
    if (rang+i-Nx<0 )
    {
        C[0][i]= 0 ;
      }
      else{
        C[0][i]=betay ;
      }
    if (rang+i+Nx<Nx*Ny)
      {
        C[4][i]=betay ;
      }
      else{
        C[4][i]= 0;
      }
      if (rang+i-1<0 || ((rang+i)%Nx)==0 )
      {
        C[1][i]= 0  ;  }
        else{
          C[1][i]= betax ;
        }
      if (rang+i+1>Nx*Ny || ((rang+i)%Nx)==Nx-1)
        {
          C[3][i]= 0 ;
        }
        else{
          C[3][i]= betax ;
        }
      C[2][i]=alpha ;
  }




  int q , r ;
  bool a,b,c,d ;
  q= Nx/size ;
  r=Nx-q*size ;
  vector < double>  z (Nx), y(Nx) ;


  //me va envoyer ses éléments aux procs qui en a besoin
  MPI_Status Status ;
  for (int k=1 ; k< q+2 ; k++)
  {
    if ( me+k < Np)
        {
          if ( k==q+1)
            {

              MPI_Send (& x[0], r , MPI_DOUBLE , me+k,0,MPI_COMM_WORLD ) ;
              MPI_Recv (& y[(k-1)*size],r , MPI_DOUBLE , me+k , 0 , MPI_COMM_WORLD, & Status );
              }
          else{

              MPI_Send (& x[0],size, MPI_DOUBLE , me+k,Np,MPI_COMM_WORLD ) ;

              MPI_Recv (& y[(k-1)*size],size, MPI_DOUBLE , me+k , 0 , MPI_COMM_WORLD, & Status );
          }
        }
    if ( me-k>-1)

    {


        if ( k==q+1)
      {

        MPI_Send (& x[0], r , MPI_DOUBLE , me-k,0,MPI_COMM_WORLD ) ;

        MPI_Recv (& z[(k-1)*size],r , MPI_DOUBLE , me-k , 0 , MPI_COMM_WORLD, & Status );


      }
        else{


        MPI_Send (& x[0],size, MPI_DOUBLE , me-k,Np,MPI_COMM_WORLD ) ;

        MPI_Recv (& z[(k-1)*size],size , MPI_DOUBLE , me-k , 0, MPI_COMM_WORLD, & Status );
      }

    }
  }


  for (int i=0 ; i < size ; i++ )
    {
      a=(i-Nx>-1);
      b=(i+Nx<size);
      c=(i>0);
      d=(i+1<size);
      //e=(rang+i)%Nx!=0);
    // prod[i]=C[0][i]*x[i]+d*C[1][i]*x[i+1]+(1-d)*C[1][i]*x[i+1]+(1-b)*C[2][i]*x[i+Nx]+b*C[2][i]*y[-size+i+Nx]+c*e*betax*x[i-1]+(1-c)*e*betax*z[0]+a*f*betay*x[i-Nx])+(1-a)*f*betay*z[i-Nx+1] ;

      prod[i]=C[2][i]*x[i]+a*C[0][i]*x[i-Nx]+(1-a)*C[0][i]*z[Nx-i-1]+c*C[1][i]*x[i-1]+(1-c)*C[1][i]*z[0]+b*C[4][i]*x[i+Nx]+(1-b)*C[4][i]*y[Nx+i-size]+d*C[3][i]*x[i+1]+(1-d)*C[3][i]*y[0] ;

    }
    x=prod ;

    bloc

    print_vector(x);

    bloc

MPI_Finalize() ;
  return 0 ;
}
