#ifndef _GC_CPP

#include <stdlib.h>
#include "GradConj.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include <mpi.h>




#define debug std::cout <<"step here100"<< std::endl;
#define siz(a) std::cout <<a.size()<< std::endl;




using namespace std;

vector<int> charge_c (int n ,int Np, int me )
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


void print_vector1(std::vector<double> x)
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

//constructeur
GradConj::GradConj(std::vector<std::vector<double>> A,std::vector<double> b, int Nx, int Ny) : Nx_(Nx), Ny_(Ny)
{

	A_=A;
	b_=b;
	std::cout<<"Classe GC initié"<<endl;
}

//matrix vector and vector vector manipulations

std::vector<double> GradConj::product(std::vector<std::vector<double>> A,std::vector<double> x, int Nx, int Ny)
{
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
}

std::vector<double> GradConj::sum(std::vector<double> x,std::vector<double> y, int sign) //-1 ou 1
{
	int n=x.size();
	int m=y.size();
	if(m!=n)
	{
		exit(0);
	}else{
		std::vector<double> z(n);
		for(int i=0; i<n; i++)
		{
			z[i]=x[i]+sign*y[i];
		}
		return z;
	}
}


std::vector<double> GradConj::MPI_sum(std::vector<double> x ,std::vector<double> y, int sign)
{
  int n=x.size();
	int m=y.size();
  int root=0;
  if(m!=n)
	{
		exit(0);
	}else{
    //à recontruire
    std::vector<double> z(n);
    std::vector<double> temp;
    std::vector<int> count(2);
    //env parallel
    int me,Np,tag,input,begin,end;
    tag=100;
    MPI_Comm_size(MPI_COMM_WORLD,&Np);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);



    if(me!=root)
    {
      count=charge_c(n,Np,me);
      begin=count[0];
      end=count[1];

      int size(end-begin+1);
      temp.resize(size);
      for(int i=0; i<size; i++)
      {
        temp[i]=x[i]+y[i];
      }

      //bloc vérif
      // siz(temp)
      // print_vector1(temp);
      // cout<<"---------------------"<<me<<" "<<begin<<" "<<end<<endl;

      MPI_Send(&temp[0],size,MPI_INT,root,tag,MPI_COMM_WORLD);
    }else{
      //travail usuel
      count=charge_c(n,Np,me);
      begin=count[0];
      end=count[1];
      int size(end-begin+1);
      temp.resize(size);

      for(int i=0; i<size; i++)
      {
        temp[i]=x[i]+y[i];
      }
      //bloc vérif
      // siz(temp)
      // print_vector1(temp);
      // cout<<"---------------------"<<me<<" "<<begin<<" "<<end<<endl;
      //bloc vérif
      // siz(temp)
      // print_vector1(temp);
      // cout<<"---------------------"<<me<<" "<<begin<<" "<<end<<endl;
      //reconstruction
      MPI_Status status;
      for(int k=0; k<Np; k++)
      {
        count=charge_c(n,Np,k);
        begin=count[0];
        end=count[1];
        int size(end-begin+1);
        temp.resize(size);
        if(k==0)
        {
          //don't do a thing
          //reconstruct
          for(int l=0; l<size; l++)
          {
            z[l]=temp[l];
          }
        }else{
          MPI_Recv(&temp[0],size,MPI_INT,k,tag,MPI_COMM_WORLD,&status);

          for(int l=0; l<size; l++)
          {
            z[begin+l]=temp[l];
          }
        }

      }
    }
    return z;
  }
}















std::vector<double> GradConj::prod_scal(std::vector<double> x,double y)
{
	int n=x.size();
  	std::vector<double> f(n,0);
	for (int i = 0; i <n; i++)
	{
    f[i]=y*x[i];
 	}
  	return f;
}
double GradConj::norm(std::vector<double> x)
{
		int n=x.size();
		double sum=0.;
		for(int i=0; i<n; i++)
		{
			sum+=x[i]*x[i];
		}
		return sqrt(sum);

}
double GradConj::dot_product(std::vector<double> x,std::vector<double> y)
{
	int n=x.size();
	int m=y.size();
	if(m!=n)
	{
		exit(0);
	}else{
		double z(0.);
		for(int i=0; i<n; i++)
		{
			z+=x[i]*y[i];
		}
		return z;
	}
};



//gradient conjugué
void GradConj::Solve(int state,std::vector<double>& u)
{
	//cout<<"**********************GC commences**********************"<<endl;
  	int n = Nx_*Ny_;
	k_=state;
	cout<<"le nombre d'itérations d'entrée "<<k_<<endl;
	std::vector<std::vector<double>> A(A_);
	std::vector<double> r(n),b(b_),p(n),temp(n);
	std::vector<double> x(n);
	for (int i = 0; i < n; i++)
	{
		x[i]=0.;

	}

	//cout<<"le second terme"<<endl;
	//print_vector1(b);

	temp=GradConj::product(A,x,Nx_,Ny_);
	//print_vector1(temp);
	r=GradConj::sum(b,temp,-1);
	//print_vector1(r);
	p= r  ;	// calcul du residu
	double alpha;
	double gamma;
	std::vector<double> rSuivant(n);
	std::vector<double> xSuivant(n);
	std::vector<double> z(n);
	int j = 0;
	double beta=GradConj::norm(r);
	int nb_iterat_=0;

	while (j<=k_)
	{

		//cout<<"________________________loop_____________"<<endl;
		z=GradConj::product(A,p,Nx_,Ny_);
		//print_vector1(z);
		alpha= (GradConj::dot_product(r,r) )  /(GradConj::dot_product(z,p));
		//debug;
		xSuivant=GradConj::sum(x,GradConj::prod_scal(p,alpha),1);
        rSuivant=GradConj::sum(r,GradConj::prod_scal(z,alpha),-1);
		gamma= GradConj::dot_product(rSuivant,rSuivant) /GradConj::dot_product(r,r);
    	p=GradConj::sum(rSuivant,prod_scal(p,gamma),1);
		x=xSuivant;
		//cout<<"----------------------------------------"<<endl;
		r=rSuivant;
		beta=GradConj::norm(r);
		nb_iterat_=nb_iterat_ +1;
		j++;
		if(beta<pow(10,-10))
		{
			break;
		}
	}
	cout<<nb_iterat_<<endl;

	cout<<"----------------Gradient conjugué------------------------"<<endl;
  u.resize(n);
  u=x;
}


//gradient conjugué version parallèle
void GradConj::MPI_Solve(int state,std::vector<double>& u)
{
	//cout<<"**********************GC commences**********************"<<endl;
  	int n = Nx_*Ny_;
	k_=state;
	cout<<"le nombre d'itérations d'entrée "<<k_<<endl;


	std::vector<std::vector<double>> A(A_);
	std::vector<double> r(n),b(b_),p(n),temp(n);
	std::vector<double> x(n);
	for (int i = 0; i < n; i++)
	{
		x[i]=0.;

	}

	//cout<<"le second terme"<<endl;
	//print_vector1(b);

	temp=GradConj::product(A,x,Nx_,Ny_);
	//print_vector1(temp);
	r=GradConj::sum(b,temp,-1);
	//print_vector1(r);
	p= r  ;	// calcul du residu
	double alpha;
	double gamma;
	std::vector<double> rSuivant(n);
	std::vector<double> xSuivant(n);
	std::vector<double> z(n);
	int j = 0;
	double beta=GradConj::norm(r);
	int nb_iterat_=0;

	while (j<=k_)
	{

		//cout<<"________________________loop_____________"<<endl;
		z=GradConj::product(A,p,Nx_,Ny_);
		//print_vector1(z);
		alpha= (GradConj::dot_product(r,r) )  /(GradConj::dot_product(z,p));
		//debug;
		xSuivant=GradConj::sum(x,GradConj::prod_scal(p,alpha),1);
        rSuivant=GradConj::sum(r,GradConj::prod_scal(z,alpha),-1);
		gamma= GradConj::dot_product(rSuivant,rSuivant) /GradConj::dot_product(r,r);
    	p=GradConj::sum(rSuivant,prod_scal(p,gamma),1);
		x=xSuivant;
		//cout<<"----------------------------------------"<<endl;
		r=rSuivant;
		beta=GradConj::norm(r);
		nb_iterat_=nb_iterat_ +1;
		j++;
		if(beta<pow(10,-10))
		{
			break;
		}
	}
	cout<<nb_iterat_<<endl;

	cout<<"----------------Gradient conjugué------------------------"<<endl;
  u.resize(n);
  u=x;
}

#define _GC_CPP
#endif
