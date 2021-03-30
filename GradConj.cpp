#ifndef _GC_CPP

#include <stdlib.h>
#include "GradConj.h"
#include <iostream>
#include <math.h>
#include <iostream>







using namespace std;

//constructeur
GradConj::GradConj(std::vector<std::vector<double>> A,std::vector<double> b): //A_(A), b_(b)
{
  int n1=A[0].size();
  int n2=A[1].size();
  int n3=A[2].size();
  int l=b.size();
  for{int i=0;i<n1;i++}
  {
    A_[0][i]=A[0][i];
  }
  for{int i=0;i<n2;i++}
  {
    A_[1][i]=A[1][i];
  }
  for{int i=0;i<n3;i++}
  {
    A_[2][i]=A[2][i];
  }
  for{int i=0;i<l;i++}
  {
    b_[i]=b[i];
  }
	cout<<"Classe GC initié"<<endl;
}

//matrix vector and vector vector manipulations
void GradConj::operator_vec(std::vector<double>&A,std::vector<double>B)
{
  for{int i=0;i<A.size();i++}
  {
    A[i]=B[i];
  }
}
void GradConj::operator_mat(std::vector<std::vector<double>>&A,std::vector<std::vector<double>>B)
{
  for{int i=0;i<A.size();i++}
  {
    operator_vec(A[i],B[i]);
  }
}

std::vector<double> GradConj::product(std::vector<std::vector<double>> A,std::vector<double> x, int Nx, int Ny) const
{
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

std::vector<double> GradConj::sum(std::vector<double> x,std::vector<double> y, int sign) const
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
			z.push_back(x[i]+sign*y[i]);
		}
		return z;
	}

}
std::vector<double> GradConj::prod_scal(std::vector<double> x,double y) const
{
	int n=x.size();
  std::vector<double> f(n,0);
	for (int i = 0; i <n; i++) {
    f[i]=y*x[i];
  }
}
double GradConj::norm(std::vector<double> x) const
{
		int n=x.size();
		double sum=0.;
		for(int i=0; i<n; i++)
		{
			sum+=x[i]*x[i];
		}
		return sqrt(sum);

}
std::vector<double> GradConj::dot_product(std::vector<double> x,std::vector<double> y) const
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
			z.push_back(x[i]*y[i]);
		}
		return z;
	}
};



//gradient conjugué
void GradConj::Solve(double& state,std::vector<double> d)const
{

	std::vector<std::vector<double>> A;
  operator_mat(A,A_);
	//int n = sqrt(A.size());
  int n=Nx_*Ny_;
	std::vector<double> x(n);
	for (int i = 0; i < n; i++)
	{
		x[i]=0.;

	}

	std::vector<double> r(n),b(b_),p(n),temp(n);
	//cout<<b<<endl;

	temp=GradConj::product(A,x);
	r=GradConj::sum(b,temp,-1);
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
		operator_mat(z,GradConj::product(A,p));
		alpha= (GradConj::dot_product(r,r) )  / (GradConj::dot_product(z,p)) ;
		operator_mat(xSuivant,GradConj::sum(x,alpha*p));
        operator_mat(temp,GradConj::prod_scal(p,alpha));
        operator_mat(rSuivant,GradConj::sum(r,temp,-1));
		gamma= GradConj::dot_product(rSuivant,rSuivant))/(GradConj::dot_product(r,r));
    operator_mat(temp,GradConj::prod_scal(p,gamma));
    operator_mat(p,GradConj::sum(rSuivant,temp));
		operator_mat(x,xSuivant);
		//cout<<x<<endl;
		//cout<<"----------------------------------------"<<endl;
		operator_mat(r,rSuivant);
		beta=GradConj::norm(r);
		nb_iterat_=nb_iterat_ +1;
		j++;
			if(beta<pow(10,-100))
			{
				break;
			}
	}
	state=beta;
	cout<<nb_iterat_<<endl;

	cout<<"----------------Gradient conjugué------------------------"<<endl;
  d.resize(n);
  operator_vec(d,x);
}

#define _GC_CPP
#endif
