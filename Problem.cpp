#ifndef _PB_CPP


#include <fstream>
#include <iostream>

#include "GradConj.h"
#include "Problem.h"




#define debug std::cout <<"step here"<< std::endl;

#define siz(a) std::cout <<a.size()<< std::endl;




void print_vector2(std::vector<double> x)
{
  int n=x.size();
  std::cout<<"le vecteur de taille "<<n<<std::endl;
  
  for (int i = 0; i<n; i++)
  {
    std::cout<<x[i]<<" ";
  }
  std::cout<<std::endl;
  std::cout<<"----------------------------------"<<std::endl;

}



// Constructeur
Problem::Problem(BC* f,int Nx , int Ny, int Nt, double Lx, double Ly, double deltat) : functions_(f)
{
	Lx_=Lx;
	Ly_=Ly;
	Nx_=Nx;
	Ny_=Ny;
	Nt_=Nt;
	deltat_=deltat;
	std::cout << "Construction de la classe Matrice" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
std::vector<std::vector<double>> Problem::Construct_Matrix()
{
	//la matrice serai 3 vecteurs diag,sdiag et ssdiag
	std::vector<std::vector<double>> A(3);
	//std::vector<double> diag(Nx_*Ny_),sdiag(Nx_*Ny_),ssdiag(Nx_*Ny_);
/* 	A[0].resize(Ny_*Nx_);
	A[1].resize(Ny_*Nx_-2);
	A[2].resize(Nx_*(Ny_-1)); */

	deltax_=Lx_/(Nx_+1.);
	deltay_=Ly_/(Ny_+1.);
	//std::cout<<"1)"<<deltax_<<std::endl;
	//std::cout<<"2)"<<deltay_<<std::endl;

	double phix=-D_*deltat_/(deltax_*deltax_);
	double phiy=-D_*deltat_/(deltay_*deltay_);
	double alpha=1-2*(phix+phiy);

	/* std::cout<<"1"<<phix<<std::endl;
	std::cout<<"2"<<phiy<<std::endl;
	std::cout<<"3"<<alpha<<std::endl; */

	for(int i=0; i<Nx_*Ny_; i++)
	{
		A[0].push_back(alpha);
	}

	for(int i=0; i<Nx_*Ny_-1;i++)
	{
		if((i+1)%Nx_==0)
		{
			A[1].push_back(0.);
		}else{
			A[1].push_back(phix);
		}
	}

	for(int i=0; i<Nx_*(Ny_-1);i++)
	{
		A[2].push_back(phiy);
	}
	A_=A;/////////
	return A;
}
void Problem::Construct_F(int cas, double t, std::vector<double>& test) //sol_==un to be used
{
	int n=Nx_*Ny_;
	F_.resize(n);
	//siz(F_);
	std::vector<double> z(n);
	for(int i=0; i<Ny_; i++)
	{
		for(int j=0; j<Nx_; j++)
		{
			F_[i]=(sol_[i*Nx_+j]+deltat_ * functions_->Source_term(j*deltay_,i*deltax_,t+deltat_,cas)); //to change
		}		
	}
	//siz(F_)
	test.resize(n);
	test=F_;
}

void Problem::Solve_problem(int cas, double tf)
{
	//initialise u_ par fonction initial etc
	int n=Nx_*Ny_;
	sol_.resize(n);
	for(int i=0; i<n; i++)
	{
		sol_[i]=functions_->Source_term(i,i,i,cas); //pour l'instant ça signifie rien
	}

	//solve Au_=f with gradconj::solve after initialsie or not (static)
	//do as much as needed
	int nb_iter=100;
	double t(0.);
	std::vector<double> test;
	while(t<tf)
	{
		A_=Problem::Construct_Matrix(); //writes twice to change
		Problem::Construct_F(cas,t,test);
		GradConj gc=GradConj(A_,F_,Nx_,Ny_);
		gc.Solve(nb_iter,sol_);
		t+=deltat_;
	}
	print_vector2(sol_);
}

std::vector<double> Problem::get_sol()
{
	std::vector<double> x(sol_.size());
	x=sol_;
	return sol_;
}

#define _PB_CPP
#endif
