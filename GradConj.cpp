#ifndef _GC_CPP


#include "GradConj.h"
#include <iostream>
#include <math.h>
#include <iostream>



 



using namespace std;

//constructeur
GradConj::GradConj(std::vector<double> A ,std::vector<double> b): A_(A), b_(b)
{	
}

//matrix vector and vector vector manipulations
std::vector<double> GradConj::product(std::vector<double> A,std::vector<double> x) const
{


}

std::vector<double> GradConj::sum(std::vector<double> x,std::vector<double> y, int sign) const
{


}

double GradConj::norm(std::vector<double> x) const
{


}
std::vector<double> GradConj::dot_product(std::vector<double>,std::vector<double>) const
{

};



//gradient conjugué
std::vector<double> GradConj::Solve(double& state)const								
{
	
	std::vector<double> A=A_;
	int n = sqrt(A.size());
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
		z=A*p;
		alpha= (r.dot(r) )  / (z.dot(p)) ;
		xSuivant= x + alpha*p;
        rSuivant=r-alpha*z;
		gamma= (rSuivant.dot(rSuivant))/(r.dot(r));
		p= rSuivant+ gamma*p;
		x=xSuivant ;
		//cout<<x<<endl;
		//cout<<"----------------------------------------"<<endl;
		r=rSuivant;
		beta=r.norm();
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
	return x;
}
	
#define _GC_CPP
#endif