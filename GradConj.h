#ifndef GRADCONJ_H
#define GRADCONJ_H

#include <vector>
#include <fstream>

class GradConj
{
private:

	std::vector<double> A_;
	std::vector<double> b_;
	int k_;
	

public:
	GradConj(std::vector<double> A ,std::vector<double> b) : A_(A), b_(b) { };
	std::vector<double> product(std::vector<double>,std::vector<double>) const;
	std::vector<double> dot_product(std::vector<double>,std::vector<double>) const;
	std::vector<double> GradConj::sum(std::vector<double> ,std::vector<double> y, int sign) const;
	double norm(std::vector<double>) const;
	std::vector<double> Solve(double&)const;
};


#endif