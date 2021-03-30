#include<vector>
#include <iostream>
#include<fstream>
#include<cmath>
#include <string>

using namespace std;

int k(int i, int j, int n)
{
    int t=i*n+j;
    return t;
}

void Build_flux_mat_and_rhs()
{
    int nx=4;
    int ny=3;
    int n=nx*ny;
    double Lx=1., Ly=1.;
    double D=1.;
    double dt=1;
    double dx=Lx/(nx+1.);
    double dy=Ly/(ny+1.);
	vector<vector<double>> _mat_flux(n,vector<double>(n));
    vector<double> _BC_RHS(n,0);
    
	for (int i = 0; i <nx; i++)
	{
		for (int j = 0; j <ny; j++)
	    {
            int rw=k(i,j,nx);

            if(i==0)
            {
                if (j==0)
                {

                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i+1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j+1,ny)]=0.;

                    _BC_RHS[k(i,j,ny)]=0.;
                }

                if (j==ny-1)
                {
                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i+1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j-1,ny)]=0.;

                    _BC_RHS[k(i,j,ny)]=0.;
                }

                if (j!=0 && j!=ny-1)
                {
                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i+1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j+1,ny)]=0.;
                    _mat_flux[rw,k(i,j-1,ny)]=0.;

                    _BC_RHS[k(i,j,ny)]=0.;
                }
            }
            if(i==nx-1)
            {
                if (j==ny-1)
                {
                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i-1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j-1,ny)]=0.;
            
                    _BC_RHS[k(i,j,ny)]=0.;
                }
                
                if (j==0)
                {
                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i-1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j+1,ny)]=0.;
                    
            
                    _BC_RHS[k(i,j,ny)]=0.;
                }

                if (j!=0 && j!=ny-1)
                {
                    _mat_flux[rw,k(i,j,ny)]=0.;
                    _mat_flux[rw,k(i-1,j,ny)]=0.;
                    _mat_flux[rw,k(i,j+1,ny)]=0.;
                    _mat_flux[rw,k(i,j-1,ny)]=0.;

                    _BC_RHS[k(i,j,ny)]=0.;
                }
            }

            if(j==0 && i!=0 && i!=nx-1)
            {
                 
                _mat_flux[rw,k(i,j,ny)]=0.;
                _mat_flux[rw,k(i+1,j,ny)]=0.;
                _mat_flux[rw,k(i-1,j,ny)]=0.;
                _mat_flux[rw,k(i,j+1,ny)]=0.;
                
                _BC_RHS[k(i,j,ny)]=0.;
                              
            }
            
            if(j==ny-1 && i!=0 && i!=nx-1)
            {
                _mat_flux[rw,k(i,j,ny)]=0.;
                _mat_flux[rw,k(i+1,j,ny)]=0.;
                _mat_flux[rw,k(i-1,j,ny)]=0.;
                _mat_flux[rw,k(i,j-1,ny)]=0.;   
            
                _BC_RHS[k(i,j,ny)]=0.;
             }
            
            if(i!=0 && i!=nx-1 && j!=0 && j!=ny-1)
            {
                _mat_flux[rw,k(i,j,ny)]=1+2*D*dt*(1/(dx*dx)+1/(dy*dy));
                _mat_flux[rw,k(i+1,j,ny)]=-dt*D/(dx*dx);
                _mat_flux[rw,k(i-1,j,ny)]=-dt*D/(dx*dx);
                _mat_flux[rw,k(i,j+1,ny)]=-dt*D/(dy*dy);
                _mat_flux[rw,k(i,j-1,ny)]=-dt*D/(dy*dy);
            }

        }
                      
	} 
    for (int i = 0; i <nx*ny; i++)
	{
		for (int j = 0; j <ny*nx; j++)
	    {
            cout<<" [ "<<  _mat_flux[i,j] <<" ] ";
        }
        cout<<endl;
    }   
}
int main()
{

    return 0;
}