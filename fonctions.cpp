#include <iostream> 
#include <math.h>
#include <string.h > 

using namespace ::std ;

double f ( double x , double y ,double t ,  double Lx , double Ly ,int cas ) 
{
  if ( cas==4 )
    return 2*(y-y*y+x-x*x);
  
  if ( cas==5)
    return sin(x)+cos(y);
  
  if ( cas ==6 ) 
    return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(x-Ly/2))*cos( pi*t/2);
    }
double g ( double x , double y, int cas  )
{
  if ( cas==4)
    return 0 ;
  if ( cas==5) 
    return sin(x)+cos(y);
  if ( cas==6)
    return 0 ;
      }
double h ( double x , double y , int cas )
{
  if ( cas==4 )
    return 0 ; 
  if ( cas==5 ) 
    return sin(x)+cos(y);
  if ( cas==6)
    return 1 ;
}
