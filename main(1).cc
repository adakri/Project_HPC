
#include <iostream>
#include <sstream>
#include <string>
#include <cmath> // For Absolute function
#include "matrice.h" // Matrix Class


using namespace std;

int main() {
    
    //                        //
    /* Start of Opening Files */
    //                        //
    
    //----------------------------------//
   /* int* l;
    l=new int[4];
    int* c;
    c=new int[4];
    int* d;
    d=new int[1];*/
   // Matrice *mat=nullptr;
    //mat=new Matrice(4,4,1);
    //Matrice *vect=new Matrice(4,4,1);
    Matrice M(4,4,1);
    Matrice B(4,1,3);
    //mat->print();
    //vect->print();
   // B.print();
   //Matrice B=mat;
   //Matrice V=vect;
   printf("Erreur ici\n");
    Matrice K=M*B ;
    
   // B=(*this)*B;
    K.print();
    //Matrice K=mat*vect;
    //print();
    return 0;
}