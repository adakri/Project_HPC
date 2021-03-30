#include "matrice.h"
#include<vector>
#include<iostream>
#include<cmath>
using namespace std;

Matrice::Matrice(int lign,int collon,double init)
{
    _lign=lign;
    _colone=collon;
    _matrice.resize(lign);
    for (int i = 0; i < _matrice.size(); i++)
    {
        _matrice[i].resize(collon,init);
    }
    _nnz=0;
}
Matrice Matrice::operator+(Matrice &B)
{
    Matrice somme(_lign,_colone,0.0);
    for (int i = 0; i < _lign; i++)
    {
        for (int j = 0; j < _colone; j++)
        {
            somme(i,j)=this->_matrice[i][j]+B(i,j);
        }

    }

    return somme;
}
Matrice Matrice::operator=(Matrice &B)
{
    for (int i = 0; i < _lign; i++)
    {
        for (int j = 0; j < _colone; j++)
        {
            B(i,j)=this->_matrice[i][j];
        }
        return B;
    }

}
Matrice Matrice::operator-(Matrice &B)
{
    Matrice difference(_lign,_colone,0.0);
    for (int i = 0; i < _lign; i++)
    {
        for (int j = 0; j < _colone; j++)
        {
            difference(i,j)=this->_matrice[i][j]-B(i,j);
        }

    }

    return difference;
}
Matrice Matrice::operator*(Matrice &B)
{
    Matrice K(_lign,_colone,0.0);

    for (int i = 0; i < _lign; i++)
    {
        for (int j = 0; j < _colone; j++)
        {
            K(i,j)=this->_matrice[i][j];
        }
    }

    Matrice prod_CSR(_lign,1,0.0);

    int a=0;
    for (int i = 0; i < _lign; i++)
    {
        int start=K.get_IA()[i];
        printf("NOP\n");
        int end=K.get_IA()[i+1];
        printf("NOP\n");
        for (int j = start; j <end; j++)
        {
            printf("NOP\n");
            prod_CSR(i,a)=prod_CSR(i,a)+K.get_AA()[j]*B(K.get_JA()[j],a);
        }

    }

    return prod_CSR;
}

double& Matrice::operator()(int &lign,int &colone)
{
    return this->_matrice[lign][colone];
}
void Matrice::print()
{
    cout << "Matrice: " << endl;
    for (int i = 0; i < _lign; i++) {
        for (int j = 0; j <_colone; j++) {
            cout << "[" << _matrice[i][j] << "] ";
        }
        cout << endl;
    }
}
void Matrice::Creuse()
{
   /* int nnz;
    printf("NOP\n");
    vector<double> AA(nnz);
    printf("NOP\n");
    vector<int> IA(_lign),JA(nnz);
    IA[0]=0.0;
    printf("NOP\n");*/

    for (int i = 0; i < _lign; i++) {
        for (int j = 0; j <_colone; j++)
        {
            printf("NOP1\n");
           if(this->_matrice[i][j]!=0)
           {
                printf("NOP2\n");
                _AA.push_back(this->_matrice[i][j]);
                printf("NOP3\n");
                _JA.push_back(j);
                _nnz++;
               // cout<<"AA  "<<this->_matrice[i][j]<<endl;
                /*nnz++;
                AA[nnz-1]=this->_matrice[i][j];
                JA[nnz-1]=j;*/
           }
        }
        _IA[i]=_nnz;
    }
    /*_AA=AA;
    _IA=IA;
    _JA=JA;
    _nnz=nnz;*/


}
