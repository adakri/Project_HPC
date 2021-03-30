#include<iostream>
#include<vector>
#include<cmath>
using std::vector;
class Matrice
{
    private:
        int _lign;
        int _colone;
        vector<vector<double> > _matrice;
        vector<double> _AA;
        vector<int> _IA;
        vector<int> _JA;
        int _nnz;
    public:
        Matrice(int,int,double);
        //Les Opérations des matrice
        Matrice operator+(Matrice &);
        Matrice operator-(Matrice &);
        Matrice operator*(Matrice &);
        //Matrice operator=(Matrice &);
        //Les opérations matrice scalaire
        Matrice operator+(double &);
        Matrice operator-(double &);
        Matrice operator*(double &);
        //Fonctions utiles
        double& operator()(int &,int &);
        void print();
        //Fonctions pour Matrice creuse (sparse)
        void Creuse();
        vector<double> get_AA(){return _AA;};
        vector<int> get_IA(){return _IA;};
        vector<int> get_JA(){return _JA;};
        int get_nnz(){return _nnz;};

};