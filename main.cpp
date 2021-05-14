#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include "Problem.h"
#include "GradConj.h"
#include "BC.h"
#include "Output.h"
#include "Readfile.h"

using namespace std;

//macros pour la présentation et le debuggage

#define bloc std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
#define SHOW(a) std::cout << #a << std::endl;

//Fonction de réprensetation de vecteurs de la librairie vector.h

void print_vector(std::vector<double> x)
{
  int n = x.size();
  cout << "le vecteur de taille " << n << endl;

  for (int i = 0; i < n; i++)
  {
    cout << x[i] << " ";
  }
  cout << endl;
  cout << "----------------------------------" << endl;
}

//Fonctions de représentation de matrice (vecteur de vecteurs) de forme abrégée (liste de vecteurs) ou verbose (matrice carrée)

void print_matrix(std::vector<std::vector<double>> A)
{
  cout << "la diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[0].size(); i++)
  {
    cout << A[0][i] << " ";
  }
  cout << " ]" << endl;
  cout << "la sur diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[1].size(); i++)
  {
    cout << A[1][i] << " ";
  }
  cout << " ]" << endl;
  cout << "la sur-sur diagonale de la matrice:" << endl;
  cout << "[ ";
  for (int i = 0; i < A[2].size(); i++)
  {
    cout << A[2][i] << " ";
  }
  cout << " ]" << endl;
  cout << endl;
}
void print_matrix_verbose(std::vector<std::vector<double>> A)
{
  cout << "version verbose de la matrice creuse" << endl;
  std::cout << "-------------------------------------------------" << std::endl;
  if (A.size() != 3)
  {
    cout << "mauvais format de matrice" << endl;
  }
  else
  {
    //Retrouvant Nx et Ny
    int Ny = A[0].size() / (A[0].size() - A[2].size());
    int Nx = A[0].size() / Ny;
    cout << Nx << " " << Ny << endl;
    for (int i = 0; i < Nx * Ny; i++)
    {
      for (int j = 0; j < Nx * Ny; j++)
      {
        if (i == j)
        {
          cout << A[0][i] << " ";
        }
        else if ((j == i + 1) || (j == i - 1))
        {
          cout << A[1][i] << " ";
        }
        else if (j == Nx + i)
        {
          cout << A[2][i];
        }
        else
        {
          cout << "0 ";
        }
      }
      cout << endl;
    }
  }
}

//fonctions d'équilibrage de charges
vector<int> charge(int n, int Np, int me)
{
  int limite = n - Np * (n / Np);
  vector<int> res(2);

  if (me < limite)
  {
    res[0] = me * (n / Np + 1);
    res[1] = res[0] + n / Np;
  }
  else
  {
    res[0] = limite * (n / Np + 1) + (me - limite) * (n / Np);
    res[1] = res[0] + (n / Np) - 1;
  }
  return res;
}

int main(int argc, char **argv)
{

  // Démarrage du chrono
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

  if (argc < 2)
  {
    cout << "Ooops, forgot to give the data_file name" << endl;
    exit(0);
  }

  //Récuperation le nom du fichier d'entrée
  const string data_file_name = argv[1];
  bloc
  cout<< data_file_name << endl;
  Readfile *Rf = new Readfile(data_file_name);
  Rf->Read_data_file();
  // ------------------------------------------------------------

  //Récupération des données du problème
  double Lx = Rf->Get_Lx(), Ly = Rf->Get_Ly(), D = Rf->Get_D(), deltat = Rf->Get_dt(), tf = Rf->Get_tfinal();
  int Nx = Rf->Get_Nx(), Ny = Rf->Get_Ny(), Nt = 4;
  //Nt ou delta t jouent le meme role

  BC bc = BC(Nx, Ny, Lx, Ly);
  Problem P = Problem(&bc, Nx, Ny, Nt, Lx, Ly, deltat);
  std::vector<std::vector<double>> B(3);
  B = P.Construct_Matrix();

  //Initialisation des différentes classes nécessaires à la résolution dans la classe problem

  //Vecteur temp et représentation
  std::vector<double> y(Nx * Ny, 2.);
  std::vector<double> g(Nx * Ny, 4.);

  GradConj mc(B, g, Nx, Ny);

  //Résolution du problème cas
  bloc int cas = Rf->Get_cas();
  P.Solve_problem(cas, tf);
  y = P.get_sol();
  bloc
  print_vector(y);

  //sauvegarde du fichier
  Output io = Output(&P);
  io.Save_sol("sol.dat");
  io.splot_solution("sol.dat");

  // Fin du chrono
  chrono::high_resolution_clock::time_point finish = chrono::high_resolution_clock::now();

  chrono::duration<double> t = chrono::duration_cast<chrono::duration<double>>(finish - start);
  // Affichage du résultat
  std::cout << "-------------------------------------------------" << std::endl;
  cout << "Cela a pris " << t.count() << " seconds" << endl;

  return 0;
}
