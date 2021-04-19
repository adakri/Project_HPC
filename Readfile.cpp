#ifndef _DATA_FILE_CPP

#include "Readfile.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <regex>


using namespace std;

#define huzzah cout<<"huzzah"<<endl;

Readfile::Readfile(std::string file_name)
: file_name_(file_name),  if_tfinal_(false), if_dt_(false),if_results_(false), if_D_(false), if_Lx_(false), if_Ly_(false),if_Nx_(false), if_Ny_(false)
{
}

void Readfile::Read_data_file()
{
    //test of opening
    ifstream data_file(file_name_.data());
    if (!data_file.is_open())
    {
        cout << "Unable to open file " << file_name_ << endl;
        exit(0);
    }
    else
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Reading data file " << file_name_ << endl;
    }

    string file_line;
    //while not empty
    while (!data_file.eof())
    {
        getline(data_file, file_line);

        if (file_line.find("#") !=std::string::npos)
        {
        // Ignore this line (comment)
        }
        else
        {
        if (file_line.find("tf") != std::string::npos)
        {
            huzzah
            data_file >> tfinal_; 
            if_tfinal_ = true;
        }
        if (file_line.find("dt") != std::string::npos)
        {
            data_file >> dt_;
            if_dt_ = true;
        }
        if (file_line.find("Nx") != std::string::npos)
        {
            data_file >> Nx_;
            if_Nx_ = true;
        }

        if (file_line.find("Ny") != std::string::npos)
        {
            data_file >> Ny_;
            if_Ny_ = true;
        }

        if (file_line.find("Lx") != std::string::npos)
        {
            data_file >> Lx_;
            if_Lx_ = true;
        }

        if (file_line.find("Ly") != std::string::npos)
        {
            data_file >> Ly_;
            if_Ly_ = true;
        }

        if (file_line.find("D") != std::string::npos)
        {
            data_file >> D_;
            if_D_ = true;
        }
        if (file_line.find("cas") != std::string::npos)
        {
            data_file >> cas_;
            if_cas_ = true;
        }
        }

    }


    if (!if_D_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        D_ = 0.;
    }
    if (!if_tfinal_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        tfinal_ = 0.;
    }
    if (!if_dt_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        dt_ = 0.;
    }

    if (!if_Lx_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        Lx_ = 0.;
    }

    if (!if_Ly_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        Ly_ = 0.;
    }

    if (!if_Nx_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        Nx_ = 0.;
    }

    if (!if_Ny_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with 0" << endl;
        Ny_ = 0.;
    }
    if (!if_cas_)
    {
        cout << "-------------------------------------------------" << endl;
        cout << "Ooooops, initialized the value with (4)" << endl;
        cas_ = 4;
    }
    
}

void Readfile::Assembel_sol_file(int Np)
{
    ofstream myfile;
    myfile.open ("global_sol.txt");
    for(int i=0; i<Np; i++)
    {
        //test of opening
        string file_name="sol00"+std::to_string(i)+".txt";
        ifstream data_file(file_name.data());
        if (!data_file.is_open())
        {
            cout << "Unable to open file " << file_name << endl;
            exit(0);
        }
        else
        {
            cout << "-------------------------------------------------" << endl;
            cout << "Reading data file " << file_name << endl;
        }
        
        while(!data_file.eof())
        {
            string line;
            getline(data_file, line);
            if(data_file.eof()==false)
            {
                myfile<<line<<"\n";  

            }else{
                myfile<<line;  
            }
        }
        data_file.close();
    }
    myfile.close();

}


#define _DATA_FILE_CPP
#endif