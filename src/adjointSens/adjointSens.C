/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    adjointSens 

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

#include <iostream>
#include <fstream>
#include "IOmanip.H"

#include <armadillo>
using namespace arma;

#include <omp.h>

void save_scalar_list(const scalarField S, const int rows, std::basic_string<char> inputName);
void save_scalar(const double val, std::basic_string<char> inputName);

int main(int argc, char *argv[])
    {
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"


    Info<< "Start Adjoint part:\n" << endl;    

    #include "p2Eqn.H"
    
    omp_set_num_threads(12);
    const bool is3D = (mesh.nGeometricD() == 3) ? 1 : 0;
    
    #include "OAD_matrices.H"
    #include "OAD_sensitivities.H"        
    
    save_scalar(Fy, runTime.timeName() + "/Fy.dat");
    save_scalar_list(DFy_Dgma.internalField(), mesh.nCells(), runTime.timeName() + "/DFy_Dgma.dat");
    DFy_Dgma.write();
    
    save_scalar(Fx, runTime.timeName() + "/Fx.dat");
    save_scalar_list(DFx_Dgma.internalField(), mesh.nCells(), runTime.timeName() + "/DFx_Dgma.dat");
    DFx_Dgma.write();

    Info<< "End Adjoint part.\n" << endl;
    return 0;
    }
    
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ My Functions:
void save_scalar_list(const scalarField S, const int rows, std::basic_string<char> inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(20);
    for(int i=0; i<rows ; i++)
        {
        writeToDisk << S[i] << std::endl;
        }
    writeToDisk.close();
    }
    
void save_scalar(const double val, std::basic_string<char> inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(16);
    writeToDisk << val << std::endl;
    writeToDisk.close();
    }

