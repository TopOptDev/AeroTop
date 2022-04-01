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
    preProc

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

void save_label_list(const List<label> L, const char* inputName);
void save_vector_list(const vectorField Sf, const int rows, const char* inputName);
int main(int argc, char *argv[])
    {
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start preprocess:\n" << endl;
    
    List<label> extData(10, 0);
    extData[0] = mesh.nCells(); // NC
    extData[1] = mesh.nFaces(); // NF
    extData[2] = mesh.neighbour().size(); // NFi    
    
    save_label_list(extData, "data/extData.dat");
    save_vector_list(mesh.C(), mesh.C().size(), "data/coords.dat");
    
    Info<< "End preprocess.\n" << endl;
    
    return 0;
    }
// ************************************************************************* //
void save_label_list(const List<label> L, const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    for(int i=0; i<L.size() ; i++)
        {
        writeToDisk << L[i] << " " << std::endl;
        }
    writeToDisk.close();
    }
        
void save_vector_list(const vectorField Sf, const int rows, const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(14);
    for(int i=0; i<rows ; i++)
        {
        writeToDisk << Sf[i].x() << " " << Sf[i].y() << " " << Sf[i].z() << std::endl;
        }
    writeToDisk.close();
    }  
