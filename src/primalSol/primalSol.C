/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    primalSol

Description
    Steady-state solver for laminar incompressible, using the SIMPLE algorithm.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include <iostream>
#include <fstream>
#include "IOmanip.H"

double DSFv(double var);
void save_list(const Foam::scalarField sF, const char* inputName);
void save_RecMat(const Foam::RectangularMatrix<scalar> Mat, const int rows , const char* inputName);
void save_SqMat(const Foam::SquareMatrix<scalar> Mat, const int rows , const char* inputName);
void save_DiagMat(const Foam::scalarDiagonalMatrix Mat, const char* inputName);
void save_scalar(const double val, std::basic_string<char> inputName); //const char* inputName);
void save_label(const unsigned int val, const char* inputName);
void save_label_list(const List<label> L, const char* inputName);
void save_vector_list(const vectorField Sf, const int rows, const char* inputName);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
    {
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ reading X.dat  	
    ifstream fileReader;
    
    // ~~~~~~~ find the size of X (design varianles)
    fileReader.open(runTime.timeName() + "/X.dat");    
    label sizeOfX = 0;
    while(!fileReader.eof()) 
	    {
        scalar S;
        if(fileReader >> S)  sizeOfX++;	    
	    }
    fileReader.close();	

    // ~~~~~~~ read X.dat and save it in sF_X (design varianles)    
    scalarField sF_X(sizeOfX, 0.0);
    fileReader.open(runTime.timeName() + "/X.dat");
    for(label i = 0; i < sizeOfX; i++)
	    {   
        fileReader >> sF_X[i];
	    }
    fileReader.close();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ assigning X.dat to gamma   	
    // if parallel:
    if(Pstream::parRun())
        {
	    labelIOList localToGlobalID
    	    (
	        IOobject
	            (
		        "cellProcAddressing",
		        mesh.facesInstance(), mesh.meshSubDir, mesh,
		        IOobject::MUST_READ, IOobject::NO_WRITE
	            )
    	    );
        
        forAll(gamma, i)
            {
            gamma[i] = sF_X[localToGlobalID[i]];          
            }
        }
    else    // if serial:
        {
        forAll(gamma, i)
            {
            gamma[i] = sF_X[i];            
            }
        }
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Original SIMPLE loop
    Info<< "Start mySolver:\n" << endl;
    
    while (simple.loop(runTime))
        {
        simple.loop(runTime);
            {
            #include "UEqn.H"
            #include "pEqn.H"
            }
        }
    
    // ~~~~~~~ save the last time name forlder
    save_label(runTime.endTime().value() , "data/cTime.dat");

    scalar FX = sum(mesh.V() * (U.component(vector::X) * gamma))/etta.value();
    scalar FY = sum(mesh.V() * (U.component(vector::Y) * gamma))/etta.value();

    Info << "Drag (Fx) = " << FX << ", and Lift (Fy) = " << FY << endl; 
    Info<< "End mySolver.\n" << endl;    
    return 0;
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ My Functions:
double DSFv(double var)
    {
    double value(0.0);
    if(var >= 0.0)
        {
        value = var;
        }
    else
        {
        value = 0.0;
        }
    return value;
    }
    
void save_list(const Foam::scalarField sF, const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(14);
    for(int i=0; i<sF.size() ; i++)
        {
        writeToDisk << sF[i] << " " << std::endl;
        }
    writeToDisk.close();
    }

void save_RecMat(const Foam::RectangularMatrix<scalar> Mat, const int rows , const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(14);
    for(int i=0; i<rows ; i++)
        {
        for(int j=0; j<Mat.size()/rows ; j++)
            {
            writeToDisk << Mat[i][j] << " ";
            }
        writeToDisk << std::endl;  
        }
    writeToDisk.close();
    }

void save_SqMat(const Foam::SquareMatrix<scalar> Mat, const int rows , const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(14);
    for(int i=0; i<rows ; i++)
        {
        for(int j=0; j<Mat.size()/rows ; j++)
            {
            writeToDisk << Mat[i][j] << " ";
            }
        writeToDisk << std::endl;  
        }
    writeToDisk.close();    
    }

void save_DiagMat(const Foam::scalarDiagonalMatrix Mat, const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(14);
    for(int i=0; i<Mat.size() ; i++)
        {
        writeToDisk << Mat[i] << " " << std::endl;
        }
    writeToDisk.close();
    }
    
void save_scalar(const double val, std::basic_string<char> inputName) //const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk.precision(32);
    writeToDisk << val << std::endl;
    writeToDisk.close();
    }

void save_label(const unsigned int val, const char* inputName)
    {
    ofstream writeToDisk;
    writeToDisk.open(inputName);
    writeToDisk << val << std::endl;
    writeToDisk.close();
    }
    
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
// ************************************************************************* //
