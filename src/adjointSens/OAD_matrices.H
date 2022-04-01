Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;
Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Starting_OAD_Section~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;
Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 0: INITIALIZATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~ constant numbers: 
const label NC( mesh.nCells() );
const label NF(is3D ? mesh.nFaces() : (mesh.nFaces()-2*NC));
const label NFi = mesh.neighbour().size();

// ~~~~~~~~~~ constrant mesh fields:
scalarField sF_wf(mesh.surfaceInterpolation::weights().internalField());
scalarField sF_df(NF, 0.0);
vectorField vF_Sf(NF, vector(0.0, 0.0, 0.0));
scalarField sF_magSf(NF, 0.0);
vectorField vF_Cf(NF, vector(0.0, 0.0, 0.0));

// ~~~~~~~~~~ assigning Sf/Cf/df:
for(label faceI = 0; faceI < NFi; faceI++)
    {
    vF_Sf[faceI] = mesh.Sf()[faceI];
    vF_Cf[faceI] = mesh.Cf()[faceI];
    sF_df[faceI] = 1.0/mesh.surfaceInterpolation::deltaCoeffs().internalField()[faceI];
    }
forAll(mesh.boundaryMesh(), patchI)
    {
    word pName = mesh.boundary()[patchI].name();
    label patch_start =  mesh.boundary()[patchI].start();
    label patch_size  = mesh.boundary()[patchI].Cf().size();
    for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
        {
        if(faceI < NF)
            {
            vF_Sf[faceI] = mesh.boundary()[patchI].Sf()[faceI-patch_start];        
            vF_Cf[faceI] = mesh.boundary()[patchI].Cf()[faceI-patch_start];
            sF_df[faceI] = 1.0/mesh.surfaceInterpolation::deltaCoeffs().boundaryField()[patchI][faceI-patch_start];
            }
        } 
    }

// ~~~~~~~~~~ assigning magSf: 
sF_magSf = mag(vF_Sf);

// ~~~~~~~~~~ Some other mesh parameters
dvec Vec_vol(NC);
sp_dmat DMat_Vec_vol(NC, NC);
sp_dmat iDMat_Vec_vol(NC, NC);
for(label i=0; i<NC; ++i)
    {
    Vec_vol(i) = mesh.V()[i];
    DMat_Vec_vol(i, i) = Vec_vol(i);
    iDMat_Vec_vol(i, i) = 1.0/Vec_vol(i);
    }

dvec Vec_Sfx(NF);
dvec Vec_Sfy(NF);
dvec Vec_Sfz(is3D ? NF : 1);
sp_dmat DMat_Vec_Sfx(NF, NF);
sp_dmat DMat_Vec_Sfy(NF, NF);
sp_dmat DMat_Vec_Sfz((is3D ? NF : 1), (is3D ? NF : 1));
for(label i=0; i<NF; ++i)
    {
    Vec_Sfx(i) = mesh.Sf()[i].x();
    Vec_Sfy(i) = mesh.Sf()[i].y();
    if(is3D) Vec_Sfz(i) = mesh.Sf()[i].z(); 
    
    DMat_Vec_Sfx(i, i) = Vec_Sfx(i);
    DMat_Vec_Sfy(i, i) = Vec_Sfy(i);
    if(is3D) DMat_Vec_Sfz(i, i) = Vec_Sfz(i);
    }

// ~~~~~~~~~~ ownerList and neighList:
List<label> ownerList(NF, 0);
List<label> neighList(NFi, 0);

for(label faceI = 0; faceI < NF; faceI++)
    {
    if(faceI < NFi)
        {
        ownerList[faceI] = mesh.owner()[faceI];
        neighList[faceI] = mesh.neighbour()[faceI];
        }
    else
        {
        ownerList[faceI] = mesh.owner()[faceI];
        }
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Logical Indicators: 
// note: 0 is for fixed DOF (fixedValue), and 1 is for free boundary DOF (zeroGradient, ...)
uvec Log_p(NF, fill::ones);
uvec Log_nu(NF, fill::ones);
uvec Log_U(NF, fill::ones);
uvec Log_phi(NF, fill::ones);

forAll(mesh.boundaryMesh(), patchI)
    {
    word pName = mesh.boundary()[patchI].name();
    label patch_start =  mesh.boundary()[patchI].start();
    label patch_size  = mesh.boundary()[patchI].Cf().size();
    if(pName == "outlet")        
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            Log_p(faceI)  = 0;
            }
    if(pName == "inlet")        
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            Log_nu(faceI) = 0;       
            Log_phi(faceI)  = 0;
            Log_U(faceI)  = 0;
            }
    if((pName == "walls"))        
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {          
            Log_phi(faceI)  = 0;
            Log_U(faceI)  = 0;
            }
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ fixed and free DOF IDs:
uvec allFaceIDs(NF, fill::zeros);
for(label i=0; i<NF; ++i)
    {
    allFaceIDs(i) = i;
    }

uvec freeDOF_phi = allFaceIDs(find(Log_phi > 0.5));
const label NFF(freeDOF_phi.size()); // number of free phi dofs (phiF) contiguous.

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ State variables (read from solver)
// ~~~~~~~~~~ assigning p/U values
dvec Vec_p(NC);
dvec Vec_p2(NC);
dvec Vec_Ux(NC); 
dvec Vec_Uy(NC); 
dvec Vec_Uz(is3D ? NC : 1);
dvec Vec_gma(NC);
dvec Vec_phi(NF);
dvec Vec_nuf(NF);

for(label i=0; i<NC; ++i)
    {
    Vec_p(i)  = p.primitiveField()[i];
    Vec_p2(i) = p2.primitiveField()[i];
    Vec_Ux(i) = U.primitiveField()[i].x();
    Vec_Uy(i) = U.primitiveField()[i].y();
    if(is3D) Vec_Uz(i) = U.primitiveField()[i].z();
    Vec_gma(i) = gamma[i];
    }

for(label faceI = 0; faceI < NFi; faceI++)
    {
    Vec_phi(faceI) = phi.primitiveField()[faceI];
    }
forAll(mesh.boundaryMesh(), patchI)
    {
    label patch_start = mesh.boundary()[patchI].start();
    label patch_size  = mesh.boundary()[patchI].Cf().size();
    for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
        {
        if(faceI < NF)
            {        
            Vec_phi(faceI) = phi.boundaryField()[patchI][faceI-patch_start];
            }
        }        
    }

for(label i=0; i<NF; ++i)
    {
    Vec_nuf(i) = nu.value();
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 1a: R1 assembly
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: forming fvc::grad(p)
// ~~~~~~~~~~ linear interpolation (cell-to-face) matrix assembly:
sp_dmat Mat_cToF_p(NF, NC);

for(label faceI = 0; faceI < NF; faceI++) 
    {
    label cellI_o = ownerList[faceI];
    if(faceI < NFi)
        {
        label cellI_n = neighList[faceI];
        Mat_cToF_p(faceI, cellI_o) = sF_wf[faceI];
        Mat_cToF_p(faceI, cellI_n) = 1.0 - sF_wf[faceI];
        }
    else
        {
        if(Log_p(faceI) == 1)
            {     
            Mat_cToF_p(faceI, cellI_o) = 1.0;            
            }
        }
    }

// ~~~~~~~~~~ grad operator matrix assembly:
sp_dmat Mat_fToGx(NC, NF);
sp_dmat Mat_fToGy(NC, NF);
sp_dmat Mat_fToGz((is3D ? NC : 1), (is3D ? NF : 1));

for(label faceI = 0; faceI < NF; faceI++)
    {
    label cellI_o = ownerList[faceI];
    scalar vol_o  = Vec_vol(cellI_o);
    vector Sf_vec = vF_Sf[faceI];

    if (faceI < NFi)
        {           
        label cellI_n = neighList[faceI]; 
        scalar vol_n  = Vec_vol(cellI_n);
                                        
        Mat_fToGx(cellI_o, faceI) += Sf_vec.x()/vol_o;
        Mat_fToGy(cellI_o, faceI) += Sf_vec.y()/vol_o;
        if(is3D) Mat_fToGz(cellI_o, faceI) += Sf_vec.z()/vol_o;
                
        Mat_fToGx(cellI_n, faceI) -= Sf_vec.x()/vol_n;
        Mat_fToGy(cellI_n, faceI) -= Sf_vec.y()/vol_n;
        if(is3D) Mat_fToGz(cellI_n, faceI) -= Sf_vec.z()/vol_n;
        }
    else
        {                
        Mat_fToGx(cellI_o, faceI) += Sf_vec.x()/vol_o;
        Mat_fToGy(cellI_o, faceI) += Sf_vec.y()/vol_o;
        if(is3D) Mat_fToGz(cellI_o, faceI) += Sf_vec.z()/vol_o;
        }
    } 

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 2: forming fvm::laplacian(nu, U)
// ~~~~~~~~~~ Mat_Lap matrix and sources assembly:
sp_dmat Mat_Lap(NC, NC);
dvec    Vec_Lap_src_x(NC, fill::zeros);
dvec    Vec_Lap_src_y(NC, fill::zeros);
dvec    Vec_Lap_src_z(is3D ? NC : 1, fill::zeros);

for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI]; 
    
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_no   = sF_df[faceI];
    scalar nuf    = Vec_nuf(faceI);     
    
    Mat_Lap(cellI_o, cellI_o) +=  nuf*Sf_mag/d_no;
    Mat_Lap(cellI_o, cellI_n) -=  nuf*Sf_mag/d_no;   

    Mat_Lap(cellI_n, cellI_n) +=  nuf*Sf_mag/d_no;
    Mat_Lap(cellI_n, cellI_o) -=  nuf*Sf_mag/d_no;
    }

forAll(mesh.boundaryMesh(), patchI) 
    {
    word pName = mesh.boundary()[patchI].name();
    if((pName == "inlet") || (pName == "walls"))
        {
        label patch_start = mesh.boundary()[patchI].start();
        label patch_size  = mesh.boundary()[patchI].Cf().size();
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            label cellI_o = ownerList[faceI];
            scalar Sf_mag = mag(vF_Sf[faceI]);
            scalar nuf    = Vec_nuf(faceI);
            scalar d_fo   = sF_df[faceI];            
            vector Uf_vec     = U.boundaryField()[patchI][faceI-patch_start];
            vector nuFlux_vec = (nuf*Sf_mag/d_fo) * Uf_vec;
            
            Mat_Lap(cellI_o, cellI_o) += nuf*Sf_mag/d_fo;
            Vec_Lap_src_x(cellI_o) += nuFlux_vec.x();
            Vec_Lap_src_y(cellI_o) += nuFlux_vec.y();
            if(is3D) Vec_Lap_src_z(cellI_o) += nuFlux_vec.z();
            }        
        }
    }    

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3: forming fvm::div(phi, U)
// ~~~~~~~~~~ Mat_div matrix and sources assembly:
sp_dmat Mat_div(NC, NC);
dvec    Vec_div_src_x(NC, fill::zeros);
dvec    Vec_div_src_y(NC, fill::zeros);
dvec    Vec_div_src_z(is3D ? NC : 1, fill::zeros);

for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];
    
    scalar vol_flux = Vec_phi(faceI);

    if(vol_flux >= 0.0)
        {
        Mat_div(cellI_o, cellI_o) +=  vol_flux;        
        Mat_div(cellI_n, cellI_o) -=  vol_flux;
        }
    else
        {
        Mat_div(cellI_n, cellI_n) -=  vol_flux;
        Mat_div(cellI_o, cellI_n) +=  vol_flux;
        }
    }

forAll(mesh.boundaryMesh(), patchI)
    {
    word  pName       = mesh.boundary()[patchI].name();
    label patch_start = mesh.boundary()[patchI].start();
    label patch_size  = mesh.boundary()[patchI].Cf().size();
    if((pName == "inlet") || (pName == "walls"))
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            label cellI_o = ownerList[faceI];
            vector Uf_vec = U.boundaryField()[patchI][faceI-patch_start];

            Vec_div_src_x(cellI_o) -= Uf_vec.x() * Vec_phi(faceI);
            Vec_div_src_y(cellI_o) -= Uf_vec.y() * Vec_phi(faceI);
            if(is3D) Vec_div_src_z(cellI_o) -= Uf_vec.z() * Vec_phi(faceI);
            }      
    if((pName == "outlet")) // important: zeroGradient
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            label cellI_o   = ownerList[faceI];
            scalar vol_flux = Vec_phi(faceI);
            
            Mat_div(cellI_o, cellI_o) += vol_flux;
            }
    }

// ~~~~~~~~~~ Mat_1 matrix assembly:
dmat Mat_1(NFF, 4, fill::zeros);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    scalar vol_flux = Vec_phi(faceI);    
    if(vol_flux >= 0.0)
        {
        Mat_1(faceI, 0) = 1.0;
        Mat_1(faceI, 2) = -1.0;
        }
    else
        {
        Mat_1(faceI, 1) = 1.0; 
        Mat_1(faceI, 3) = -1.0;
        }
    }
for(label faceI = NFi; faceI < NFF; faceI++)
    {
    Mat_1(faceI, 0) = 1.0;
    }



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 4: forming BVP
sp_dmat DMat_BVP(NC, NC);
for(label i=0; i<NC; ++i)
    {
    DMat_BVP(i, i) = (gamma[i]/etta.value()) * Vec_vol(i);
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 5: Evaluating R1

// ~~~~~~~~~~ computing grad(p), using pf:
dvec Vec_Gp_x = Mat_fToGx * (Mat_cToF_p * Vec_p);
dvec Vec_Gp_y = Mat_fToGy * (Mat_cToF_p * Vec_p);
dvec Vec_Gp_z; if(is3D)
     Vec_Gp_z = Mat_fToGz * (Mat_cToF_p * Vec_p);

// ~~~~~~~~~~ computing R1_x,y,z:
dvec Vec_R1_x = (Mat_div + Mat_Lap + DMat_BVP) * Vec_Ux - (Vec_div_src_x + Vec_Lap_src_x) + (DMat_Vec_vol * Vec_Gp_x);
dvec Vec_R1_y = (Mat_div + Mat_Lap + DMat_BVP) * Vec_Uy - (Vec_div_src_y + Vec_Lap_src_y) + (DMat_Vec_vol * Vec_Gp_y);
dvec Vec_R1_z; if(is3D) 
     Vec_R1_z = (Mat_div + Mat_Lap + DMat_BVP) * Vec_Uz - (Vec_div_src_z + Vec_Lap_src_z) + (DMat_Vec_vol * Vec_Gp_z);
      
Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< endl;
cout  << "*** " << "Inf_norm(Vec_R1_x): " << norm(Vec_R1_x, "inf") << endl;
cout  << "*** " << "Inf_norm(Vec_R1_y): " << norm(Vec_R1_y, "inf") << endl;
cout  << "*** " << "Inf_norm(Vec_R1_z): " << norm(Vec_R1_z, "inf") << endl;
Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< endl;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 2: R2 assembly
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: forming rAU
// ~~~~~~~~~~ rAU:
dvec Vec_rAU(NC);

for(label i=0; i<NC; i++)
    {
    scalar Vc = Vec_vol(i);
    scalar Ac = Mat_div(i, i) + Mat_Lap(i, i) + DMat_BVP(i, i);

    Vec_rAU(i) = Vc/Ac;
    }

// ~~~~~~~~~~ diagonal rAU:
sp_dmat DMat_rAU(NC, NC);
for(label i=0; i<NC; ++i)
    {
    DMat_rAU(i, i) = Vec_rAU(i);
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 2: forming rAUf
// ~~~~~~~~~~ Mat_cToF_rAU matrix assembly:
sp_dmat  Mat_cToF_rAU(NF, NC);
for(label faceI = 0; faceI < NF; faceI++) 
    {
    label cellI_o   = ownerList[faceI];    
    if(faceI < NFi) // Internal faces
        {
        label cellI_n = neighList[faceI];

        Mat_cToF_rAU(faceI, cellI_o) = sF_wf[faceI];
        Mat_cToF_rAU(faceI, cellI_n) = 1.0 - sF_wf[faceI];        
        }
    else // All boundary faces
        {
        Mat_cToF_rAU(faceI, cellI_o) = 1.0;
        }
    }

// ~~~~~~~~~~ rAUf:
dvec Vec_rAUf = Mat_cToF_rAU * Vec_rAU;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3: forming HbyA
// ~~~~~~~~~~ H_bar:
sp_dmat Mat_H_bar = Mat_div + Mat_Lap + DMat_BVP;
for(label i=0; i<NC; ++i)
    {
    Mat_H_bar(i, i) = 0.0;
    }

// ~~~~~~~~~~ H:
dvec Vec_Hx = ((Vec_div_src_x + Vec_Lap_src_x) - Mat_H_bar * Vec_Ux)/ Vec_vol;
dvec Vec_Hy = ((Vec_div_src_y + Vec_Lap_src_y) - Mat_H_bar * Vec_Uy)/ Vec_vol;
dvec Vec_Hz; if(is3D)
     Vec_Hz = ((Vec_div_src_z + Vec_Lap_src_z) - Mat_H_bar * Vec_Uz)/ Vec_vol;


// ~~~~~~~~~~ diag(H):
sp_dmat DMat_Hx(NC, NC);
sp_dmat DMat_Hy(NC, NC);
sp_dmat DMat_Hz((is3D ? NC : 1), (is3D ? NC : 1));
for(label i=0; i<NC; ++i)
    {
    DMat_Hx(i, i) = Vec_Hx(i);
    DMat_Hy(i, i) = Vec_Hy(i); 
    if(is3D) DMat_Hz(i, i) = Vec_Hz(i);
    }

// ~~~~~~~~~~ HbyA:
dvec Vec_HbyA_x = (Vec_Hx % Vec_rAU);
dvec Vec_HbyA_y = (Vec_Hy % Vec_rAU);
dvec Vec_HbyA_z; if(is3D)
     Vec_HbyA_z = (Vec_Hz % Vec_rAU);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 4: forming HbyAf
// ~~~~~~~~~~ Mat_cToF_HbyA matrix assembly:
sp_dmat Mat_cToF_HbyA(NF, NC);
for(label faceI = 0; faceI < NF; faceI++) 
    {
    label cellI_o = ownerList[faceI];
    if(faceI < NFi)
        {
        label cellI_n = neighList[faceI];

        Mat_cToF_HbyA(faceI, cellI_o) = sF_wf[faceI];
        Mat_cToF_HbyA(faceI, cellI_n) = 1.0 - sF_wf[faceI];        
        }
    else
        {
        if(Log_U(faceI) == 1)
            {     
            Mat_cToF_HbyA(faceI, cellI_o) = 1.0;
            }   
        }
    }

// ~~~~~~~~~~ HbyAf_b:
dvec Vec_HbyAfb_x(NF, fill::zeros);
dvec Vec_HbyAfb_y(NF, fill::zeros);
dvec Vec_HbyAfb_z(is3D ? NF : 1, fill::zeros);
forAll(mesh.boundaryMesh(), patchI)
    {
    word pName = mesh.boundary()[patchI].name();
    label patch_start =  mesh.boundary()[patchI].start();
    label patch_size  = mesh.boundary()[patchI].Cf().size();
    if((pName == "inlet") || (pName == "walls"))
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            vector Uf_vec = U.boundaryField()[patchI][faceI-patch_start];
            Vec_HbyAfb_x(faceI) = Uf_vec.x();
            Vec_HbyAfb_y(faceI) = Uf_vec.y(); 
            if(is3D) Vec_HbyAfb_z(faceI) = Uf_vec.z();   
            }
    }

// ~~~~~~~~~~ HbyAf:
dvec Vec_HbyAf_x = (Mat_cToF_HbyA * Vec_HbyA_x) + Vec_HbyAfb_x;
dvec Vec_HbyAf_y = (Mat_cToF_HbyA * Vec_HbyA_y) + Vec_HbyAfb_y;
dvec Vec_HbyAf_z; if(is3D) 
     Vec_HbyAf_z = (Mat_cToF_HbyA * Vec_HbyA_z) + Vec_HbyAfb_z;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 5: forming phiHbyA: fvc::flux(HbyA)
// ~~~~~~~~~~ phiHbyA
dvec Vec_phiHbyA;
if(is3D)
    Vec_phiHbyA = DMat_Vec_Sfx * Vec_HbyAf_x + DMat_Vec_Sfy * Vec_HbyAf_y + DMat_Vec_Sfz * Vec_HbyAf_z;
else
    Vec_phiHbyA = DMat_Vec_Sfx * Vec_HbyAf_x + DMat_Vec_Sfy * Vec_HbyAf_y;

// ~~~~~~~~~~ Mat_div_phiHbyA  matrix assembly:
sp_dmat Mat_div_phiHbyA(NC, NF);
for(label faceI = 0; faceI < NF; faceI++)
    {
    label cellI_o  = ownerList[faceI];
    vector Sf_vec  = vF_Sf[faceI];    
    vector r_o     = mesh.C()[cellI_o];
    vector r_f     = vF_Cf[faceI];
    vector r_of(r_f - r_o);    
    scalar SGN((r_of & Sf_vec)/(mag(Sf_vec)*mag(r_of)));
    
    if(faceI < NFi)
        {
        label cellI_n = neighList[faceI];
        Mat_div_phiHbyA(cellI_o, faceI) += SGN / Vec_vol(cellI_o);
        Mat_div_phiHbyA(cellI_n, faceI) -= SGN / Vec_vol(cellI_n);
        }
    else
        {
        Mat_div_phiHbyA(cellI_o, faceI) += SGN / Vec_vol(cellI_o);
        }
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 6: forming fvc::div(phiHbyA)
//dvec Vec_div_phiHbyA = Mat_div_phiHbyA * Vec_phiHbyA;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 7: forming fvm::laplacian(rAU, p)
// ~~~~~~~~~~ Mat_Lap_p matrix and sources assembly:
sp_dmat Mat_Lap_p(NC, NC);

for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];                   
    
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_no = sF_df[faceI];
    scalar rAUf = Vec_rAUf(faceI);        
    
    Mat_Lap_p(cellI_o, cellI_o) +=  rAUf*Sf_mag/d_no;
    Mat_Lap_p(cellI_o, cellI_n) -=  rAUf*Sf_mag/d_no;

    Mat_Lap_p(cellI_n, cellI_n) +=  rAUf*Sf_mag/d_no;
    Mat_Lap_p(cellI_n, cellI_o) -=  rAUf*Sf_mag/d_no;
    } 
forAll(mesh.boundaryMesh(), patchI)
    {
    word pName = mesh.boundary()[patchI].name();
    if((pName == "outlet"))
        {
        label patch_start = mesh.boundary()[patchI].start();
        label patch_size  = mesh.boundary()[patchI].Cf().size();
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            label cellI_o = ownerList[faceI];
            scalar Sf_mag = sF_magSf[faceI];
            scalar d_fo   = sF_df[faceI]; 
            scalar rAUf   = Vec_rAUf(faceI); 

            Mat_Lap_p(cellI_o, cellI_o) +=  rAUf*Sf_mag/d_fo;
            }        
        }   
    }

// ~~~~~~~~~~ R2:
dvec Vec_R2 = -(Mat_Lap_p * Vec_p2) - (DMat_Vec_vol * (Mat_div_phiHbyA * Vec_phiHbyA));

Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< endl;
cout  << "*** " << "Inf_norm(Vec_R2): " << norm(Vec_R2, "inf") << endl;
Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< endl;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Part 3, R3 assembly
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: forming pEqn.flux()
// ~~~~~~~~~~ M_pFlux matrix assembly:
sp_dmat Mat_pFlux(NF, NC);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI]; 
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_no   = sF_df[faceI];
    scalar rAUf   = Vec_rAUf(faceI); 

    Mat_pFlux(faceI, cellI_o) = -(rAUf*Sf_mag/d_no);
    Mat_pFlux(faceI, cellI_n) =  (rAUf*Sf_mag/d_no);
    }
forAll(mesh.boundaryMesh(), patchI)
    {
    word pName = mesh.boundary()[patchI].name();
    if((pName == "outlet"))
        {
        label patch_start = mesh.boundary()[patchI].start();
        label patch_size  = mesh.boundary()[patchI].Cf().size();
        for(label faceI = patch_start; faceI < patch_start+patch_size; faceI++)
            {
            label cellI_o = ownerList[faceI];            
            scalar Sf_mag = sF_magSf[faceI];       
            scalar rAUf   = Vec_rAUf(faceI);
            scalar d_fo   = sF_df[faceI];
            
            Mat_pFlux(faceI, cellI_o) = -(rAUf*Sf_mag/d_fo);
            }        
        }   
    }

// ~~~~~~~~~~ R3:
dvec Vec_R3 = Vec_phi - (Vec_phiHbyA - Mat_pFlux * Vec_p2);

Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
cout  << "*** " << "Inf_norm(Vec_R3): " << norm(Vec_R3, "inf") << endl;
Info  << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;

