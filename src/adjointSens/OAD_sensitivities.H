// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 1: R1 partial derivatives
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: dR1 WRT U
// ~~~~~~~~~~ Note: dR1_dU = dR1x_dUx = dR1y_dUy = dR1z_dUz
sp_dmat Mat_dR1_dU = Mat_div + Mat_Lap + DMat_BVP;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 2: dR1 WRT p
sp_dmat Mat_dR1x_dp = DMat_Vec_vol * (Mat_fToGx * Mat_cToF_p);
sp_dmat Mat_dR1y_dp = DMat_Vec_vol * (Mat_fToGy * Mat_cToF_p);
sp_dmat Mat_dR1z_dp; if(is3D)
        Mat_dR1z_dp = DMat_Vec_vol * (Mat_fToGz * Mat_cToF_p);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3: dR1 WRT gamma
sp_dmat DMat_dR1x_dgma(NC, NC);
sp_dmat DMat_dR1y_dgma(NC, NC);
sp_dmat DMat_dR1z_dgma(is3D ? NC : 1, is3D ? NC : 1);
for(label i=0; i<NC; ++i)
    {
    DMat_dR1x_dgma(i, i) = (Vec_vol(i) / etta.value()) * Vec_Ux(i);
    DMat_dR1y_dgma(i, i) = (Vec_vol(i) / etta.value()) * Vec_Uy(i);    
    if(is3D) DMat_dR1z_dgma(i, i) = (Vec_vol(i) / etta.value()) * Vec_Uz(i);
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 4: dR1 WRT phi
sp_dmat Mat_dR1x_dphiF(NC, NFF);
sp_dmat Mat_dR1y_dphiF(NC, NFF);
sp_dmat Mat_dR1z_dphiF(is3D ? NC : 1, is3D ? NFF : 1);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];

    dmat G = Mat_1.row(faceI);
    
    Mat_dR1x_dphiF(cellI_o, faceI) = G(0) * Vec_Ux(cellI_o) + G(1) * Vec_Ux(cellI_n);
    Mat_dR1x_dphiF(cellI_n, faceI) = G(2) * Vec_Ux(cellI_o) + G(3) * Vec_Ux(cellI_n);

    Mat_dR1y_dphiF(cellI_o, faceI) = G(0) * Vec_Uy(cellI_o) + G(1) * Vec_Uy(cellI_n);
    Mat_dR1y_dphiF(cellI_n, faceI) = G(2) * Vec_Uy(cellI_o) + G(3) * Vec_Uy(cellI_n);

    if(is3D) Mat_dR1z_dphiF(cellI_o, faceI) = G(0) * Vec_Uz(cellI_o) + G(1) * Vec_Uz(cellI_n);
    if(is3D) Mat_dR1z_dphiF(cellI_n, faceI) = G(2) * Vec_Uz(cellI_o) + G(3) * Vec_Uz(cellI_n);
    }
for(label faceI = NFi; faceI < NFF; faceI++)
    {
    label cellI_o = ownerList[faceI];

    Mat_dR1x_dphiF(cellI_o, faceI) = Vec_Ux(cellI_o);
    Mat_dR1y_dphiF(cellI_o, faceI) = Vec_Uy(cellI_o);
    if(is3D) Mat_dR1z_dphiF(cellI_o, faceI) = Vec_Uz(cellI_o);
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 2: R2 partial derivatives
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~ constant factors:
sp_dmat Mat_Gx = DMat_Vec_vol * Mat_div_phiHbyA * DMat_Vec_Sfx * Mat_cToF_HbyA;
sp_dmat Mat_Gy = DMat_Vec_vol * Mat_div_phiHbyA * DMat_Vec_Sfy * Mat_cToF_HbyA;
sp_dmat Mat_Gz; if(is3D)
        Mat_Gz = DMat_Vec_vol * Mat_div_phiHbyA * DMat_Vec_Sfz * Mat_cToF_HbyA;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: dR2 WRT p
sp_dmat Mat_dR2_dp = -Mat_Lap_p;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 2: dR2 WRT Ux, Uy, Uz
sp_dmat Mat_dR2_dUx = Mat_Gx * DMat_rAU * iDMat_Vec_vol * Mat_H_bar;
sp_dmat Mat_dR2_dUy = Mat_Gy * DMat_rAU * iDMat_Vec_vol * Mat_H_bar;
sp_dmat Mat_dR2_dUz; if(is3D)
        Mat_dR2_dUz = Mat_Gz * DMat_rAU * iDMat_Vec_vol * Mat_H_bar;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3: dR2 WRT gamma
// ~~~~~~~~~~ Note: R2_I  = -(Mat_Lap_p * Vec_p2)
// ~~~~~~~~~~ Note: R2_II = -(DMat_Vec_vol * (Mat_div_phiHbyA * Vec_phiHbyA))

sp_dmat Mat_dR2_I_drAUfF(NC, NFF);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];    
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_fo   = sF_df[faceI];
    
    Mat_dR2_I_drAUfF(cellI_o, faceI) = (Sf_mag/d_fo) * (Vec_p2(cellI_o) - Vec_p2(cellI_n));
    Mat_dR2_I_drAUfF(cellI_n, faceI) = (Sf_mag/d_fo) * (Vec_p2(cellI_n) - Vec_p2(cellI_o));
    }
for(label faceI = NFi; faceI < NFF; faceI++)
    {
    label cellI_o = ownerList[faceI];
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_fo   = sF_df[faceI];        
    Mat_dR2_I_drAUfF(cellI_o, faceI) = (Sf_mag/d_fo) * Vec_p2(cellI_o);
    }

// ~~~~~~~~~~ drAUf_drAU:
sp_dmat Mat_drAUf_drAU = Mat_cToF_rAU.rows(0, NFF-1);

sp_dmat Mat_dR2_I_drAU = Mat_dR2_I_drAUfF * Mat_drAUf_drAU;
sp_dmat Mat_drAU_dgma = (-1.0/etta.value()) * DMat_rAU * DMat_rAU;

// ~~~~~~~~~~ dR2_I WRT gamma
sp_dmat Mat_dR2_I_dgma = -1.0 * Mat_dR2_I_drAU * Mat_drAU_dgma;

// ~~~~~~~~~~ dR2_II WRT gamma
sp_dmat Mat_dR2_II_drAU;
if(is3D)
    Mat_dR2_II_drAU = -(Mat_Gx * DMat_Hx + Mat_Gy * DMat_Hy + Mat_Gz * DMat_Hz);
else
    Mat_dR2_II_drAU = -(Mat_Gx * DMat_Hx + Mat_Gy * DMat_Hy);

sp_dmat Mat_dR2_II_dgma = Mat_dR2_II_drAU * Mat_drAU_dgma;

// ~~~~~~~~~~ dR2 WRT gamma
sp_dmat Mat_dR2_dgma = Mat_dR2_I_dgma + Mat_dR2_II_dgma;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3: dR2 WRT phi
sp_dmat Mat_drAU_dphiF(NC, NFF);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];
    
    double G_oo = Mat_1(faceI, 0);
    double G_nn = Mat_1(faceI, 3);
    
    Mat_drAU_dphiF(cellI_o, faceI) = -G_oo/Vec_vol(cellI_o) * pow(Vec_rAU(cellI_o), 2);
    Mat_drAU_dphiF(cellI_n, faceI) = -G_nn/Vec_vol(cellI_n) * pow(Vec_rAU(cellI_n), 2);
    }
for(label faceI = NFi; faceI < NFF; faceI++)
    {
    label cellI_o = ownerList[faceI];
    double G_oo = Mat_1(faceI, 0);

    Mat_drAU_dphiF(cellI_o, faceI) = -G_oo/Vec_vol(cellI_o) * pow(Vec_rAU(cellI_o), 2);
    }

// ~~~~~~~~~~ dR2_I WRT phiF
sp_dmat Mat_dR2_I_dphiF = -1.0 * Mat_dR2_I_drAU * Mat_drAU_dphiF;

// ~~~~~~~~~~ dH WRT phi
sp_dmat Mat_dHx_dphiF(NC, NFF);
sp_dmat Mat_dHy_dphiF(NC, NFF);
sp_dmat Mat_dHz_dphiF(is3D ? NC : 1, is3D ? NFF : 1);
for(label faceI = 0; faceI < NFi; faceI++)
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI];
    
    dmat G = Mat_1.row(faceI);
    
    Mat_dHx_dphiF(cellI_o, faceI) = -1.0 * G(1) * Vec_Ux(cellI_n) / Vec_vol(cellI_o);
    Mat_dHx_dphiF(cellI_n, faceI) = -1.0 * G(2) * Vec_Ux(cellI_o) / Vec_vol(cellI_n);

    Mat_dHy_dphiF(cellI_o, faceI) = -1.0 * G(1) * Vec_Uy(cellI_n) / Vec_vol(cellI_o);
    Mat_dHy_dphiF(cellI_n, faceI) = -1.0 * G(2) * Vec_Uy(cellI_o) / Vec_vol(cellI_n);

    if(is3D) Mat_dHz_dphiF(cellI_o, faceI) = -1.0 * G(1) * Vec_Uz(cellI_n)/ Vec_vol(cellI_o); 
    if(is3D) Mat_dHz_dphiF(cellI_n, faceI) = -1.0 * G(2) * Vec_Uz(cellI_o)/ Vec_vol(cellI_n);
    }

sp_dmat Mat_dR2_II_dphiF;
if(is3D)
    {
    Mat_dR2_II_dphiF = - Mat_Gx * (DMat_Hx * Mat_drAU_dphiF + DMat_rAU * Mat_dHx_dphiF) 
                       - Mat_Gy * (DMat_Hy * Mat_drAU_dphiF + DMat_rAU * Mat_dHy_dphiF)
                       - Mat_Gz * (DMat_Hz * Mat_drAU_dphiF + DMat_rAU * Mat_dHz_dphiF);
    }
else
    {
    Mat_dR2_II_dphiF = - Mat_Gx * (DMat_Hx * Mat_drAU_dphiF + DMat_rAU * Mat_dHx_dphiF) 
                       - Mat_Gy * (DMat_Hy * Mat_drAU_dphiF + DMat_rAU * Mat_dHy_dphiF);
    }


// ~~~~~~~~~~ dR2 WRT phiF
sp_dmat Mat_dR2_dphiF = Mat_dR2_I_dphiF + Mat_dR2_II_dphiF;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 3c: R3F partial derivatives
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 1: dR3 WRT p
sp_dmat Mat_dR3_dp = Mat_pFlux;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 2: dR3 WRT Ux, Uy, Uz
sp_dmat Mat_dR3_dUx = DMat_Vec_Sfx * Mat_cToF_HbyA * iDMat_Vec_vol * DMat_rAU * Mat_H_bar;
sp_dmat Mat_dR3_dUy = DMat_Vec_Sfy * Mat_cToF_HbyA * iDMat_Vec_vol * DMat_rAU * Mat_H_bar;
sp_dmat Mat_dR3_dUz; if(is3D)
        Mat_dR3_dUz = DMat_Vec_Sfz * Mat_cToF_HbyA * iDMat_Vec_vol * DMat_rAU * Mat_H_bar;
        
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STEP 3-4: dR3 WRT gamma/phi
// ~~~~~~~~~~ computing dR3_III_dphi:

sp_dmat Mat_dR3_III_drAUfF(NF, NFF);
for(label faceI = 0; faceI < NFi; faceI++) 
    {
    label cellI_o = ownerList[faceI];
    label cellI_n = neighList[faceI]; 
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_no   = sF_df[faceI];

    Mat_dR3_III_drAUfF(faceI, faceI) = (Sf_mag/d_no) * (-Vec_p2(cellI_o) + Vec_p2(cellI_n));
    }
for(label faceI = NFi; faceI < NFF; faceI++)
    {
    label cellI_o = ownerList[faceI];
    scalar Sf_mag = sF_magSf[faceI];
    scalar d_no   = sF_df[faceI];
    
    Mat_dR3_III_drAUfF(faceI, faceI) = -1.0 * (Sf_mag/d_no) * Vec_p2(cellI_o);
    }

sp_dmat Mat_dR3_III_drAU  = Mat_dR3_III_drAUfF * Mat_drAUf_drAU;
sp_dmat Mat_dR3_III_dphiF = Mat_dR3_III_drAU  * Mat_drAU_dphiF;

// ~~~~~~~~~~ dR3_III_dgma
sp_dmat Mat_dR3_III_dgma = Mat_dR3_III_drAU * Mat_drAU_dgma;

// ~~~~~~~~~~ dR3_II_dphiF
sp_dmat Mat_dR3_II_drAU;
if(is3D)
    {
    Mat_dR3_II_drAU = -DMat_Vec_Sfx * Mat_cToF_HbyA * DMat_Hx
                      -DMat_Vec_Sfy * Mat_cToF_HbyA * DMat_Hy
                      -DMat_Vec_Sfz * Mat_cToF_HbyA * DMat_Hz;    
    }
else
    {
    Mat_dR3_II_drAU = -DMat_Vec_Sfx * Mat_cToF_HbyA * DMat_Hx
                      -DMat_Vec_Sfy * Mat_cToF_HbyA * DMat_Hy;
    }

sp_dmat Mat_dR3_II_dphiF;
if(is3D)
    {
    Mat_dR3_II_dphiF = -DMat_Vec_Sfx * Mat_cToF_HbyA * (DMat_Hx * Mat_drAU_dphiF + DMat_rAU * Mat_dHx_dphiF)
                       -DMat_Vec_Sfy * Mat_cToF_HbyA * (DMat_Hy * Mat_drAU_dphiF + DMat_rAU * Mat_dHy_dphiF)
                       -DMat_Vec_Sfz * Mat_cToF_HbyA * (DMat_Hz * Mat_drAU_dphiF + DMat_rAU * Mat_dHz_dphiF);
    }
else
    {
    Mat_dR3_II_dphiF = -DMat_Vec_Sfx * Mat_cToF_HbyA * (DMat_Hx * Mat_drAU_dphiF + DMat_rAU * Mat_dHx_dphiF)
                       -DMat_Vec_Sfy * Mat_cToF_HbyA * (DMat_Hy * Mat_drAU_dphiF + DMat_rAU * Mat_dHy_dphiF);
    }

// ~~~~~~~~~~ dR3_II_dgma
sp_dmat Mat_dR3_II_dgma = Mat_dR3_II_drAU * Mat_drAU_dgma;

// ~~~~~~~~~~ dR3_I_dphiF
sp_dmat DMat_dR3_I_dphi   = speye<sp_dmat>(NF, NF);
sp_dmat DMat_dR3_I_dphiF  = DMat_dR3_I_dphi.cols(0, NFF-1);

// ~~~~~~~~~~ dR3_dphiF
sp_dmat Mat_dR3_dphi = DMat_dR3_I_dphiF + Mat_dR3_II_dphiF + Mat_dR3_III_dphiF;

// ~~~~~~~~~~ dR3_dgma
sp_dmat Mat_dR3_dgma = Mat_dR3_II_dgma + Mat_dR3_III_dgma;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dR3F_d():
sp_dmat Mat_dR3F_dUx   = Mat_dR3_dUx.rows(0, NFF-1);
sp_dmat Mat_dR3F_dUy   = Mat_dR3_dUy.rows(0, NFF-1);
sp_dmat Mat_dR3F_dUz; if(is3D)
        Mat_dR3F_dUz   = Mat_dR3_dUz.rows(0, NFF-1);

sp_dmat Mat_dR3F_dp    = Mat_dR3_dp.rows(0, NFF-1);
sp_dmat Mat_dR3F_dphi  = Mat_dR3_dphi.rows(0, NFF-1);
sp_dmat Mat_dR3F_dgma  = Mat_dR3_dgma.rows(0, NFF-1);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Residual Matrix:
sp_dmat Mat_R(is3D ? 4*NC+NFF : 3*NC+NFF, is3D ? 4*NC+NFF : 3*NC+NFF);

const label lRow = is3D ? 4*NC+NFF-1 : 3*NC+NFF-1;
const label lCol = lRow;

if(is3D)
    {
    Mat_R.submat(0     , 0     , 1*NC-1 , 1*NC-1) = trans(Mat_dR1_dU); 
    Mat_R.submat(3*NC  , 0     , 4*NC-1 , 1*NC-1) = trans(Mat_dR1x_dp);
    Mat_R.submat(4*NC  , 0     , lRow   , 1*NC-1) = trans(Mat_dR1x_dphiF);

    Mat_R.submat(1*NC  , 1*NC  , 2*NC-1 , 2*NC-1) = trans(Mat_dR1_dU);
    Mat_R.submat(3*NC  , 1*NC  , 4*NC-1 , 2*NC-1) = trans(Mat_dR1y_dp);
    Mat_R.submat(4*NC  , 1*NC  , lRow   , 2*NC-1) = trans(Mat_dR1y_dphiF);

    Mat_R.submat(2*NC  , 2*NC  , 3*NC-1 , 3*NC-1) = trans(Mat_dR1_dU);
    Mat_R.submat(3*NC  , 2*NC  , 4*NC-1 , 3*NC-1) = trans(Mat_dR1z_dp);
    Mat_R.submat(4*NC  , 2*NC  , lRow   , 3*NC-1) = trans(Mat_dR1z_dphiF);

    Mat_R.submat(0     , 3*NC  , 1*NC-1 , 4*NC-1) = trans(Mat_dR2_dUx);
    Mat_R.submat(1*NC  , 3*NC  , 2*NC-1 , 4*NC-1) = trans(Mat_dR2_dUy);
    Mat_R.submat(2*NC  , 3*NC  , 3*NC-1 , 4*NC-1) = trans(Mat_dR2_dUz);
    Mat_R.submat(3*NC  , 3*NC  , 4*NC-1 , 4*NC-1) = trans(Mat_dR2_dp);
    Mat_R.submat(4*NC  , 3*NC  , lRow   , 4*NC-1) = trans(Mat_dR2_dphiF);

    Mat_R.submat(0     , 4*NC  , 1*NC-1 , lCol  ) = trans(Mat_dR3F_dUx);
    Mat_R.submat(1*NC  , 4*NC  , 2*NC-1 , lCol  ) = trans(Mat_dR3F_dUy);
    Mat_R.submat(2*NC  , 4*NC  , 3*NC-1 , lCol  ) = trans(Mat_dR3F_dUz);
    Mat_R.submat(3*NC  , 4*NC  , 4*NC-1 , lCol  ) = trans(Mat_dR3F_dp);
    Mat_R.submat(4*NC  , 4*NC  , lRow   , lCol  ) = trans(Mat_dR3F_dphi);
    }
else
    {
    Mat_R.submat(0     , 0     , 1*NC-1 , 1*NC-1) = trans(Mat_dR1_dU); 
    Mat_R.submat(2*NC  , 0     , 3*NC-1 , 1*NC-1) = trans(Mat_dR1x_dp);
    Mat_R.submat(3*NC  , 0     , lRow   , 1*NC-1) = trans(Mat_dR1x_dphiF);

    Mat_R.submat(1*NC  , 1*NC  , 2*NC-1 , 2*NC-1) = trans(Mat_dR1_dU);
    Mat_R.submat(2*NC  , 1*NC  , 3*NC-1 , 2*NC-1) = trans(Mat_dR1y_dp);
    Mat_R.submat(3*NC  , 1*NC  , lRow   , 2*NC-1) = trans(Mat_dR1y_dphiF);

    Mat_R.submat(0     , 2*NC  , 1*NC-1 , 3*NC-1) = trans(Mat_dR2_dUx);
    Mat_R.submat(1*NC  , 2*NC  , 2*NC-1 , 3*NC-1) = trans(Mat_dR2_dUy);
    Mat_R.submat(2*NC  , 2*NC  , 3*NC-1 , 3*NC-1) = trans(Mat_dR2_dp);
    Mat_R.submat(3*NC  , 2*NC  , lRow   , 3*NC-1) = trans(Mat_dR2_dphiF);

    Mat_R.submat(0     , 3*NC  , 1*NC-1 , lCol  ) = trans(Mat_dR3F_dUx);
    Mat_R.submat(1*NC  , 3*NC  , 2*NC-1 , lCol  ) = trans(Mat_dR3F_dUy);
    Mat_R.submat(2*NC  , 3*NC  , 3*NC-1 , lCol  ) = trans(Mat_dR3F_dp);
    Mat_R.submat(3*NC  , 3*NC  , lRow   , lCol  ) = trans(Mat_dR3F_dphi);    
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PART 4: Sensetivities
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ forces:
scalar Fx = sum(mesh.V() * (U.component(vector::X) * gamma)) / etta.value();
scalar Fy = sum(mesh.V() * (U.component(vector::Y) * gamma)) / etta.value();
scalar Fz = sum(mesh.V() * (U.component(vector::Z) * gamma)) / etta.value();
Info << "Fx = " << Fx << nl << "Fy = " << Fy << nl << "Fz = " << Fz << endl; 

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ partial derivatives:
dvec Vec_dFx_dUx(NC, fill::zeros);
dvec Vec_dFy_dUy(NC, fill::zeros); 
dvec Vec_dFx_dgma(NC, fill::zeros); 
dvec Vec_dFy_dgma(NC, fill::zeros);
for(label i=0; i<NC; i++)
    {
    Vec_dFx_dUx(i)  = gamma[i]  * (Vec_vol(i) / etta.value());
    Vec_dFy_dUy(i)  = Vec_dFx_dUx(i);
    
    Vec_dFx_dgma(i) = Vec_Ux(i) * (Vec_vol(i) / etta.value());
    Vec_dFy_dgma(i) = Vec_Uy(i) * (Vec_vol(i) / etta.value());
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Obj/Const func values:
// Note: J = Fy; C = Fx;

// ~~~~~~~~~~~~~~~~~~~~ Fx:
dvec Vec_dC_dgma = Vec_dFx_dgma;
dvec Vec_dC_dUx  = Vec_dFx_dUx;
dvec Vec_dC_dUy(NC, fill::zeros);

// ~~~~~~~~~~~~~~~~~~~~ Fy:
dvec Vec_dJ_dgma = Vec_dFy_dgma;
dvec Vec_dJ_dUx(NC, fill::zeros);
dvec Vec_dJ_dUy  = Vec_dFy_dUy;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Adjoint r.h.s. vectors:
dvec Vec_R_src_J(is3D ? 4*NC+NFF : 3*NC+NFF, fill::zeros);
dvec Vec_R_src_C(is3D ? 4*NC+NFF : 3*NC+NFF, fill::zeros);

// apply objective partials on r.h.s
Vec_R_src_J.subvec(0, NC-1)    = -1.0 * Vec_dJ_dUx;
Vec_R_src_J.subvec(NC, 2*NC-1) = -1.0 * Vec_dJ_dUy;

Vec_R_src_C.subvec(0, NC-1)    = -1.0 * Vec_dC_dUx;
Vec_R_src_C.subvec(NC, 2*NC-1) = -1.0 * Vec_dC_dUy;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Matlab as solver:
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save data for Matlab solver:
sp_mat::const_iterator it     = Mat_R.begin();
sp_mat::const_iterator it_end = Mat_R.end();
ofstream outfile;

outfile.open("data/I.bin", ios::binary | ios::out);
for(it = Mat_R.begin(); it != it_end; ++it)
    {
    int I = it.row();
    outfile.write((char*) &I, sizeof(int));
    }
outfile.close();

outfile.open("data/J.bin", ios::binary | ios::out);
for(it = Mat_R.begin(); it != it_end; ++it)
    {
    int J = it.col();
    outfile.write((char*) &J, sizeof(int));
    }
outfile.close();

outfile.open("data/VAL.bin", ios::binary | ios::out);
for(it = Mat_R.begin(); it != it_end; ++it)
    {
    double val = (*it);
    outfile.write((char*) &val, sizeof(double));
    }
outfile.close();

outfile.open("data/Vec_R_src_J.bin", ios::binary | ios::out);
for(unsigned int i=0; i< Vec_R_src_J.size(); ++i)
    {
    double val = Vec_R_src_J(i);   
    outfile.write((char*) &val, sizeof(double));
    }
outfile.close();

outfile.open("data/Vec_R_src_C.bin", ios::binary | ios::out);
for(unsigned int i=0; i< Vec_R_src_C.size(); ++i)
    {
    double val = Vec_R_src_C(i);   
    outfile.write((char*) &val, sizeof(double));
    }
outfile.close();

// ~~~~~~~~~~ call Solve Adjoint in Matlab
system ("/path/to/matlab -nodisplay -nosplash -nodesktop -r \"run(\'adjointSolver.m\'); exit;\" ");

// ~~~~~~~~~~ read the solution:
dvec Lambda_J, Lambda_C;
Lambda_J.load("data/Lambda_J.bin");
Lambda_C.load("data/Lambda_C.bin");

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ decompose adjoint vector:
dvec Lambda_J_1x(NC), Lambda_J_1y(NC), Lambda_J_1z(is3D ? NC : 1), Lambda_J_2(NC), Lambda_J_3F(NFF);
dvec Lambda_C_1x(NC), Lambda_C_1y(NC), Lambda_C_1z(is3D ? NC : 1), Lambda_C_2(NC), Lambda_C_3F(NFF);

if(is3D)
    {
    Lambda_J_1x = Lambda_J.subvec(0   , 1*NC-1);
    Lambda_J_1y = Lambda_J.subvec(1*NC, 2*NC-1);
    Lambda_J_1z = Lambda_J.subvec(2*NC, 3*NC-1);
    Lambda_J_2  = Lambda_J.subvec(3*NC, 4*NC-1);
    Lambda_J_3F = Lambda_J.subvec(4*NC, lRow);
    
    Lambda_C_1x = Lambda_C.subvec(0   , 1*NC-1);
    Lambda_C_1y = Lambda_C.subvec(1*NC, 2*NC-1);
    Lambda_C_1z = Lambda_C.subvec(2*NC, 3*NC-1);
    Lambda_C_2  = Lambda_C.subvec(3*NC, 4*NC-1);
    Lambda_C_3F = Lambda_C.subvec(4*NC, lRow);    
    }
else
    {
    Lambda_J_1x = Lambda_J.subvec(0   , 1*NC-1);
    Lambda_J_1y = Lambda_J.subvec(1*NC, 2*NC-1);
    Lambda_J_2  = Lambda_J.subvec(2*NC, 3*NC-1);
    Lambda_J_3F = Lambda_J.subvec(3*NC, lRow);

    Lambda_C_1x = Lambda_C.subvec(0   , 1*NC-1);
    Lambda_C_1y = Lambda_C.subvec(1*NC, 2*NC-1);
    Lambda_C_2  = Lambda_C.subvec(2*NC, 3*NC-1);
    Lambda_C_3F = Lambda_C.subvec(3*NC, lRow);    
    }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Total Lagrangian Derivetives:
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmat DLJ_Dgma, DLC_Dgma;

if(is3D)
    {
    DLJ_Dgma = trans(Vec_dJ_dgma) + trans(Lambda_J_1x) * DMat_dR1x_dgma
                                  + trans(Lambda_J_1y) * DMat_dR1y_dgma
                                  + trans(Lambda_J_1z) * DMat_dR1z_dgma
                                  + trans(Lambda_J_2)  * Mat_dR2_dgma
                                  + trans(Lambda_J_3F) * Mat_dR3F_dgma;

    DLC_Dgma = trans(Vec_dC_dgma) + trans(Lambda_C_1x) * DMat_dR1x_dgma
                                  + trans(Lambda_C_1y) * DMat_dR1y_dgma
                                  + trans(Lambda_C_1z) * DMat_dR1z_dgma
                                  + trans(Lambda_C_2)  * Mat_dR2_dgma
                                  + trans(Lambda_C_3F) * Mat_dR3F_dgma;                                  
    }
else
    {
    DLJ_Dgma = trans(Vec_dJ_dgma) + trans(Lambda_J_1x) * DMat_dR1x_dgma
                                  + trans(Lambda_J_1y) * DMat_dR1y_dgma
                                  + trans(Lambda_J_2)  * Mat_dR2_dgma
                                  + trans(Lambda_J_3F) * Mat_dR3F_dgma;

    DLC_Dgma = trans(Vec_dC_dgma) + trans(Lambda_C_1x) * DMat_dR1x_dgma
                                  + trans(Lambda_C_1y) * DMat_dR1y_dgma
                                  + trans(Lambda_C_2)  * Mat_dR2_dgma
                                  + trans(Lambda_C_3F) * Mat_dR3F_dgma;                                  
    }

// ~~~~~~~~~~ Save fields:
for(label i=0; i<NC; ++i)
    {
    DFy_Dgma[i]    = DLJ_Dgma(0, i);
    DFx_Dgma[i]    = DLC_Dgma(0, i);    
    }
    
Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;
Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ENDING_OAD_Section~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;
Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl; 
