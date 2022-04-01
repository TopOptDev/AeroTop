// ~~~~~~~~~~ create a copy of p:
volScalarField p_tmp = p;

// ~~~~~~~~~~ Momentum solver (last step)
Info << nl << "~~~~~~~~~~ Solving for U (final step): " << nl << endl;
fvVectorMatrix UEqn_ref
(
    fvm::div(phi, U)
    -fvm::laplacian(nu, U)
    +fvm::Sp(gamma/etta, U)
);
solve(UEqn_ref == -fvc::grad(p));

// ~~~~~~~~~~ Pressure solver (last step)
Info << nl << "~~~~~~~~~~ Solving for p2 (final step): " << nl << endl;
volScalarField rAU_ref(1.0/UEqn_ref.A());
volVectorField HbyA_ref(constrainHbyA(rAU_ref*UEqn_ref.H(), U, p)); 
surfaceScalarField phiHbyA_ref("phiHbyA_ref", fvc::flux(HbyA_ref));

fvScalarMatrix pEqn_ref
(
    fvm::laplacian(rAU_ref, p) == fvc::div(phiHbyA_ref)
);
pEqn_ref.setReference(pRefCell, pRefValue);
pEqn_ref.solve();

// ~~~~~~~~~~ update p/p2 values 
p2 = p;
p = p_tmp;

