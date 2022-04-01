// Momentum predictor
fvVectorMatrix UEqn   //use tmp<fvVectorMatrix> for better performance
(
    fvm::div(phi, U)
    -fvm::laplacian(nu, U)
    +fvm::Sp(gamma/etta, U)
 ==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

solve(UEqn == -fvc::grad(p));

fvOptions.correct(U);

