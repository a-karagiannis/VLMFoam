// ------------------------------------------------------------------------------------------------------------- //
// --- ENERGY
// ------------------------------------------------------------------------------------------------------------- //

surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
      + fvc::dotInterpolate(mesh.Sf(), tauMC)
    )
  & (a_pos*U_pos + a_neg*U_neg)
);

solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  - fvc::div(sigmaDotU)
  ==
  Source_e
);

// Update the energy variable for the new density
e = max(rhoE/rho - 0.5*magSqr(U), Emin);
e.correctBoundaryConditions();
thermo.correct();
rhoE.boundaryFieldRef() ==
    rho.boundaryField()*
    (
        e.boundaryField() + 0.5*magSqr(U.boundaryField())
    );

if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
    );
    thermo.correct();
    rhoE = rho*(e + 0.5*magSqr(U));
}

