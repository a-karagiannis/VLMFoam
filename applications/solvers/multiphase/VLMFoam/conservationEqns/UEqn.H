// ------------------------------------------------------------------------------------------------------------- //
// --- MOMENTUM
// ------------------------------------------------------------------------------------------------------------- //
volVectorField Source_U("Source_U", -U*Source_Y);
volVectorField Source_U_linear_explicit("Source_U_linear_explicit", SourceMomentumFactor*rhoU);

solve
(
	fvm::ddt(rhoU) 
	+ fvc::div(phiUp) 
	== 
	Source_U 
	- Source_U_linear_explicit 
	+ fvm::SuSp(SourceMomentumFactor, rhoU)
);

// Update U value for the new density
U.ref() =
    rhoU()
   /rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

// Correct the U value for diffusive transport, if the flow is viscous
if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(muEff, U)
      - fvc::div(tauMC)
    );
    rhoU = rho*U;
}

