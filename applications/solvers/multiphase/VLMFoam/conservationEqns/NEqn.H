        // ------------------------------------------------------------------------------------------------------------------
        // --- DROPLET NUMBER
        // ------------------------------------------------------------------------------------------------------------------

solve(fvm::ddt(rhoN) + fvc::div(phiN) == J);
N.ref() =
    rhoN()
   /rho();
N.correctBoundaryConditions();
rhoN.boundaryFieldRef() == rho.boundaryField()*N.boundaryField();
