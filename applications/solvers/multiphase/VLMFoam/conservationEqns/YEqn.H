// --------------------------------------------------------------- //
// LIQUID MASS FRACTION //
// --------------------------------------------------------------- //

volScalarField Source_Y_growth_linear_explicit("Source_Y_growth_linear_explicit", Source_Y_growth_active_coeff*rhoY);
solve
    (
	fvm::ddt(rhoY) 
	+ fvc::div(phiY) 
	== 
	Source_Y 
	- Source_Y_growth_linear_explicit 
	+ fvm::SuSp(Source_Y_growth_active_coeff, rhoY)
    );

Y.ref() =
    rhoY()
   /rho();

Y.correctBoundaryConditions();

rhoY.boundaryFieldRef() == rho.boundaryField()*Y.boundaryField();
