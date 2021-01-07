solve(
fvm::ddt(rho) 
+ fvc::div(phi) 
== 
-(
Source_Y_nucleation 
+ Source_Y_growth_active_coeff*rhoY)
);

