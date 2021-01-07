// ----------------------- //
// Misc properties
// ----------------------- //

	// Value of Pi
    	const double pi = 3.1415926535897;

    	// Boltzmann constant
    	const dimensionedScalar k_B("k_B", dimensionSet(1,2,-2,-1,0,0,0), 1.38064852e-23);

	// Molecular (note: not molar) mass of H2O 
    	const dimensionedScalar m_gas("m_gas", dimensionSet(1,0,0,0,0,0,0), 2.988E-26);	

	// Minimum physically meaningful radius of a H2O liduid (see documentation)
    	const dimensionedScalar r_droplet_minimum("r_droplet_minimum", dimensionSet(0,1,0,0,0,0,0), 2.0*Foam::pow(10.0, -10.0));	

	// Specific gas constant
	const dimensionedScalar R("R", dimensionSet(0,2,-2,-1,0,0,0), 461.52);

	// Properties at the critical point of water
	const dimensionedScalar T_crit("T_crit", dimensionSet(0,0,0,1,0,0,0), 647.096); // K
	const dimensionedScalar p_crit("p_crit", dimensionSet(1,-1,-2,0,0,0,0), 22064000.0); // Pa
	const dimensionedScalar rho_crit("rho_crit", dimensionSet(1,-3,0,0,0,0,0), 322.0); // kg/m3

	// Vukalovich terms
	const scalar ee = 63.2;  // Called "e" in the paper, renamed here to resolve variable naming conflicts
	const scalar bb = 0.00085;
	const scalar C  = 0.3977E+06;
	const scalar G  = 47.053;
	const scalar m1 = 1.968;
	const scalar m2 = 2.957;
	const scalar w1 = (3.0+2.0*m1)/2.0;
	const scalar w2 = (3.0*m2-4.0*m1)/2.0;
	const scalar vukal_aux_term3 = C*C*G*G;

// ----------------------- //
// Saturation properties
// ----------------------- //

	// Saturation pressure
	    const scalar psat_0 = -7.85951783;
	    const scalar psat_1 = 1.84408259;
	    const scalar psat_2 = -11.7866497;
	    const scalar psat_3 = 22.6807411;
	    const scalar psat_4 = -15.9618719;
	    const scalar psat_5 = 1.80122502;

	// Saturation temperature
	    const scalar Tsat_0 = 0.11670521452767E+04;
	    const scalar Tsat_1 = -0.72421316703206E+06;
	    const scalar Tsat_2 = -0.17073846940092E+02;
	    const scalar Tsat_3 = 0.12020824702470E+05;
	    const scalar Tsat_4 = -0.32325550322333E+07;
	    const scalar Tsat_5 = 0.14915108613530E+02;
	    const scalar Tsat_6 = -0.48232657361591E+04;
	    const scalar Tsat_7 = 0.40511340542057E+06;
	    const scalar Tsat_8 = -0.23855557567849;
	    const scalar Tsat_9 = 0.65017534844798E+03;

	// Saturation density
	    const scalar rhosat_0 = -2.03150240;
	    const scalar rhosat_1 = -2.68302940;
	    const scalar rhosat_2 = -5.38626492;
	    const scalar rhosat_3 = -17.2991605;
	    const scalar rhosat_4 = -44.7586581;
	    const scalar rhosat_5 = -63.9201063;

	// Saturated enthalpy
    	    const dimensionedScalar hsat_ref("hsat_ref", dimEnergy/dimMass, 1000.0);
	    const scalar hsat_0 = -1135.905627715;
	    const scalar hsat_1 = -5.65134998e-08;
	    const scalar hsat_2 = 2690.66631;
	    const scalar hsat_3 = 127.287297;
	    const scalar hsat_4 = -135.003439;
	    const scalar hsat_5 = 0.981825814;

// ----------------------- //
// Liquid properties
// ----------------------- //

	// Liquid density
	    const scalar  rhol_0 = 1.99274064;
	    const scalar  rhol_1 = 1.09965342;
	    const scalar  rhol_2 = - 0.510839303;
	    const scalar  rhol_3 = - 1.75493479;
	    const scalar  rhol_4 = - 45.5170352;
	    const scalar  rhol_5 = - 6.74694450E+05;

	// Liquid surface tension
	    const scalar sigma_B = 0.2358;
	    const scalar sigma_b = -0.625;
	    const scalar sigma_mu = 1.256;

// ----------------------- //
// Dimensioned factors that are later used to remove dimensionality from some thermophysical properties
// This is done because "for" loops are used in the code and:
//	1) OpenFOAM is sensitive to dimensional consistency (and rightly so)
//	2) OpenFOAM does not carry a parameter's dimensionality into a loop - this is inherent to how the software (and C++) is structured
// So the factors below are used to nondimensionalise some parameters that are used in loops
// ----------------------- //

    dimensionedScalar unit_Kelvin("unit_Kelvin", dimTemperature, 1.0); // To resolve dimensional inconsistencies
    dimensionedScalar unit_MPa("unit_MPa", dimPressure, 1000000.0); // To resolve dimensional inconsistencies
    dimensionedScalar unit_specific_heat("unit_specific_heat", dimensionSet(0, 2, -2, -1, 0, 0, 0), 1.0); // To resolve dimensional inconsistencies
