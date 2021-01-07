	// Auxilliary nondimensional term
	scalar Tterm	 = max(0.0, 1.0 - T[cellI]/T_crit.value()); // max used to avoid negatives

if (T[cellI] > 273.16)
{

// ----------------------- //
// Density
// ----------------------- //

    	rho_liquid[cellI] = rho_crit.value()*(1.0 + rhol_0*Foam::pow(Tterm, 1.0/3.0) 
				+ rhol_1*Foam::pow(Tterm, 2.0/3.0) 
				+ rhol_2*Foam::pow(Tterm, 5.0/3.0) 
				+ rhol_3*Foam::pow(Tterm, 16.0/3.0) 
				+ rhol_4*Foam::pow(Tterm, 43.0/3.0) 
				+ rhol_5*Foam::pow(Tterm, 110.0/3.0));

// ----------------------- //
// Surface tension
// ----------------------- //

    	surface_tension[cellI] = sigma_B*Foam::pow(Tterm, sigma_mu)*(1.0 + sigma_b*Tterm);

    	dsdT[cellI] = sigma_B*sigma_mu*Foam::pow(Tterm, sigma_mu-1.0)*(-1.0/T_crit.value())*(1.0 + sigma_b*Tterm) 
							+ sigma_B*Foam::pow(Tterm, sigma_mu)*(-sigma_b/T_crit.value());


// ----------------------- //
// Total specific enthalpy
// ----------------------- //

	// From IAPWS
	dpsdTs[cellI] = - psat[cellI]/T[cellI]*(Foam::log(psat[cellI]/p_crit.value()) 
						+ psat_0 + 1.5*psat_1*Foam::pow(Tterm, 0.5) 
						+ 3.0*psat_2*Foam::pow(Tterm, 2.0) 
						+ 3.5*psat_3*Foam::pow(Tterm, 2.5) 
						+ 4.0*psat_4*Foam::pow(Tterm, 3.0) 
						+ 7.5*psat_5*Foam::pow(Tterm, 6.5));

    	h_liquid[cellI] = hsat[cellI] + 1.0/rho_liquid[cellI]*T[cellI]*dpsdTs[cellI]; //  - 1.0/rho_liquid[cellI] - We have ignored the specific volume of water (v_g >> v_f)


	// Latent heat release upon droplet condensation
	h_fg[cellI] = h[cellI] - h_liquid[cellI];

	// Degree of supercooling
	T_sc[cellI] = Tsat[cellI] - T[cellI];

	// Supersaturation ratio
	S_sat[cellI] = p[cellI]/psat[cellI];

	// Droplet enthalpy (contribution of surface energy is commented out because it is very small)
	h_droplet[cellI] = h_liquid[cellI] + 0.5*magSqr(U[cellI]);
		  	// +3.0/rho_liquid[cellI]/r_droplet_actual[cellI]*(surface_tension[cellI] 
			// - T_droplet[cellI]*dsdT[cellI])
}
else if (T[cellI] < 273.16 && T[cellI] >= 173.16) // Different set of thermophysical properties for temperature range below the freezing point
{						    // Only (presumably) valid for extremely rapid expansions (e.g., supersonic micronozzles)
	rho_liquid[cellI] 	= (0.08*Foam::tanh((T[cellI]-225.0)/46.2) 
					+ 0.7415*Foam::pow(((T_crit.value() 
					- T[cellI])/T_crit.value()), 0.33) + 0.32)*1000.0;
	surface_tension[cellI] 	= 96.6635e-03 + 9.133e-06*T[cellI] - 2.75e-07*T[cellI]*T[cellI]; 
	dsdT[cellI] 		= 9.133e-06 - 5.5e-07*T[cellI];
	h_fg[cellI] 		= (
				   6763.22/T[cellI]/T[cellI] - 4.21/T[cellI] + 0.000367
				   + (0.0415*(53.878 - 1331.22/T[cellI] + 0.014025*T[cellI] - 9.44523*Foam::log(T[cellI])))/
					((Foam::cosh(0.0415*(T[cellI] - 218.8))) * (Foam::cosh(0.0415*(T[cellI] - 218.8))))
				   + (0.014025 + 1331.22/T[cellI]/T[cellI] - 9.44523/T[cellI])*Foam::tanh(0.0415*(T[cellI] - 218.8))
				  ) 
					* R.value() * T[cellI] * T[cellI]; 
	h_liquid[cellI] 	=  h[cellI] - h_fg[cellI]; // Below the freezing point, this is by convention negative

	// Degree of supercooling
	T_sc[cellI] = Tsat[cellI] - T[cellI];

	// Supersaturation ratio
	S_sat[cellI] = p[cellI]/psat[cellI];

	// Droplet enthalpy (contribution of surface energy is commented out because it is very small)
	h_droplet[cellI] = h_liquid[cellI] + 0.5*magSqr(U[cellI]);
		  	// +3.0/rho_liquid[cellI]/r_droplet_actual[cellI]*(surface_tension[cellI] 
			// - T_droplet[cellI]*dsdT[cellI])
}
