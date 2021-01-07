if (T[cellI] >= 273.16)
{

	scalar Tratio	 = min(1.0, T[cellI]/T_crit.value()); // min used to avoid ratios >1
	scalar Tterm	 = max(0.0, 1.0 - T[cellI]/T_crit.value()); // max used to avoid negatives
	scalar invTterm  = max(0.0, T_crit.value()/T[cellI] - 1.0); // max used to avoid negatives


// Saturation pressure at local temperature [Pa]

    	psat[cellI] = p_crit.value()*Foam::exp(T_crit.value()/T[cellI]*
				(psat_0*Foam::pow(Tterm, 1.0) 
				+ psat_1*Foam::pow(Tterm, 1.5) 
				+ psat_2*Foam::pow(Tterm, 3.0) 
				+ psat_3*Foam::pow(Tterm, 3.5) 
				+ psat_4*Foam::pow(Tterm, 4.0) 
				+ psat_5*Foam::pow(Tterm, 7.5)));

// Saturation temperature at local pressure [K]
	scalar beta_Ts	 = Foam::pow(p[cellI]/unit_MPa.value(), 0.25);
	scalar E_Ts	 = beta_Ts*beta_Ts + Tsat_2*beta_Ts + Tsat_5;
	scalar F_Ts	 = Tsat_0*beta_Ts*beta_Ts + Tsat_3*beta_Ts + Tsat_6;
	scalar G_Ts	 = Tsat_1*beta_Ts*beta_Ts + Tsat_4*beta_Ts + Tsat_7;
	scalar D_Ts	 = 2.0*G_Ts/(-F_Ts - Foam::pow((F_Ts*F_Ts - 4.0*E_Ts*G_Ts), 0.5));
    
	Tsat[cellI] = ((Tsat_9 + D_Ts - Foam::pow(((Tsat_9 + D_Ts)*(Tsat_9 + D_Ts) - 4.0*(Tsat_8 + Tsat_9*D_Ts)), 0.5))/2.0);

// Saturation density at local temperature [kg/m3]
    	rhosat[cellI] = rho_crit.value()*
				Foam::exp(rhosat_0*Foam::pow(Tterm, 2.0/6.0) 
				+ rhosat_1*Foam::pow(Tterm, 4.0/6.0) 
				+ rhosat_2*Foam::pow(Tterm, 8.0/6.0)
				+ rhosat_3*Foam::pow(Tterm, 18.0/6.0) 
				+ rhosat_4*Foam::pow(Tterm, 37.0/6.0) 
				+ rhosat_5*Foam::pow(Tterm, 71.0/6.0));

// Enthalpy at saturated temperature [J/kg]
    	hsat[cellI] = hsat_ref.value()*(hsat_0 
				+ hsat_1*Foam::pow(Tratio, -19.0) 
				+ hsat_2*Tratio
				+ hsat_3*Foam::pow(Tratio, 4.5) 
				+ hsat_4*Foam::pow(Tratio, 5.0) 
				+ hsat_5*Foam::pow(Tratio, 54.5));

}
else // Extra temperature check for T > 173.16 is added in the main code
{
	// Saturation properties of water vapour in the metastable region below the triple point
    	scalar Tratio	 = min(1.0, T[cellI]/T_crit.value()); // min used to avoid ratios >1
	scalar Tterm	 = max(0.0, 1.0 - T[cellI]/T_crit.value()); // max used to avoid negatives

	// Saturation pressure at local temperature [Pa]
	psat[cellI] = Foam::exp(
			54.842763 - 6763.22/T[cellI] 
			- 4.210*Foam::log(T[cellI]) + 0.000367*T[cellI]
			+ Foam::tanh(0.0415*(T[cellI] - 218.8))
			* (53.878 - 1331.22/T[cellI] - 9.44523*Foam::log(T[cellI]) + 0.014025*T[cellI])
			);

	// Saturation temperature at local pressure [K]
	Tsat[cellI] =  2.08395627e+02
			+ 7.59174915e+00*Foam::log(p[cellI])
			+ 2.88687479e-01*Foam::pow(Foam::log(p[cellI]), 2.0)
			+ 2.03780654e-03*Foam::pow(Foam::log(p[cellI]), 3.0)
			+ 2.02354582e-03*Foam::pow(Foam::log(p[cellI]), 4.0);

	// Saturated density at local temperature [K]
	rhosat[cellI] = rho_crit.value()*Foam::exp(rhosat_0*Foam::pow(Tterm, 2.0/6.0) + rhosat_1*Foam::pow(Tterm, 4.0/6.0) + rhosat_2*Foam::pow(Tterm, 8.0/6.0)
			+ rhosat_3*Foam::pow(Tterm, 18.0/6.0) + rhosat_4*Foam::pow(Tterm, 37.0/6.0) + rhosat_5*Foam::pow(Tterm, 71.0/6.0));
}

