
	// To accelerate computations by reducing the number of times that the funtion pow() is executed, 
	// we use the auxilliary terms below and express the exponents as in the example below:
	// T^{-2.0w1 - 2.0} -> aux_term2*aux_term2*(1.0/T/T)

	volScalarField aux_term1("aux_term1", pow(Tdimless, -w1));
	volScalarField aux_term2("aux_term2", pow(Tdimless, -w2));

	volScalarField phi1("phi1", C*G/(pow(Tdimless, w1)));
	volScalarField phi2("phi2", 1.0 -22.7/(pow(Tdimless, w2)));

	// Virial coefficients in the formula of Vukalovich
	const scalar B1 = 1.0;

	volScalarField B2("B2", -ee/(G*Tdimless) - phi1 + bb);

	volScalarField B3("B3", -bb*phi1 + 4.0*phi1*phi1*phi2);

	volScalarField B4("B4", 32.0*bb*phi1*phi1*phi2);
/*
	// First derivatives of coefficients w.r.t. temperature T
	volScalarField dB2("dB2", ee/G*(1.0/Tdimless/Tdimless) + C*G*w1*aux_term1*(1.0/Tdimless));

	volScalarField dB3("dB3", bb*w1*C*G*aux_term1*(1.0/Tdimless) - 8.0*w1*vukal_aux_term3*aux_term1*aux_term1*(1.0/Tdimless)*(1.0 - 22.7*aux_term2)
					+ 4.0*vukal_aux_term3*aux_term1*aux_term1*(22.7*w2*aux_term2*(1.0/Tdimless)));

	volScalarField dB4("dB4", -32.0*bb*vukal_aux_term3*2.0*w1*aux_term1*aux_term1*(1.0/Tdimless)*(1.0-22.7*aux_term2) 
					+ 32.0*bb*vukal_aux_term3*aux_term1*aux_term1*(22.7*w2*aux_term2*(1.0/Tdimless)));

	// Second derivatives of coefficients w.r.t. temperature T	
	volScalarField d2B2("d2B2", -2.0*ee/G*(1.0/Tdimless/Tdimless/Tdimless) - C*G*w1*(w1 + 1.0)*aux_term1*(1.0/Tdimless/Tdimless));

	volScalarField d2B3("d2B3", -bb*w1*(w1 + 1.0)*C*G*aux_term1*(1.0/Tdimless/Tdimless)
					+ 8.0*w1*(2.0*w1 + 1.0)*vukal_aux_term3*aux_term1*aux_term1*(1.0/Tdimless/Tdimless)*(1.0 - 22.7*aux_term2)
					- 2.0*8.0*w1*vukal_aux_term3*aux_term1*aux_term1*(1.0/Tdimless)*(22.7*w2*aux_term2*(1.0/Tdimless))
					- 4.0*vukal_aux_term3*aux_term1*aux_term1*(22.7*w2*(w2+1.0)*aux_term2*(1.0/Tdimless/Tdimless)));

	volScalarField d2B4("d2B4", + 32.0*bb*vukal_aux_term3*2.0*w1*(2.0*w1 + 1.0)*aux_term1*aux_term1*(1.0/Tdimless/Tdimless)*(1.0 - 22.7*aux_term2)
					- 2.0*32*bb*vukal_aux_term3*2.0*w1*aux_term1*aux_term1*(1.0/Tdimless)*(22.7*w2*aux_term2*(1.0/Tdimless))
					- 32.0*bb*vukal_aux_term3*aux_term1*aux_term1*(22.7*w2*(w2+1.0)*aux_term2*(1.0/Tdimless/Tdimless)));
*/
