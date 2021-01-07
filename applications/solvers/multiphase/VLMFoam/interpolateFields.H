        // -----------------------------------------------------------
        // --- Directed interpolation of primitive fields onto faces
        // -----------------------------------------------------------
        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        surfaceScalarField rhoY_pos(interpolate(rhoY, pos, Y.name()));
	      surfaceScalarField rhoY_neg(interpolate(rhoY, neg, Y.name()));

        surfaceScalarField rhoN_pos(interpolate(rhoN, pos, N.name()));
	      surfaceScalarField rhoN_neg(interpolate(rhoN, neg, N.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField Y_pos("Y_pos", rhoY_pos/rho_pos);
	      surfaceScalarField Y_neg("Y_neg", rhoY_neg/rho_neg);

        surfaceScalarField N_pos("N_pos", rhoN_pos/rho_pos);
	      surfaceScalarField N_neg("N_neg", rhoN_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
