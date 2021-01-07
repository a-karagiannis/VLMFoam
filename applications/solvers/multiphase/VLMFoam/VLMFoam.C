/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    VLMFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor, augmented to account for the non-equilibrium
    homogeneous condensation of the gaseous phase, as modelled for steam
    by Gerber and Kermani (2004).

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "readUserInputs/readNucleationCoefficients.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readUserInputs/readFluxScheme.H"

    #include "phaseChangeThermodynamics/constantsForEmpiricalEqns.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    // Limit the minimum energy to prevent numerical instabilities at the start
    #include "phaseChangeThermodynamics/boundEnergy.H"

    // Store cell volumes (useful later in the code)
    cellVolume.ref() = mesh.V(); // Not applicable to cells at boundaries


    Info<< "\nStarting time loop\n" << endl;


    while (runTime.run())
    {

        #include "interpolateFields.H"

        // -----------------------------------------------------------
        // --- Thermophysical properties
        // -----------------------------------------------------------
        Cp = thermo.Cp();              // Specific heat at constant pressure
        Cv = thermo.Cv();              // Specific heat at constant volume
        gamma = thermo.gamma();        // Ratio of specific heats
        c = Foam::sqrt(gamma*rPsi);    // Speed of sound

        // Evaluate convective fluxes
        #include "computeFluxes.H"

        // -----------------------------------------------------------
        // --- Transport properties
        // -----------------------------------------------------------

        // Effective viscosity (laminar + turbulent contribution)
        volScalarField muEff("muEff", turbulence->muEff());

        // Effective thermal conductivity
        volScalarField kappaEff("kappaEff", turbulence->kappaEff());

        // Navier-Stokes stress tensor
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        // SOLVE CONSERVATION LAWS * * * * * * * * * * * * * * * * * //

        // Auxilliary terms for the source terms
        SourceMomentumFactor = -(Source_Y_nucleation +
                                  Source_Y_growth_active_coeff*rhoY)/rho;

        Source_e = -h_droplet*(Source_Y_nucleation +
                                Source_Y_growth_active_coeff*rhoY);

        // Mass conservation
        #include "conservationEqns/rhoEqn.H"

        // Momentum conservation
        #include "conservationEqns/UEqn.H"

        // Energy conservation
        #include "conservationEqns/EEqn.H"

        // Droplet number conservation
        #include "conservationEqns/NEqn.H"

        // Liquid mass fraction conservation
        #include "conservationEqns/YEqn.H"

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // ---------------------------------------------------------------
        // --- Update thermodynamic properties
        // ---------------------------------------------------------------

        p.ref() = rho()/psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

        h = e + p/rho; // Specific enthalpy of gas

        turbulence->correct();

        // Transport properties
          // Plain ones, not "effective" ones as above, because the
          // nucleation model does not concern itself with turbulence artefacts
        volScalarField kappa("kappa", thermo.kappa());	// Thermal conductivity
        volScalarField mu("mu", thermo.mu());         	// Dynamic viscosity

        // Flow similarity parameters
        volScalarField Prandtl("Prandtl", thermo.Cp()*thermo.mu()
                                          /thermo.kappa());	// Prandtl number



        // -------------------------------------------------------------- //
        // -------------------------------------------------------------- //
      	// --- UPDATE THE CONSERVATION SOURCE TERMS
      	// -------------------------------------------------------------- //
        // -------------------------------------------------------------- //

        // --- Loop over all cells in the internal field
        // (i.e. cells at the boundaries not included, coming later)
        forAll(T, cellI)
        {
          #include "phaseChangeThermodynamics/saturationProperties.H"

          // Since the conservation equations for N and Y are not directly
          // coupled, it may happen numerically that N>0 while Y=0 (or, rarely,
          // vice-versa).
          // To remedy this, here, we delete N or Y if there is an
          // inconsistency between them
          if (N[cellI]*rho[cellI]*cellVolume[cellI] > 1.0
                          && Y[cellI] > Foam::pow(10.0, -12.0))
          {
          }
          else
          {
            N[cellI] = 0.0;
            Y[cellI] = 0.0;
          }

          // Only bother with condensation checks if temperature below...
              // ...both the critical value and saturation value and also
              // higher than a minimum threshold of 273.15 for the
              // thermophysical relations
          if (T[cellI] < T_crit.value()
              && T[cellI] < Tsat[cellI]
              && T[cellI] >= 173.16)
          {
            // --- UPDATE LIQUID-RELATED THERMODYNAMIC PROPERTIES IN THE CELL

                // External file for liquid-related properties
                #include "phaseChangeThermodynamics/liquidThermophysicalProperties.H"

        //******************************************************************//
        //******************************************************************//

                // Only bother with condensation calculations if there is...
                // ...at least one liquid droplet in the current control volume
                if (N[cellI]*rho[cellI]*cellVolume[cellI] > 1.0
                                && Y[cellI] > Foam::pow(10.0, -12.0))
                {
                    // Mean radius resulting from existing droplets in the cell
                    r_droplet_actual[cellI] = Foam::pow((3.0*Y[cellI]/
                                              4.0/pi/rho_liquid[cellI]/
                                              N[cellI]), 1.0/3.0);

                    // Total interfacial surface area
                    beta[cellI] = N[cellI] *
                                  4.0*pi*r_droplet_actual[cellI]*
                                  r_droplet_actual[cellI];
                }
                else // Delete the droplets in the cell
                {
                    r_droplet_actual[cellI] = 0.0;
                    beta[cellI]             = 0.0;
                    N[cellI]                = 0.0;
                    Y[cellI]                = 0.0;
                }

                // Only bother with nucleation if saturation ratio > 1
                if (S_sat[cellI] > 1.0 + SMALL)
                {
                  // Change in Gibbs free energy
                  DeltaG[cellI] =  R.value()*T[cellI]*(Foam::log(S_sat[cellI]));

                  // Critical radius
                  r_droplet_critical[cellI] = 2.0*surface_tension[cellI]/
                                                rho_liquid[cellI]/DeltaG[cellI];
                }
                else // No supersaturation, set radius to 0, for convenience
                {
                    r_droplet_critical[cellI] = 0.0;
                }

                if (r_droplet_actual[cellI] >= r_droplet_critical[cellI]
                  && r_droplet_critical[cellI] > r_droplet_minimum.value())
                {

                  // Droplet temperature
                  T_droplet[cellI] = Tsat[cellI] - (Tsat[cellI] - T[cellI]) *
                                      r_droplet_critical[cellI]/
                                      r_droplet_actual[cellI];

                  // Correction factor for the coefficient of
                  // thermal convection from gas to droplet
                  v_corr[cellI] = R.value()*Tsat[cellI]/h_fg[cellI] *
                                  (alpha_coeff - 0.5 - (2.0 - q_c)/(2.0*q_c) *
                                  (gamma[cellI]+1.0)/(2.0*(gamma[cellI] - 1.0))*
                                  R.value()*Tsat[cellI]/h_fg[cellI]);
                  if (v_corr[cellI] >= 1.0) // For numerical artifacts
                  {
                       v_corr[cellI] = 1.0;
                  }

                  // Mean free path
                  mean_free_path[cellI] = 1.5*mu[cellI]/p[cellI] *
                                          Foam::sqrt(R.value()*T[cellI]);

                  // Knudsen number based on the droplet diameter
                  Knudsen_droplet[cellI] = mean_free_path[cellI]/
                                            (2.0*r_droplet_actual[cellI]);

                  // Coefficient of thermal convection from gas to droplet
                  lambda_g[cellI] = kappa[cellI] *1.0/
                                    (r_droplet_actual[cellI]
                                          + (1.0-v_corr[cellI])*
                                    (2.0*Foam::sqrt(8.0*pi)/1.5)*
                                    (gamma[cellI]/(gamma[cellI] + 1.0))*
                                    mean_free_path[cellI]/2.0/Prandtl[cellI]);

                  // Rate of droplet growth
                  drdt[cellI] = (lambda_g[cellI] *
                                (T_droplet[cellI] - T[cellI])) /
                                (h_fg[cellI]*rho_liquid[cellI]);

                  Source_Y_growth_active_coeff[cellI]   =
                  3.0*drdt[cellI]*1.0/r_droplet_actual[cellI];
                }
                else
                {
                  drdt[cellI] = 0.0;
                  Source_Y_growth_active_coeff[cellI] = 0.0;
                }


        //****************************************************************//
        //****************************************************************//

                if (r_droplet_critical[cellI] > r_droplet_minimum.value())
                {
                  // Kantrowitz correction factor for the nucleation rate
                  etaKantrowitz[cellI] = 2.0*(gamma[cellI]-1.0)/
                                          (gamma[cellI]+1.0)*
                                          (h_fg[cellI]/(R.value()*T[cellI]))*
                                          (h_fg[cellI]/
                                            (R.value()*T[cellI]) - 0.5);
                  if (etaKantrowitz[cellI] < 0.0 )
                  {
                    Info << "WARNING: Kantrowitz correction factor is negative, assumed zero" << endl;
                    etaKantrowitz[cellI] = 0.0;
                  }

                  // Nucleation rate of critical radius droplets,
                  // per unit volume of vapour
                  J[cellI] = q_c/(1.0 + etaKantrowitz[cellI])*
                              Foam::sqrt(2.0*surface_tension[cellI]/pi/
                                    m_gas.value()/m_gas.value()/m_gas.value())*
                              rho[cellI]*rho[cellI]/rho_liquid[cellI]*
                              Foam::exp(-4.0*pi*r_droplet_critical[cellI]*
                                r_droplet_critical[cellI]*surface_tension[cellI]
                                /3.0/k_B.value()/T[cellI]);

                  // Correction by Wolk and Strey
		              // J[cellI] = J[cellI]*Foam::exp(-27.56 + 6500.0/T[cellI]);

                  // If there is at least 1 critically sized molecule generated
                  // in the cell volume over this timestep
                  if (J[cellI]*cellVolume[cellI]*
                        runTime.time().deltaT().value() > 1.0)
                  {
                  }
                  else
                  {
                    J[cellI] = 0.0;
                  }
                }
                else // Value of critical radius unphysically small
                {
                  r_droplet_critical[cellI] = 0.0;
                  J[cellI] = 0.0;
                }

        //******************************************************************//
        //*****************************************************************//

                // Contribution to the liquid mass from the
                // nucleation of new droplets
               Source_Y_nucleation[cellI] = J[cellI]*rho_liquid[cellI]*
                                                4.0/3.0*pi*
                                                r_droplet_critical[cellI]*
                                                r_droplet_critical[cellI]*
                                                r_droplet_critical[cellI];

                // Contribution to the liquid mass from the
                // growth of existing droplets
                Source_Y_growth[cellI]    = rho_liquid[cellI]*beta[cellI]
                                                  *drdt[cellI]*rho[cellI];

                // Total (nucleation + growth) liquid mass source term
                Source_Y[cellI]                = Source_Y_growth[cellI] +
                                                  Source_Y_nucleation[cellI];

              }
              else // No supercooling, relevant source terms set to zero
              {
                  J[cellI]                              = 0.0;
                  Source_Y[cellI]                       = 0.0;
                  Source_Y_nucleation[cellI]            = 0.0;
                  Source_Y_growth[cellI]                = 0.0;
                  Source_Y_growth_active_coeff[cellI]   = 0.0;
                  SourceMomentumFactor[cellI]           = 0.0;
                  Source_e[cellI]                       = 0.0;
              }         // End of loop with the tempperature checks
        }

        // ---------------------------------------------------------------- //
        // -------------------- END CONDENSATION MODEL -------------------- //
        // ---------------------------------------------------------------- //

        // -------------------------------------------------------------------
        // --- Update timestep
        // -------------------------------------------------------------------

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
