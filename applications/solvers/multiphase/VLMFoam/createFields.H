#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& e = thermo.he();
volScalarField& p = thermo.p();

///////
// --- Fields for the gaseous phase --- //
///////

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


///////
// --- Fields for the liquid phase --- //
///////

volScalarField N	// Number of droplets per unit vapour mass
(
    IOobject
    (
        "N",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Y	// Liquid mass fraction
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rhoN
(
    IOobject
    (
        "rhoN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*N
);

volScalarField rhoY
(
    IOobject
    (
        "rhoY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Y
);

volVectorField rhoUN
(
    IOobject
    (
        "rhoUN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*N
);

volVectorField rhoUY
(
    IOobject
    (
        "rhoUY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*Y
);

surfaceScalarField phiN("phiN", fvc::flux(rhoUN));

surfaceScalarField phiY("phiY", fvc::flux(rhoUY)); 

	// Specific heats:
volScalarField Cp
	(
	    IOobject
	    (
		"Cp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	mesh,
	dimensionedScalar("Cp", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0)
);


volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("Cv", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0)
);

volScalarField gamma
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("gamma", dimless, 0)
);

volScalarField h // Enthalpy
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	e + p/rho
);

volScalarField hsat // Enthalpy at saturation conditions
(
    IOobject
    (
        "hsat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("hsat", dimEnergy/dimMass, 0)
);

volScalarField DeltaG	// Change in Gibbs free energy
(
    IOobject
    (
        "DeltaG",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("DeltaG", dimEnergy/dimMass, 0)
);

volScalarField c	// Speed of sound
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("c", dimLength/dimTime, 0)
);

// ----------------------------- //
// --- Saturation properties --- //
// ----------------------------- //

volScalarField psat	
(
    IOobject
    (
        "psat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("psat", dimPressure, 0)
);

volScalarField Tsat	
(
    IOobject
    (
        "Tsat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tsat", dimTemperature, 0)
);

volScalarField rhosat	
(
    IOobject
    (
        "rhosat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhosat", dimDensity, 0)
);

volScalarField S_sat	
(
    IOobject
    (
        "S_sat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("S_sat", dimless, 0)
);

volScalarField T_sc	
(
    IOobject
    (
        "T_sc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T_sc", dimTemperature, 0)
);

volScalarField dpsdTs	
(
    IOobject
    (
        "dpsdTs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dpsdTs", dimPressure/dimTemperature, 0)
);

// ------------------------------------- //
// Other condensation related properties
// ------------------------------------- //

volScalarField h_fg // Change in enthalpy upon condensation
(
    IOobject
    (
	"h_fg",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("h_fg", dimEnergy/dimMass, 0)
);

volScalarField r_droplet_critical // Critical droplet radius
(
    IOobject
    (
	"r_droplet_critical",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("r_droplet_critical", dimLength, 0)
);

volScalarField r_droplet_actual // Current droplet radius
(
    IOobject
    (
	"r_droplet_actual",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("r_droplet_actual", dimLength, 0)
);

volScalarField T_droplet // Temperature of the droplet
(
    IOobject
    (
	"T_droplet",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("T_droplet", dimTemperature, 0)
);

volScalarField mean_free_path // Mean free path of the gas at local conditions
(
    IOobject
    (
	"mean_free_path",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mean_free_path", dimLength, 0)
);

volScalarField Knudsen_droplet // Knudsen number based on the diameter of the droplet
(
    IOobject
    (
	"Knudsen_droplet",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Knudsen_droplet", dimless, 0)
);

volScalarField lambda_g // Coefficient of convective heat transfer from between the liquid and the gas [W/m2/K]
(
    IOobject
    (
	"lambda_g",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda_g", dimPower/(pow(dimLength, 2)*dimTemperature), 0)
);

volScalarField v_corr // Young's correction factor for the droplet growth model
(
    IOobject
    (
	"v_corr",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("v_corr", dimless, 0)
);

volScalarField etaKantrowitz // Kantrowitz's correction factor on the nucleation rate for non-isothermal effects
(
    IOobject
    (
	"etaKantrowitz",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("etaKantrowitz", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
);

volScalarField J  // Nucleation rate, i.e., number of droplets per unit volume
(
    IOobject
    (
	"J",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("J", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0)
);

volScalarField drdt // Rate of droplet growth in the radial direction
(
    IOobject
    (
	"drdt",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("drdt", dimLength/dimTime, 0)
);

volScalarField beta // Total surface area of droplets, per unit mass of vapour [m2/kg]
(
    IOobject
    (
	"beta",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("beta", dimArea/dimMass, 0)
);

// ------------------------------------- //
// Liquid properties
// ------------------------------------- //

volScalarField surface_tension
(
    IOobject
    (
	"surface_tension",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("surface_tension", dimForce/dimLength, 0)
);

volScalarField dsdT // Rate of change of surface tension with temperature
(
    IOobject
    (
	"dsdT",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dsdT", dimForce/dimLength/dimTemperature, 0)
);

volScalarField rho_liquid // Density of the liquid
(
    IOobject
    (
	"rho_liquid",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_liquid", dimDensity, 0)
);

volScalarField h_liquid // Specific enthalpy energy of the liquid phase
(
    IOobject
    (
	"h_liquid",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("h_liquid", dimEnergy/dimMass, 0)
);

volScalarField h_droplet // Enthalpy of the droplet (Bulk liquid enthalpy + surface energy + kinetic energy)
(
    IOobject
    (
	"h_droplet",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("h_droplet", dimEnergy/dimMass, 0)
);

volScalarField e_liquid // Specific internal energy of the liquid phase [kJ/kg]
(
    IOobject
    (
	"e_liquid",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e_liquid", dimEnergy/dimMass, 0)
);

volScalarField Cp_liquid // Specific heat capacity of the liquid
(
    IOobject
    (
	"Cp_liquid",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Cp_liquid", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0)
);




// ------------------------------------- //
// Source terms for the liquid phase
// ------------------------------------- //

	// Source terms are linearised into passive and active components to promote convergence, according to Picard's rule
	// Active components are, in turn, split into explicit active components and implicit active components -> designated as "active_explicit" and "active_coeff", respectively
	// The passive terms are underrelaxed according to a user-specified factor -> Fields for "old" and "new" values designated below
	// The nucleation term cannot be linearised due to its nature (only has a current-timestep component). It is still underrelaxed, though
		// - Source_Y -> Mass transfer source term
		// - SourceMomentum -> Momentum transfer source term
		// - Source_e -> Energy transfer source term


	// Mass terms

volScalarField Source_Y
(
    IOobject
    (
	"Source_Y",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Source_Y", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
);

volScalarField Source_Y_growth_passive
(
    IOobject
    (
	"Source_Y_growth_passive",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Source_Y_growth_passive", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
);

volScalarField Source_Y_growth_active_coeff 
(
    IOobject
    (
	"Source_Y_growth_active_coeff",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Source_Y_growth_active_coeff", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0)
);

volScalarField Source_Y_growth_active_coeff_new
(
    IOobject
    (
	"Source_Y_growth_active_coeff_new",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Source_Y_growth_active_coeff_new", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0)
);

volScalarField Source_Y_growth_active_explicit 
(
    IOobject
    (
	"Source_Y_growth_active_explicit",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    Source_Y_growth_active_coeff*rho
);

volScalarField Source_Y_nucleation 
(
    IOobject
    (
	"Source_Y_nucleation",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   Source_Y_growth_passive // Same units
);

volScalarField Source_Y_growth 
(
    IOobject
    (
	"Source_Y_growth",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   Source_Y_growth_passive // Same units
);

	// Momentum terms

volVectorField SourceMomentum
(
    IOobject
    (
	"SourceMomentum",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("SourceMomentum", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
);

volScalarField SourceMomentumFactor
(
    IOobject
    (
	"SourceMomentumFactor",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("SourceMomentumFactor", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0)
);

volVectorField SourceMomentum_passive 
(
    IOobject
    (
	"SourceMomentum_passive",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("SourceMomentum_passive", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
);

volVectorField SourceMomentum_active_explicit // Old time step values of U, for the Picard linearisation of the source term
(
    IOobject
    (
	"Source_Y_growth_active_explicit",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
   rho_liquid*drdt*beta*rho*U // Same units
);

volScalarField SourceMomentum_active_coeff // Auxilliary term for the numerical implementation of the source term
(
    IOobject
    (
	"SourceMomentum_active_coeff",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    rho_liquid*drdt*beta // Same units
);


	// Energy terms

volScalarField Source_e
(
    IOobject
    (
	"Source_e",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    e_liquid*Source_Y_growth_passive // Same units
);


// ------------------------------------ //
// Auxilliary dimensionless parameters
// ------------------------------------ //

	// Dimensionless temperature, pressure and density of the gas

volScalarField Tdimless
(
    IOobject
    (
        "Tdimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tdimless", dimless, 0)
);

volScalarField pdimless
(
    IOobject
    (
        "pdimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pdimless", dimless, 0)
);

volScalarField rhodimless
(
    IOobject
    (
        "rhodimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhodimless", dimless, 0)
);


// ------------------------------------ //
// --- Other parameters
// ------------------------------------ //

volScalarField cellVolume  // For computing the volume of each cell in the mesh
(
	IOobject
	(
		"cellVolume",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("cellVolume", dimVolume, 0)
);
