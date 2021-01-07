// Set a minimum limit to the energy
// to prevent numerical instabilities 
// at the early stages of the simulation
    dimensionedScalar Cp_bound("Cp_bound", dimEnergy/dimMass/dimTemperature, 1850.0); // J/kg/K
    dimensionedScalar Tref("Tref", dimTemperature, 300.0);  // K
    dimensionedScalar Tbound("Tbound", dimTemperature, 40.0); // K
    dimensionedScalar pbound("pbound", dimPressure, 1.0); // Pa
    dimensionedScalar Href("Href", dimEnergy/dimMass, 2551013.479); // J/kg
    dimensionedScalar rhobound("rhobound", dimMass/dimVolume, 5.41688e-05); // pbound/461.52/Tbound
    dimensionedScalar Emin("Emin", dimEnergy/dimMass, 2051552.668); // Href + Cp_bound*(Tbound - Tref) - pbound/rhobound
