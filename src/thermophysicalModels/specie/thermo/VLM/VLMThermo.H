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

Class
    Foam::VLMThermo

Description
    Thermodynamics package for the properties of metastable H2O steam, according to 
	"Gerber and Kermani (2004) - A pressure based Eulerian–Eulerian multi-phase model for non-
	equilibrium condensation in transonic steam flow, vol. 47, nr. 10-11, pp. 2217-2231"

SourceFiles
    VLMThermoI.H
    VLMThermo.C

\*---------------------------------------------------------------------------*/

#ifndef VLMThermo_H
#define VLMThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class VLMThermo;

template<class EquationOfState>
inline VLMThermo<EquationOfState> operator+
(
    const VLMThermo<EquationOfState>&,
    const VLMThermo<EquationOfState>&
);

template<class EquationOfState>
inline VLMThermo<EquationOfState> operator*
(
    const scalar,
    const VLMThermo<EquationOfState>&
);

template<class EquationOfState>
inline VLMThermo<EquationOfState> operator==
(
    const VLMThermo<EquationOfState>&,
    const VLMThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const VLMThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class VLMThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class VLMThermo
:
    public EquationOfState
{

public:

    // Public data

private:

    // Private data

    // Private Member Functions

public:

    // Constructors

        //- Construct from components
        inline VLMThermo
        (
            const EquationOfState& st
        );

        //- Construct from dictionary
        VLMThermo(const dictionary& dict);

        //- Construct as a named copy
        inline VLMThermo(const word&, const VLMThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VLM<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;

        // Access

        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume [J/(kg K)]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Ratio of specific heats (real-gas version) [-]
            inline scalar gamma(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            //inline scalar S(const scalar p, const scalar T) const;

        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const VLMThermo&);


    // Friend operators

        friend VLMThermo operator+ <EquationOfState>
        (
            const VLMThermo&,
            const VLMThermo&
        );

        friend VLMThermo operator* <EquationOfState>
        (
            const scalar,
            const VLMThermo&
        );

        friend VLMThermo operator== <EquationOfState>
        (
            const VLMThermo&,
            const VLMThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const VLMThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VLMThermoI.H"

#ifdef NoRepository
    #include "VLMThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
