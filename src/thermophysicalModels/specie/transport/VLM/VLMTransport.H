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
    Foam::VLMTransport

Description
    Transport package using a truncated version of the standard IAPWS formula for the viscosity of steam.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    VLMTransportI.H
    VLMTransport.C

\*---------------------------------------------------------------------------*/

#ifndef VLMTransport_H
#define VLMTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class VLMTransport;

template<class Thermo>
inline VLMTransport<Thermo> operator+
(
    const VLMTransport<Thermo>&,
    const VLMTransport<Thermo>&
);

template<class Thermo>
inline VLMTransport<Thermo> operator*
(
    const scalar,
    const VLMTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const VLMTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class VLMTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class VLMTransport
:
    public Thermo
{
    // Private data

	// Coefficients for the mu and kappa formulas
	scalar H0_;
	scalar H1_;
	scalar H2_;
	scalar H3_;
	scalar muref_;
	scalar L0_;
	scalar L1_;
	scalar L2_;
	scalar L3_;
	scalar L4_;
	scalar kapparef_;
	scalar Tref_;

    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline VLMTransport
        (
            const Thermo& t,
		const scalar& H0,
		const scalar& H1,
		const scalar& H2,
		const scalar& H3,
		const scalar& muref,
		const scalar& L0,
		const scalar& L1,
		const scalar& L2,
		const scalar& L3,
		const scalar& L4,
		const scalar& kapparef,
		const scalar& Tref
        );

        //- Construct as named copy
        inline VLMTransport(const word&, const VLMTransport&);

        //- Construct from dictionary
        VLMTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        VLMTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<VLMTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<VLMTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VLMTransport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const VLMTransport&);

        inline void operator+=(const VLMTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend VLMTransport operator+ <Thermo>
        (
            const VLMTransport&,
            const VLMTransport&
        );

        friend VLMTransport operator* <Thermo>
        (
            const scalar,
            const VLMTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const VLMTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VLMTransportI.H"

#ifdef NoRepository
    #include "VLMTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
