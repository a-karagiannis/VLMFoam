/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::VLMEoS

Description
    Reduced Vukalovich equation of state for H2O steam.

SourceFiles
    VLMEoSI.H
    VLMEoS.C

\*---------------------------------------------------------------------------*/

#ifndef VLMEoS_H
#define VLMEoS_H

#include "scalarList.H"
#include "autoPtr.H"
#include "scalar.H"
#include "FixedList.H"

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class VLMEoS;

template<class Specie>
inline VLMEoS<Specie> operator+
(
    const VLMEoS<Specie>&,
    const VLMEoS<Specie>&
);

template<class Specie>
inline VLMEoS<Specie> operator*
(
    const scalar,
    const VLMEoS<Specie>&
);

template<class Specie>
inline VLMEoS<Specie> operator==
(
    const VLMEoS<Specie>&,
    const VLMEoS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const VLMEoS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class VLMEoS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class VLMEoS
:
    public Specie
{
    // Private data

	// Empirical constants - All units SI
	scalar ee_;
	scalar bb_; // In the paper called "b", changed here to "bb" to resolve notation conflicts
	scalar C_;
	scalar G_;
	scalar m1_;
	scalar m2_;

	scalar w1_;
	scalar w2_;

	scalar vukal_aux_term3_;

public:

    // Public data

public:

    // Constructors

        //- Construct from components
        inline VLMEoS
	(
	const Specie& sp,
	const scalar& ee,
	const scalar& bb,
	const scalar& C,
	const scalar& G,
	const scalar& m1,
	const scalar& m2,
	const scalar& w1,
	const scalar& w2,
	const scalar& vukal_aux_term3
	);

        //- Construct from Istream
        VLMEoS(Istream&);

        //- Construct from dictionary
        VLMEoS(const dictionary& dict);

        //- Construct as named copy
        inline VLMEoS(const word& name, const VLMEoS&);

        //- Construct and return a clone
        inline autoPtr<VLMEoS> clone() const;

        // Selector from Istream
        inline static autoPtr<VLMEoS> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<VLMEoS> New
        (
            const dictionary& dict
        );

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VLMEoS<" + word(Specie::typeName_()) + '>';
        }

	// Virial coefficient and its (first and second) derivatives

            inline scalar B2(scalar T) const;

            inline scalar dB2(scalar T) const;

            //inline scalar d2B2(scalar T) const;	


        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K)]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return specific heat ratio [-]
            inline scalar gamma(scalar p, scalar T) const;

    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const VLMEoS&);
        inline void operator-=(const VLMEoS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend VLMEoS operator+ <Specie>
        (
            const VLMEoS&,
            const VLMEoS&
        );

        friend VLMEoS operator* <Specie>
        (
            const scalar s,
            const VLMEoS&
        );

        friend VLMEoS operator== <Specie>
        (
            const VLMEoS&,
            const VLMEoS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const VLMEoS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VLMEoSI.H"

#ifdef NoRepository
    #include "VLMEoS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
