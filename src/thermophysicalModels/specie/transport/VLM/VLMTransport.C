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

\*---------------------------------------------------------------------------*/

#include "VLMTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::VLMTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return readScalar(dict.subDict("transport").lookup(coeffName));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::VLMTransport<Thermo>::VLMTransport(const dictionary& dict)
:
    Thermo(dict),
	H0_(1.67752),
	H1_(2.20462),
	H2_(0.6366564),
	H3_(-0.241605),
	muref_(1e-06), // Pa.s
	L0_(2.443221e-03),
	L1_(1.323095e-02),
	L2_(6.770357e-03),
	L3_(-3.454586e-03),
	L4_(4.096266e-04),
	kapparef_(0.001), // W/(m.K)
	Tref_(647.096) // Kelvin
{}


template<class Thermo>
Foam::VLMTransport<Thermo>::VLMTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::VLMTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VLMTransport<Thermo>& kst
)
{
    kst.write(os);
    return os;
}


// ************************************************************************* //
