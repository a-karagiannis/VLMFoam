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

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::VLMTransport<Thermo>::VLMTransport
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
)
:
    Thermo(t),
	H0_(H0),
	H1_(H1),
	H2_(H2),
	H3_(H3),
	muref_(muref),
	L0_(L0),
	L1_(L1),
	L2_(L2),
	L3_(L3),
	L4_(L4),
	kapparef_(kapparef),
	Tref_(Tref)
{}


template<class Thermo>
inline Foam::VLMTransport<Thermo>::VLMTransport
(
    const word& name,
    const VLMTransport& kst
)
:
    Thermo(name, kst),
	H0_(kst.H0_),
	H1_(kst.H1_),
	H2_(kst.H2_),
	H3_(kst.H3_),
	muref_(kst.muref_),
	L0_(kst.L0_),
	L1_(kst.L1_),
	L2_(kst.L2_),
	L3_(kst.L3_),
	L4_(kst.L4_),
	kapparef_(kst.kapparef_),
	Tref_(kst.Tref_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::VLMTransport<Thermo>>
Foam::VLMTransport<Thermo>::clone() const
{
    return autoPtr<VLMTransport<Thermo>>
    (
        new VLMTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::VLMTransport<Thermo>>
Foam::VLMTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<VLMTransport<Thermo>>
    (
        new VLMTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::VLMTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
	if (T >= 273.16)
	{
	    const scalar Tbar = T/Tref_;
	    return muref_*(100.0*Foam::sqrt(Tbar))/(H0_ + H1_/Tbar + H2_/Tbar/Tbar + H3_/Tbar/Tbar/Tbar);
	}
	else
	{
		return (3.51068999e-06
			+ 3.27641137e-09 * T
			+ 7.63396116e-11 * T*T
			- 6.52022078e-14 * T*T*T
			+ 1.96761568e-17 * T*T*T*T
			);
	}
}


template<class Thermo>
inline Foam::scalar Foam::VLMTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
	if (T > 273.16)
	{
		const scalar Tbar = T/Tref_;

		return kapparef_*Foam::sqrt(Tbar)/(L0_ + L1_/Tbar + L2_/Tbar/Tbar + L3_/Tbar/Tbar/Tbar + L4_/Tbar/Tbar/Tbar/Tbar);
	}
	else
	{
		return (7.24966752e-03
			- 1.58112357e-06 * T
			+ 1.54438048e-07 * T*T
			- 8.44198169e-11 * T*T*T
			+ 2.01111364e-14 * T*T*T*T
			);
	}
}


template<class Thermo>
inline Foam::scalar Foam::VLMTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/(this->Cp(p, T));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::VLMTransport<Thermo>::operator=
(
    const VLMTransport<Thermo>& kst
)
{
    Thermo::operator=(kst);
}


template<class Thermo>
inline void Foam::VLMTransport<Thermo>::operator+=
(
    const VLMTransport<Thermo>& kst
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(kst);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = kst.Y()/this->Y();
    }
}


template<class Thermo>
inline void Foam::VLMTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::VLMTransport<Thermo> Foam::operator+
(
    const VLMTransport<Thermo>& kst1,
    const VLMTransport<Thermo>& kst2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(kst1) + static_cast<const Thermo&>(kst2)
    );

    if (mag(t.Y()) < small)
    {
        return VLMTransport<Thermo>
        (
            t,
            0
        );
    }
    else
    {
        scalar Y1 = kst1.Y()/t.Y();
        scalar Y2 = kst2.Y()/t.Y();

        return VLMTransport<Thermo>
        (
            t
        );
    }
}


template<class Thermo>
inline Foam::VLMTransport<Thermo> Foam::operator*
(
    const scalar s,
    const VLMTransport<Thermo>& kst
)
{
    return VLMTransport<Thermo>
    (
        s*static_cast<const Thermo&>(kst)
    );
}


// ************************************************************************* //
