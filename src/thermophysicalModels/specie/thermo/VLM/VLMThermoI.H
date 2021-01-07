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

#include "VLMThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::Cp(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::Cv(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::gamma
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::gamma(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::H(p, T);
}

template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::VLMThermo<EquationOfState>::Hc() const
{
    return 0;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::VLMThermo<EquationOfState>::operator+=
(
    const VLMThermo<EquationOfState>& yt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(yt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = yt.Y()/this->Y();
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::VLMThermo<EquationOfState> Foam::operator+
(
    const VLMThermo<EquationOfState>& yt1,
    const VLMThermo<EquationOfState>& yt2
)
{
    EquationOfState eofs = yt1;
    eofs += yt2;

    if (mag(eofs.Y()) < small)
    {
        return VLMThermo<EquationOfState>
        (
            eofs
        );
    }
    else
    {
        const scalar Y1 = yt1.Y()/eofs.Y();
        const scalar Y2 = yt2.Y()/eofs.Y();

        return VLMThermo<EquationOfState>
        (
            eofs
        );
    }
}

template<class EquationOfState>
inline Foam::VLMThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const VLMThermo<EquationOfState>& yt
)
{
    return VLMThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(yt)
    );
}


template<class EquationOfState>
inline Foam::VLMThermo<EquationOfState> Foam::operator==
(
    const VLMThermo<EquationOfState>& yt1,
    const VLMThermo<EquationOfState>& yt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(yt1)
     == static_cast<const EquationOfState&>(yt2)
    );

    return VLMThermo<EquationOfState>
    (
        eofs
    );
}

// ************************************************************************* //
