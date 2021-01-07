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

\*---------------------------------------------------------------------------*/

#include "VLMEoS.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::VLMEoS<Specie>::VLMEoS(Istream& is)
:
    Specie(is) 
{
    is.check("VLMEoS<Specie>::VLMEoS(Istream& is)");
}

template<class Specie>
Foam::VLMEoS<Specie>::VLMEoS(const dictionary& dict)
:
    Specie(dict),
    ee_(63.2),
    bb_(0.00085),
    C_(0.3977e+06),
    G_(47.053),
    m1_(1.968),
    m2_(2.957),
    w1_((3.0+2.0*m1_)/2.0),
    w2_((3.0*m2_ - 4.0*m1_)/2.0),
    vukal_aux_term3_(1.0) // The nominal value is assigned below

{
    vukal_aux_term3_ = C_*C_*G_*G_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::VLMEoS<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VLMEoS<Specie>& pg
)
{
    os  << static_cast<const Specie&>(pg);
    os.check("Ostream& operator<<(Ostream& os, const VLMEoS<Specie>& st)");    
    pg.write(os);
    return os;
}


// ************************************************************************* //
