/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free_ software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free_ Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See_ the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see_ <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// For the thermophysical model
#include "VLMEoS.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::VLMEoS<Specie>::VLMEoS
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
)
:
    Specie(sp),
	 ee_(ee),
	 bb_(bb),
	 C_(C),
	 G_(G),
	 m1_(m1),
	 m2_(m2),
	 w1_(w1),
	 w2_(w2),
	 vukal_aux_term3_(vukal_aux_term3)
{}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::VLMEoS<Specie>::VLMEoS
(
    const word& name,
    const VLMEoS& pg
)
:
    Specie(name, pg),
	 ee_(pg.ee_),
	 bb_(pg.bb_),
	 C_(pg.C_),
	 G_(pg.G_),
	 m1_(pg.m1_),
	 m2_(pg.m2_),
	 w1_(pg.w1_),
	 w2_(pg.w2_),
	 vukal_aux_term3_(pg.vukal_aux_term3_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::VLMEoS <Specie>>
Foam::VLMEoS<Specie>::clone() const
{
    return autoPtr<VLMEoS<Specie>>
    (
        new VLMEoS<Specie>(*this)
    );
}

template<class Specie>
inline Foam::autoPtr<Foam::VLMEoS<Specie> >
Foam::VLMEoS<Specie>::New(Istream& is)
{
    return autoPtr<VLMEoS<Specie> >(new VLMEoS<Specie>(is));
}

template<class Specie>
inline Foam::autoPtr<Foam::VLMEoS<Specie>>
Foam::VLMEoS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<VLMEoS<Specie>>
    (
        new VLMEoS<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::B2
(
    scalar T
) const
{
		return (-ee_/(G_*T) - C_*G_/(Foam::pow(T, w1_)) + bb_);
}

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::dB2
(
    scalar T
) const
{
    return (ee_/G_*(1.0/T/T) + C_*G_*w1_*pow(T, -w1_)*(1.0/T));
}

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::rho
(
    scalar p,
    scalar T
) const
{
	if (T >= 240.0)
	{
		return (2.0*p/this->R()/T) / (1.0 + Foam::sqrt(1.0 + 4.0*p*this->B2(T)/this->R()/T));
	}
	else
	{
		return p/(this->R()*T);
	}
}


template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::H(scalar p, scalar T) const
{

	if (T >= 240.0)
	{
		return (p/this->rho(p,T) - this->R()*T*T*this->rho(p,T)*this->dB2(T)
					+ 1000.0*(1.111177*T + 0.000355878*T*T - 6991.96/T + 2070.54)); // 1000.0 factor to bring it from kJ/kg to J/kg
	}
	else
	{
		return (1848.85620*T + 3.7214822*Foam::sinh(0.02059888*T) + 1995458.2016437165);
								// 1995458.2016437165 = - A*Tref_ - Tcrit*B/C*Foam::sinh(C*Tref_/Tcrit) + Href_
	}							// 3.721482183904635 = Tcrit*B/C
								// 0.02059888 = C/Tcrit
}								// From JANAF: Href = 2551013.479 J/kg at Tref = 300.0 K and p = 1bar

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::Cv
(
    const scalar p,
    const scalar T
) const
{
	if (T >= 240.0)
	{
    const scalar rho = this->rho(p,T);
		const scalar d2B2 = -2.0*ee_/G_*(1.0/T/T/T) - C_*G_*w1_*(w1_ + 1.0)*pow(T, -w1_)*(1.0/T/T);
		return
		(
			- this->R()*(2.0*T*(rho*this->dB2(T)) + T*T*(rho*d2B2))
			+ 1000.0*(1.111177 + 0.00071*T + 6991.96/T/T)
		); // 1000 factor to bring it from kJ/kg to J/kg
	}
	else
	{
		return 1848.85620 + 0.076658367*Foam::cosh(0.02059888*T) - this->R();
	}
}

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::Cp
(
    const scalar p,
    const scalar T
) const
{
	if (T >= 240.0)
	{
		const scalar rho = this->rho(p, T);
		const scalar B2 = this->B2(T);
		return(
			this->Cv(p, T) - T*Foam::pow(((this->R() * (rho + rho*rho*B2)) +
				    (this->R() * T * (rho*rho*this->dB2(T)))), 2.0) /
			     (-this->R() * T * (rho*rho + 2.0*rho*rho*rho*B2))
			);
	}
	else
	{
		return this->Cv(p, T) + this->R();
	}
}

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::psi
(
    scalar p,
    scalar T
) const
{
	if (T >= 240.0)
	{
    		return 1.0/(this->Z(p, T)*this->R()*T);
	}
	else
	{
		return 1.0/(this->R()*T);
	}
}


template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::Z
(
    scalar p,
    scalar T
) const
{
	if (T >= 240.0)
	{
    		return (1.0 + this->B2(T)*this->rho(p,T));
	}
	else
	{
		return 1.0;
	}
}


template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
	if (T >= 240.0)
	{
		return this->Cp(p, T) - this->Cv(p, T);
	}
	else
	{
		return this->R();
	}

}

template<class Specie>
inline Foam::scalar Foam::VLMEoS<Specie>::gamma
(
    scalar p,
    scalar T
) const
{
	if (T >= 240.0)
	{
		const scalar rho = this->rho(p, T);

		return (rho/p * this->Cp(p, T)/this->Cv(p, T) *
			this->R()*T*(1.0 + 2.0*rho*this->B2(T)));
	}
	else
	{
		return this->Cp(p, T)/this->Cv(p, T);
	}
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::VLMEoS<Specie>::operator+=
(
    const VLMEoS<Specie>& pg
)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::VLMEoS<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::VLMEoS<Specie> Foam::operator+
(
    const VLMEoS<Specie>& pg1,
    const VLMEoS<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return VLMEoS<Specie>
        (
            sp
        );
    }
    else
    {
        return VLMEoS<Specie>
        (
            sp
        );
    }
}


template<class Specie>
Foam::VLMEoS<Specie> Foam::operator*
(
    const scalar s,
    const VLMEoS<Specie>& pg
)
{
    return VLMEoS<Specie>
    (
        s*static_cast<const Specie&>(pg)
    );
}


template<class Specie>
Foam::VLMEoS<Specie> Foam::operator==
(
    const VLMEoS<Specie>& pg1,
    const VLMEoS<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    return VLMEoS<Specie>
    (
        sp
    );
}


// ************************************************************************* //
