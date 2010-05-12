/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "specieThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// construct from components
template<class thermo>
inline specieThermo<thermo>::specieThermo
(
    const thermo& sp
)
:
    thermo(sp)
{}


// return the temperature corresponding to the value of the
// thermodynamic property f, given the function f = F(T) and dF(T)/dT
template<class thermo>
inline scalar specieThermo<thermo>::T
(
    scalar f,
    scalar T0,
    scalar (specieThermo<thermo>::*F)(const scalar) const,
    scalar (specieThermo<thermo>::*dFdT)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew = Test - ((this->*F)(Test) - f)/(this->*dFdT)(Test);

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "specieThermo<thermo>::T(scalar f, scalar T0, "
                "scalar (specieThermo<thermo>::*F)(const scalar) const, "
                "scalar (specieThermo<thermo>::*dFdT)(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded.  "
                << "T0 = " << T0 << " Test = " << Test << " Tnew = " << Tnew
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy
template<class thermo>
inline specieThermo<thermo>::specieThermo
(
    const word& name,
    const specieThermo& st
)
:
    thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate and return derived properties
// (These functions need not provided in derived types)

// mole specific properties

// Heat capacities [J/(kmol K)]
template<class thermo>
inline scalar specieThermo<thermo>::cv(const scalar T) const
{
    return this->cp(T) - this->RR;
}

// gamma = cp/cv []
template<class thermo>
inline scalar specieThermo<thermo>::gamma(const scalar T) const
{
    scalar CP = this->cp(T);
    return CP/(CP - this->RR);
}

// Internal energy [J/kmol]
template<class thermo>
inline scalar specieThermo<thermo>::e(const scalar T) const
{
    return this->h(T) - this->RR*(T - this->Tstd);
}

// Gibbs free energy [J/kmol]
template<class thermo>
inline scalar specieThermo<thermo>::g(const scalar T) const
{
    return this->h(T) - T*this->s(T);
}

// Helmholtz free energy [J/kmol]
template<class thermo>
inline scalar specieThermo<thermo>::a(const scalar T) const
{
    return this->e(T) - T*this->s(T);
}


// mass specific properties

// Heat capacity at constant pressure [J/(kg K)]
template<class thermo>
inline scalar specieThermo<thermo>::Cp(const scalar T) const
{
    return this->cp(T)/this->W();
}

// Heat capacity at constant pressure [J/(kg K)]
template<class thermo>
inline scalar specieThermo<thermo>::Cv(const scalar T) const
{
    return this->cv(T)/this->W();
}

// Enthalpy [J/kg]
template<class thermo>
inline scalar specieThermo<thermo>::H(const scalar T) const
{
    return this->h(T)/this->W();
}

// Entropy [J/kg]
template<class thermo>
inline scalar specieThermo<thermo>::S(const scalar T) const
{
    return this->s(T)/this->W();
}

// Internal energy [J/kg]
template<class thermo>
inline scalar specieThermo<thermo>::E(const scalar T) const
{
    return this->e(T)/this->W();
}

// Gibbs free energy [J/kg]
template<class thermo>
inline scalar specieThermo<thermo>::G(const scalar T) const
{
    return this->g(T)/this->W();
}

// Helmholtz free energy [J/kg]
template<class thermo>
inline scalar specieThermo<thermo>::A(const scalar T) const
{
    return this->a(T)/this->W();
}


// Equilibrium reaction thermodynamics

// Equilibrium constant []
template<class thermo>
inline scalar specieThermo<thermo>::K(const scalar T) const
{
    scalar arg = -this->nMoles()*this->g(T)/(this->RR*T);

    if (arg < 600.0)
    {
        return ::exp(arg);
    }
    else
    {
        return VGREAT;
    }
}

//- Equilibrium constant [] i.t.o. partial pressures
//  = PIi(pi/Pstd)^nui
//  For low pressures (where the gas mixture is near perfect) Kp = K
template<class thermo>
inline scalar specieThermo<thermo>::Kp(const scalar T) const
{
    return K(T);
}


//- Equilibrium constant i.t.o. concentration
//  For low pressures (where the gas mixture is near perfect)
//  Kc = Kp(pstd/(RR*T))^nu
template<class thermo>
inline scalar specieThermo<thermo>::Kc(const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/(this->RR*T), this->nMoles());
    }
}


//- Equilibrium constant [] i.t.o. mole-fractions
//  For low pressures (where the gas mixture is near perfect) 
//  Kx = Kp(pstd/p)^nui
template<class thermo>
inline scalar specieThermo<thermo>::Kx(const scalar T, const scalar p) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/p, this->nMoles());
    }
}


//- Equilibrium constant [] i.t.o. number of moles
//  For low pressures (where the gas mixture is near perfect) 
//  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
template<class thermo>
inline scalar specieThermo<thermo>::Kn
(
    const scalar T,
    const scalar p,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(n*this->Pstd/p, this->nMoles());
    }
}


// Energy->temperature  inversion functions

// Temperature from Enthalpy given an initial temperature T0
template<class thermo>
inline scalar specieThermo<thermo>::TH(const scalar h, const scalar T0) const
{
    return T(h, T0, &specieThermo<thermo>::H, &specieThermo<thermo>::Cp);
}


// Temperature from internal energy given an initial temperature T0
template<class thermo>
inline scalar specieThermo<thermo>::TE(const scalar e, const scalar T0) const
{
    return T(e, T0, &specieThermo<thermo>::E, &specieThermo<thermo>::Cv);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void specieThermo<thermo>::operator+=(const specieThermo<thermo>& st)
{
    thermo::operator+=(st);
}

template<class thermo>
inline void specieThermo<thermo>::operator-=(const specieThermo<thermo>& st)
{
    thermo::operator-=(st);
}

template<class thermo>
inline void specieThermo<thermo>::operator*=(const scalar s)
{
    thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline specieThermo<thermo> operator+
(
    const specieThermo<thermo>& st1,
    const specieThermo<thermo>& st2
)
{
    return specieThermo<thermo>
    (
        static_cast<const thermo&>(st1) + static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline specieThermo<thermo> operator-
(
    const specieThermo<thermo>& st1,
    const specieThermo<thermo>& st2
)
{
    return specieThermo<thermo>
    (
        static_cast<const thermo&>(st1) - static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline specieThermo<thermo> operator*
(
    const scalar s,
    const specieThermo<thermo>& st
)
{
    return specieThermo<thermo>
    (
        s*static_cast<const thermo&>(st)
    );
}


template<class thermo>
inline specieThermo<thermo> operator==
(
    const specieThermo<thermo>& st1,
    const specieThermo<thermo>& st2
)
{
    return st2 - st1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
