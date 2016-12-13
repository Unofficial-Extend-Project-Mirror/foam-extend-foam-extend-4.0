/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#include "realGasSpecieThermo.H"

template<class thermo>
inline Foam::realGasSpecieThermo<thermo>::realGasSpecieThermo
(
    const thermo& sp
)
:
    thermo(sp)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


//CL: using two one dimensional newton solvers in a row
template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::T
(
    scalar f,
    scalar &T0,
    scalar p,
    scalar &rho0,
    scalar (realGasSpecieThermo<thermo>::*F)(const scalar,const scalar) const,
    scalar (realGasSpecieThermo<thermo>::*dFdT)(const scalar,const scalar) const
) const
{
    scalar Test ;
    scalar Tnew = T0;
    scalar rhoOld;
    scalar rho = rho0;
    scalar Ttol = T0*tol_();
    scalar rhotol = rho0*tol_();
    label iter = 0;
    label i;

    do
    {
        Test = Tnew;
        rhoOld = rho;
        rho = this->rho(p, Test, rhoOld);
        i = 0;

        do
        {
            //CL: using a stabilizing newton solver
            //CL: if the solve is diverging, the step is reduced until the solver converges
            Tnew = Test - ((this->*F)(rho, Test) - f)/(this->*dFdT)(rho, Test)/(pow(2, i));
            i++;
        } while
        (
            (i<20)
            &&
            (
                mag((this->*F)(rho, Tnew) - f)
                >
                mag((this->*F)(rho, Test) - f)
            )
        );

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "realGasSpecieThermo<thermo>::T(scalar f, scalar T0, scalar p, scalar rho0, "
                "scalar (realGasSpecieThermo<thermo>::*F)(const scalar) const, "
                "scalar (realGasSpecieThermo<thermo>::*dFdT)(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }
    } while
    //CL: both fields must converge
    (
        (mag(mag(Tnew) - mag(Test)) > Ttol)
        ||
        (mag(mag(rho) - mag(rhoOld)) > rhotol)
    );

    rho0=rho;
    T0=Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::realGasSpecieThermo<thermo>::realGasSpecieThermo
(
    const word& name,
    const realGasSpecieThermo& st
)
:
    thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::gamma(const scalar rho,  const scalar T ) const
{
    return -1/(rho*this->p(rho, T))*this->cp(rho, T)/this->cv(rho, T)*this->dpdv(rho, T);
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::g(const scalar rho,  const scalar T ) const
{
    return this->h(rho, T) - T*this->s(rho, T);
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::a(const scalar rho, const scalar T ) const
{
    return this->e(rho,T ) - T*this->s(rho, T);
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::Cp( const scalar rho, const scalar T) const
{
    return this->cp(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::Cv( const scalar rho, const scalar T) const
{
    return this->cv(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::H(const scalar rho, const scalar T) const
{
    return this->h(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::S(const scalar rho, const scalar T) const
{
    return this->s(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::E(const scalar rho, const scalar T) const
{
    return this->e(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::G(const scalar rho,  const scalar T) const
{
    return this->g(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::A(const scalar rho,  const scalar T) const
{
    return this->a(rho, T)/this->W();
}

//CL:- Return compressibility drho/dp at h=constant [s^2/m^2]
//CL:- using Bridgeman's Table
template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::psiH
(
    const scalar rho,
    const scalar T

) const
{
    scalar beta=this->isobarExpCoef(rho, T);

    return
    -(
        (T*beta*beta-beta)/this->Cp(rho, T)
      - this->isothermalCompressiblity(rho, T)*rho
    );
}

//CL:- Return compressibility drho/dp at e=constant [s^2/m^2]
//CL:- using Bridgeman's Table
template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::psiE
(
    const scalar rho,
    const scalar T

) const
{
    scalar V = 1/rho;
    scalar cp = this->Cp(rho, T);
    scalar beta = this->isobarExpCoef(rho, T);

    return
    -(
        (
            T*beta*beta*V - this->isothermalCompressiblity(rho, T)*cp
        )
        /
        (
            cp*V - beta*this->p(rho, T)*V*V
        )
    );
}

//CL:- Returns drho/dH at p=constant
//CL:- using Bridgeman's Table
template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::drhodH
(
    const scalar rho,
    const scalar T

) const
{
    return -(rho*this->isobarExpCoef(rho,T))/this->Cp(rho,T);
}

//CL:- Returns drho/dE at p=constant
//CL:- using Bridgeman's Table
template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::drhodE
(
    const scalar rho,
    const scalar T

) const
{
    scalar beta=this->isobarExpCoef(rho,T);

    return -(rho*beta)/(this->Cp(rho,T)-beta*this->p(rho,T)/rho);
}

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::TH
(
    const scalar h,
    scalar &T0,
    const scalar p,
    scalar &rho0
) const
{
   T(h, T0, p,rho0, &realGasSpecieThermo<thermo>::H, &realGasSpecieThermo<thermo>::Cp);
}


template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::TE
(
    const scalar e,
    scalar &T0,
    const scalar p,
    scalar &rho0
)  const
{
    T(e, T0, p,rho0, &realGasSpecieThermo<thermo>::E, &realGasSpecieThermo<thermo>::Cv);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator+=
(
    const realGasSpecieThermo<thermo>& st
)
{
    thermo::operator+=(st);
}

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator-=
(
    const realGasSpecieThermo<thermo>& st
)
{
    thermo::operator-=(st);
}

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator*=(const scalar s)
{
    thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator+
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return realGasSpecieThermo<thermo>
    (
        static_cast<const thermo&>(st1) + static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator-
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return realGasSpecieThermo<thermo>
    (
        static_cast<const thermo&>(st1) - static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator*
(
    const scalar s,
    const realGasSpecieThermo<thermo>& st
)
{
    return realGasSpecieThermo<thermo>
    (
        s*static_cast<const thermo&>(st)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator==
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
