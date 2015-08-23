/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandRealGasTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo>::sutherlandRealGasTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts)
{}


template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo>::sutherlandRealGasTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo>::sutherlandRealGasTransport
(
    const word& name,
    const sutherlandRealGasTransport& st
)
:
    Thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealGasTransport<Thermo> >
Foam::sutherlandRealGasTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandRealGasTransport<Thermo> >
    (
        new sutherlandRealGasTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealGasTransport<Thermo> >
Foam::sutherlandRealGasTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<sutherlandRealGasTransport<Thermo> >
    (
        new sutherlandRealGasTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealGasTransport<Thermo> >
Foam::sutherlandRealGasTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandRealGasTransport<Thermo> >
    (
        new sutherlandRealGasTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandRealGasTransport<Thermo>::mu(const scalar T) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


// CL: for real gas thermo
// Thermal conductivity [W/mK]
template<class thermo>
inline Foam::scalar Foam::sutherlandRealGasTransport<thermo>::kappa
(
    const scalar rho,
    const scalar T
) const
{
    scalar Cv_ = this->Cv(rho,T);
    return mu(T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}


// CL: for real gas thermo
// Thermal diffusivity for enthalpy [kg/ms]
template<class thermo>
inline Foam::scalar Foam::sutherlandRealGasTransport<thermo>::alpha
(
    const scalar rho,
    const scalar T
) const
{
    scalar Cv_ = this->Cv(rho,T);
    scalar Cp_ = this->Cp(rho,T);

    scalar deltaT = T - specie::Tstd;
    scalar CpBar =
        (deltaT*(this->H(rho, T) - this->H(this->rhostd(),specie::Tstd)) + Cp_)/(sqr(deltaT) + 1);

    return mu(T)*Cv_*(1.32 + 1.77*this->R()/Cv_)/CpBar;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo>&
Foam::sutherlandRealGasTransport<Thermo>::operator=
(
    const sutherlandRealGasTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    As_ = st.As_;
    Ts_ = st.Ts_;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo> Foam::operator+
(
    const sutherlandRealGasTransport<Thermo>& st1,
    const sutherlandRealGasTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandRealGasTransport<Thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo> Foam::operator-
(
    const sutherlandRealGasTransport<Thermo>& st1,
    const sutherlandRealGasTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandRealGasTransport<Thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandRealGasTransport<Thermo>& st
)
{
    return sutherlandRealGasTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_
    );
}


template<class Thermo>
inline Foam::sutherlandRealGasTransport<Thermo> Foam::operator==
(
    const sutherlandRealGasTransport<Thermo>& st1,
    const sutherlandRealGasTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
