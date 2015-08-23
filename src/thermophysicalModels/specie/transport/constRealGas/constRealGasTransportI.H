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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constRealGasTransport<Thermo>::constRealGasTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::constRealGasTransport<Thermo>::constRealGasTransport
(
    const word& name,
    const constRealGasTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constRealGasTransport<Thermo>::mu(const scalar) const
{
    return mu_;
}


// CL: for real gas thermo
// Thermal conductivity [W/mK]
template<class thermo>
inline Foam::scalar Foam::constRealGasTransport<thermo>::kappa(const scalar rho,const scalar T) const
{
    return this->Cp(rho,T)*mu(T)*rPr_;
}


// CL: for real gas thermo
// Thermal diffusivity for enthalpy [kg/ms]
template<class thermo>
inline Foam::scalar Foam::constRealGasTransport<thermo>::alpha(const scalar rho,const scalar T) const
{
    scalar Cp_ = this->Cp(rho,T);

    scalar deltaT = T - specie::Tstd;
    scalar CpBar =
        (deltaT*(this->H(rho,T) - this->H(this->rhostd(),specie::Tstd)) + Cp_)/(sqr(deltaT) + 1);

    return Cp_*mu(T)*rPr_/CpBar;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constRealGasTransport<Thermo>& Foam::constRealGasTransport<Thermo>::operator=
(
    const constRealGasTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constRealGasTransport<Thermo> Foam::operator+
(
    const constRealGasTransport<Thermo>& ct1,
    const constRealGasTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constRealGasTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ + molr2*ct2.mu_,
        molr1*ct1.rPr_ + molr2*ct2.rPr_
    );
}


template<class Thermo>
inline Foam::constRealGasTransport<Thermo> Foam::operator-
(
    const constRealGasTransport<Thermo>& ct1,
    const constRealGasTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constRealGasTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ - molr2*ct2.mu_,
        molr1*ct1.rPr_ - molr2*ct2.rPr_
    );
}


template<class Thermo>
inline Foam::constRealGasTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constRealGasTransport<Thermo>& ct
)
{
    return constRealGasTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.rPr_
    );
}


template<class Thermo>
inline Foam::constRealGasTransport<Thermo> Foam::operator==
(
    const constRealGasTransport<Thermo>& ct1,
    const constRealGasTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
