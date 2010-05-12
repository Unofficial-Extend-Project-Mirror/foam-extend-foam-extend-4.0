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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from components
template<class equationOfState>
inline hConstThermo<equationOfState>::hConstThermo
(
    const equationOfState& st,
    const scalar cp,
    const scalar hf
)
:
    equationOfState(st),
    CP(cp),
    Hf(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy
template<class equationOfState>
inline hConstThermo<equationOfState>::hConstThermo
(
    const word& name,
    const hConstThermo& ct
)
:
    equationOfState(name, ct),
    CP(ct.CP),
    Hf(ct.Hf)
{}


// Construct and return a clone
template<class equationOfState>
inline autoPtr<hConstThermo<equationOfState> > hConstThermo<equationOfState>::
clone() const
{
    return autoPtr<hConstThermo<equationOfState> >
    (
        new hConstThermo<equationOfState>(*this)
    );
}


// Selector from Istream
template<class equationOfState>
inline autoPtr<hConstThermo<equationOfState> > hConstThermo<equationOfState>::
New(Istream& is)
{
    return autoPtr<hConstThermo<equationOfState> >
    (
        new hConstThermo<equationOfState>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Heat capacity at constant pressure [J/(kmol K)]
template<class equationOfState>
inline scalar hConstThermo<equationOfState>::cp(const scalar) const
{
    return CP*this->W();
}


//- Enthalpy [J/kmol]
template<class equationOfState>
inline scalar hConstThermo<equationOfState>::h(const scalar T) const
{
    return (CP*T + Hf)*this->W();
}


//- Entropy [J/(kmol K)]
template<class equationOfState>
inline scalar hConstThermo<equationOfState>::s(const scalar T) const
{
    notImplemented("scalar hConstThermo<equationOfState>::s(const scalar T) const");
    return T;
}

/*
//- Temperature from Enthalpy given an initial temperature T0
template<class equationOfState>
inline scalar hConstThermo<equationOfState>::TH(const scalar h, const scalar T0) const
{
    return (h - Hf)/Cp(T0);
}


//- Temperature from internal energy given an initial temperature T0
template<class equationOfState>
inline scalar hConstThermo<equationOfState>::TE(const scalar e, const scalar T0) const
{
    return (e - Hf)/Cv(T0);
}
*/

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline hConstThermo<equationOfState>& hConstThermo<equationOfState>::operator=
(
    const hConstThermo& ct
)
{
    equationOfState::operator=(ct);

    CP = ct.CP;
    Hf = ct.Hf;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline hConstThermo<equationOfState> operator+
(
    const hConstThermo<equationOfState>& ct1,
    const hConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
        + static_cast<const equationOfState&>(ct2)
    );

    return hConstThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.CP + ct2.nMoles()/eofs.nMoles()*ct2.CP,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf + ct2.nMoles()/eofs.nMoles()*ct2.Hf
    );
}


template<class equationOfState>
inline hConstThermo<equationOfState> operator-
(
    const hConstThermo<equationOfState>& ct1,
    const hConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      - static_cast<const equationOfState&>(ct2)
    );

    return hConstThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.CP - ct2.nMoles()/eofs.nMoles()*ct2.CP,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf - ct2.nMoles()/eofs.nMoles()*ct2.Hf
    );
}


template<class equationOfState>
inline hConstThermo<equationOfState> operator*
(
    const scalar s,
    const hConstThermo<equationOfState>& ct
)
{
    return hConstThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.CP,
        ct.Hf
    );
}


template<class equationOfState>
inline hConstThermo<equationOfState> operator==
(
    const hConstThermo<equationOfState>& ct1,
    const hConstThermo<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
