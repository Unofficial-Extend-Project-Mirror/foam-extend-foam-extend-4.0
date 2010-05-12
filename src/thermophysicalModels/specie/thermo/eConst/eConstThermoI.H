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
inline eConstThermo::eConstThermo
(
    const specieThermo& st,
    const scalar cv,
    const scalar hf
)
:
    specieThermo(st),
    CV(cv),
    Hf(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct as named copy
inline eConstThermo::eConstThermo(const word& name, const eConstThermo& ct)
:
    specieThermo(name, ct),
    CV(ct.CV),
    Hf(ct.Hf)
{}


// Construct and return a clone
template<class equationOfState>
inline autoPtr<eConstThermo<equationOfState> > eConstThermo<equationOfState>::
clone() const
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(*this)
    );
}


// Selector from Istream
template<class equationOfState>
inline autoPtr<eConstThermo<equationOfState> > eConstThermo<equationOfState>::
New(Istream& is)
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Heat capacity at constant pressure [J/(kmol K)]
inline scalar eConstThermo::cp(const scalar) const
{
    return CV*W() + RR;
}


//- Enthalpy [J/kmol]
inline scalar eConstThermo::h(const scalar T) const
{
    return cp(T)*T + Hf*W();
}


//- Entropy [J/(kmol K)]
inline scalar eConstThermo::s(const scalar T) const
{
    notImplemented("scalar eConstThermo::s(const scalar T) const");
    return T;
}


//- Temperature from Enthalpy given an initial temperature T0
inline scalar eConstThermo::TH(const scalar h, const scalar T0) const
{
    return (h - Hf)/Cp(T0);
}


//- Temperature from internal energy given an initial temperature T0
inline scalar eConstThermo::TE(const scalar e, const scalar) const
{
    return (e - Hf)/CV;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline  eConstThermo& eConstThermo::operator=
(
    const eConstThermo& ct
)
{
    specieThermo::operator=(ct);

    CV = ct.CV;
    Hf = ct.Hf;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline eConstThermo operator+
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    specieThermo st(((const specieThermo&)ct1) + ((const specieThermo&)ct2));

    return eConstThermo
    (
        st,
        ct1.nMoles()/st.nMoles()*ct1.CV + ct2.nMoles()/st.nMoles()*ct2.CV,
        ct1.nMoles()/st.nMoles()*ct1.Hf + ct2.nMoles()/st.nMoles()*ct2.Hf
    );
}


inline eConstThermo operator-
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    specieThermo st(((const specieThermo&)ct1) - ((const specieThermo&)ct2));

    return eConstThermo
    (
        st,
        ct1.nMoles()/st.nMoles()*ct1.CV - ct2.nMoles()/st.nMoles()*ct2.CV,
        ct1.nMoles()/st.nMoles()*ct1.Hf - ct2.nMoles()/st.nMoles()*ct2.Hf
    );
}


inline eConstThermo operator*
(
    const scalar s,
    const eConstThermo& ct
)
{
    return eConstThermo
    (
        s*((const specieThermo&)ct),
        ct.CV,
        ct.Hf
    );
}


inline eConstThermo operator==
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    return ct2 - ct1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
