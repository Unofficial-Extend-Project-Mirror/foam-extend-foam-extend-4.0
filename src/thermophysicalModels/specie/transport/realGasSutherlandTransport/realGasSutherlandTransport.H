/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

Class
    Foam::realGasSutherlandTransport

Description

    Copy of the sutherlandTransport class, modified to be used with real gas properties

    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    @f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    @f]

SourceFiles
    realGasSutherlandTransportI.H
    realGasSutherlandTransport.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef realGasSutherlandTransport_H
#define realGasSutherlandTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class realGasSutherlandTransport;

template<class thermo>
inline realGasSutherlandTransport<thermo> operator+
(
    const realGasSutherlandTransport<thermo>&,
    const realGasSutherlandTransport<thermo>&
);

template<class thermo>
inline realGasSutherlandTransport<thermo> operator-
(
    const realGasSutherlandTransport<thermo>&,
    const realGasSutherlandTransport<thermo>&
);

template<class thermo>
inline realGasSutherlandTransport<thermo> operator*
(
    const scalar,
    const realGasSutherlandTransport<thermo>&
);

template<class thermo>
inline realGasSutherlandTransport<thermo> operator==
(
    const realGasSutherlandTransport<thermo>&,
    const realGasSutherlandTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const realGasSutherlandTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class realGasSutherlandTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class realGasSutherlandTransport
:
    public thermo
{
    // Private data

        // Sutherland's coefficients
        scalar As, Ts;


    // Private member functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline realGasSutherlandTransport
        (
            const thermo& t,
            const scalar as,
            const scalar ts
        );

        //- Construct from two viscosities
        inline realGasSutherlandTransport
        (
            const thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline realGasSutherlandTransport(const word&, const realGasSutherlandTransport&);

        //- Construct from Istream
        realGasSutherlandTransport(Istream&);

        //- Construct and return a clone
        inline autoPtr<realGasSutherlandTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<realGasSutherlandTransport> New(Istream& is);


    // Member functions 


        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;


        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar rho,const scalar T) const;

        //- Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar rho,const scalar T) const;



    // Member operators

        inline realGasSutherlandTransport& operator=
        (
            const realGasSutherlandTransport&
        );


    // Friend operators

        friend realGasSutherlandTransport operator+ <thermo>
        (
            const realGasSutherlandTransport&,
            const realGasSutherlandTransport&
        );

        friend realGasSutherlandTransport operator- <thermo>
        (
            const realGasSutherlandTransport&,
            const realGasSutherlandTransport&
        );

        friend realGasSutherlandTransport operator* <thermo>
        (
            const scalar,
            const realGasSutherlandTransport&
        );

        friend realGasSutherlandTransport operator== <thermo>
        (
            const realGasSutherlandTransport&,
            const realGasSutherlandTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const realGasSutherlandTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "realGasSutherlandTransportI.H"

#ifdef NoRepository
#   include "realGasSutherlandTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
