/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Class
    Foam::sutherlandTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    @f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    @f]

SourceFiles
    sutherlandTransportI.H
    sutherlandTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandTransport_H
#define sutherlandTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class sutherlandTransport;

template<class thermo>
inline sutherlandTransport<thermo> operator+
(
    const sutherlandTransport<thermo>&,
    const sutherlandTransport<thermo>&
);

template<class thermo>
inline sutherlandTransport<thermo> operator-
(
    const sutherlandTransport<thermo>&,
    const sutherlandTransport<thermo>&
);

template<class thermo>
inline sutherlandTransport<thermo> operator*
(
    const scalar,
    const sutherlandTransport<thermo>&
);

template<class thermo>
inline sutherlandTransport<thermo> operator==
(
    const sutherlandTransport<thermo>&,
    const sutherlandTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class sutherlandTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class sutherlandTransport
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
        inline sutherlandTransport
        (
            const thermo& t,
            const scalar as,
            const scalar ts
        );

        //- Construct from two viscosities
        inline sutherlandTransport
        (
            const thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline sutherlandTransport(const word&, const sutherlandTransport&);

        //- Construct from Istream
        sutherlandTransport(Istream&);

        //- Construct and return a clone
        inline autoPtr<sutherlandTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<sutherlandTransport> New(Istream& is);


    // Member functions

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar T) const;

        //- Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar T) const;


    // Member operators

        inline sutherlandTransport& operator=
        (
            const sutherlandTransport&
        );


    // Friend operators

        friend sutherlandTransport operator+ <thermo>
        (
            const sutherlandTransport&,
            const sutherlandTransport&
        );

        friend sutherlandTransport operator- <thermo>
        (
            const sutherlandTransport&,
            const sutherlandTransport&
        );

        friend sutherlandTransport operator* <thermo>
        (
            const scalar,
            const sutherlandTransport&
        );

        friend sutherlandTransport operator== <thermo>
        (
            const sutherlandTransport&,
            const sutherlandTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const sutherlandTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandTransportI.H"

#ifdef NoRepository
#   include "sutherlandTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
