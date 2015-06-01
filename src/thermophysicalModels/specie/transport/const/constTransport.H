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
    Foam::constTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constTransportI.H
    constTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constTransport_H
#define constTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class constTransport;

template<class thermo>
inline constTransport<thermo> operator+
(
    const constTransport<thermo>&,
    const constTransport<thermo>&
);

template<class thermo>
inline constTransport<thermo> operator-
(
    const constTransport<thermo>&,
    const constTransport<thermo>&
);

template<class thermo>
inline constTransport<thermo> operator*
(
    const scalar,
    const constTransport<thermo>&
);

template<class thermo>
inline constTransport<thermo> operator==
(
    const constTransport<thermo>&,
    const constTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const constTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class constTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class constTransport
:
    public thermo
{
    // Private data

        //- Constant viscosity and reciprocal Prandtl Number.
        scalar Mu, rPr;


    // Private member functions

        //- Construct from components
        inline constTransport
        (
            const thermo& t,
            const scalar nu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline constTransport(const word&, const constTransport&);

        //- Construct from Istream
        constTransport(Istream&);

        //- Construct and return a clone
        inline autoPtr<constTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<constTransport> New(Istream& is);


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

        inline constTransport& operator=
        (
            const constTransport&
        );


    // Friend operators

        friend constTransport operator+ <thermo>
        (
            const constTransport&,
            const constTransport&
        );

        friend constTransport operator- <thermo>
        (
            const constTransport&,
            const constTransport&
        );

        friend constTransport operator* <thermo>
        (
            const scalar,
            const constTransport&
        );

        friend constTransport operator== <thermo>
        (
            const constTransport&,
            const constTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const constTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constTransportI.H"

#ifdef NoRepository
#   include "constTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
