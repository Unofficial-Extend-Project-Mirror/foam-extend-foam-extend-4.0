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

Class
    Foam::constRealGasTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constRealGasTransportI.H
    constRealGasTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constRealGasTransport_H
#define constRealGasTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class constRealGasTransport;

template<class Thermo>
inline constRealGasTransport<Thermo> operator+
(
    const constRealGasTransport<Thermo>&,
    const constRealGasTransport<Thermo>&
);

template<class Thermo>
inline constRealGasTransport<Thermo> operator-
(
    const constRealGasTransport<Thermo>&,
    const constRealGasTransport<Thermo>&
);

template<class Thermo>
inline constRealGasTransport<Thermo> operator*
(
    const scalar,
    const constRealGasTransport<Thermo>&
);

template<class Thermo>
inline constRealGasTransport<Thermo> operator==
(
    const constRealGasTransport<Thermo>&,
    const constRealGasTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constRealGasTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class constRealGasTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constRealGasTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline constRealGasTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline constRealGasTransport(const word&, const constRealGasTransport&);

        //- Construct from Istream
        constRealGasTransport(Istream&);


    // Member functions

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar rho, const scalar T) const;

        //- Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar rho, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar T) const;


    // Member operators

        inline constRealGasTransport& operator=
        (
            const constRealGasTransport&
        );


    // Friend operators

        friend constRealGasTransport operator+ <Thermo>
        (
            const constRealGasTransport&,
            const constRealGasTransport&
        );

        friend constRealGasTransport operator- <Thermo>
        (
            const constRealGasTransport&,
            const constRealGasTransport&
        );

        friend constRealGasTransport operator* <Thermo>
        (
            const scalar,
            const constRealGasTransport&
        );

        friend constRealGasTransport operator== <Thermo>
        (
            const constRealGasTransport&,
            const constRealGasTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constRealGasTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constRealGasTransportI.H"

#ifdef NoRepository
#   include "constRealGasTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
