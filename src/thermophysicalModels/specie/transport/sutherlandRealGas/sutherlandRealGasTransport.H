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
    Foam::sutherlandRealGasTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    sutherlandRealGasTransportI.H
    sutherlandRealGasTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandRealGasTransport_H
#define sutherlandRealGasTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class sutherlandRealGasTransport;

template<class Thermo>
inline sutherlandRealGasTransport<Thermo> operator+
(
    const sutherlandRealGasTransport<Thermo>&,
    const sutherlandRealGasTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealGasTransport<Thermo> operator-
(
    const sutherlandRealGasTransport<Thermo>&,
    const sutherlandRealGasTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealGasTransport<Thermo> operator*
(
    const scalar,
    const sutherlandRealGasTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealGasTransport<Thermo> operator==
(
    const sutherlandRealGasTransport<Thermo>&,
    const sutherlandRealGasTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandRealGasTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class sutherlandRealGasTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sutherlandRealGasTransport
:
    public Thermo
{
    // Private data

        // Sutherland's coefficients
        scalar As_, Ts_;


    // Private Member Functions

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
        inline sutherlandRealGasTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts
        );

        //- Construct from two viscosities
        inline sutherlandRealGasTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline sutherlandRealGasTransport(const word&, const sutherlandRealGasTransport&);

        //- Construct from Istream
        sutherlandRealGasTransport(Istream&);

        //- Construct and return a clone
        inline autoPtr<sutherlandRealGasTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<sutherlandRealGasTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<sutherlandRealGasTransport> New(const dictionary& dict);


    // Member functions

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar rho,const scalar T) const;

        //- Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar rho,const scalar T) const;


    // Member operators

        inline sutherlandRealGasTransport& operator=
        (
            const sutherlandRealGasTransport&
        );


    // Friend operators

        friend sutherlandRealGasTransport operator+ <Thermo>
        (
            const sutherlandRealGasTransport&,
            const sutherlandRealGasTransport&
        );

        friend sutherlandRealGasTransport operator- <Thermo>
        (
            const sutherlandRealGasTransport&,
            const sutherlandRealGasTransport&
        );

        friend sutherlandRealGasTransport operator* <Thermo>
        (
            const scalar,
            const sutherlandRealGasTransport&
        );

        friend sutherlandRealGasTransport operator== <Thermo>
        (
            const sutherlandRealGasTransport&,
            const sutherlandRealGasTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sutherlandRealGasTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandRealGasTransportI.H"

#ifdef NoRepository
#   include "sutherlandRealGasTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
