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

Class
    Foam::realGasConstTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    realGasConstTransportI.H
    realGasConstTransport.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany
\*---------------------------------------------------------------------------*/

#ifndef realGasConstTransport_H
#define realGasConstTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class realGasConstTransport;

template<class thermo>
inline realGasConstTransport<thermo> operator+
(
    const realGasConstTransport<thermo>&,
    const realGasConstTransport<thermo>&
);

template<class thermo>
inline realGasConstTransport<thermo> operator-
(
    const realGasConstTransport<thermo>&,
    const realGasConstTransport<thermo>&
);

template<class thermo>
inline realGasConstTransport<thermo> operator*
(
    const scalar,
    const realGasConstTransport<thermo>&
);

template<class thermo>
inline realGasConstTransport<thermo> operator==
(
    const realGasConstTransport<thermo>&,
    const realGasConstTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const realGasConstTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class realGasConstTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class realGasConstTransport
:
    public thermo
{
    // Private data

        //- Constant viscosity and reciprocal Prandtl Number.
        scalar Mu, rPr;


    // Private member functions

        //- Construct from components
        inline realGasConstTransport
        (
            const thermo& t,
            const scalar nu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline realGasConstTransport(const word&, const realGasConstTransport&);

        //- Construct from Istream
        realGasConstTransport(Istream&);


    // Member functions

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar rho,const scalar T) const;

        //- Thermal diffusivity for enthalpy [kg/ms]
        inline scalar alpha(const scalar rho,const scalar T) const;




    // Member operators

        inline realGasConstTransport& operator=
        (
            const realGasConstTransport&
        );


    // Friend operators

        friend realGasConstTransport operator+ <thermo>
        (
            const realGasConstTransport&,
            const realGasConstTransport&
        );

        friend realGasConstTransport operator- <thermo>
        (
            const realGasConstTransport&,
            const realGasConstTransport&
        );

        friend realGasConstTransport operator* <thermo>
        (
            const scalar,
            const realGasConstTransport&
        );

        friend realGasConstTransport operator== <thermo>
        (
            const realGasConstTransport&,
            const realGasConstTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const realGasConstTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "realGasConstTransportI.H"

#ifdef NoRepository
#   include "realGasConstTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
