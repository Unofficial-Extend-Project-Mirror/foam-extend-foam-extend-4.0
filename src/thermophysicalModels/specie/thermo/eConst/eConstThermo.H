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
    Foam::eConstThermo

Description
    Constant properties thermodynamics package derived from the basic
    thermo package data type specieThermo.

SourceFiles
    eConstThermoI.H
    eConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef eConstThermo_H
#define eConstThermo_H

#include "specieThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class eConstThermo;

template<class equationOfState>
inline eConstThermo<equationOfState> operator+
(
    const eConstThermo<equationOfState>&,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
inline eConstThermo<equationOfState> operator-
(
    const eConstThermo<equationOfState>&,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
inline eConstThermo<equationOfState> operator*
(
    const scalar,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
inline eConstThermo<equationOfState> operator==
(
    const eConstThermo<equationOfState>&,
    const eConstThermo<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const eConstThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                        Class eConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class eConstThermo
:
    public equationOfState
{
    // Private data

        scalar Cv_;
        scalar Hf_;


    // Private member functions

        //- Construct from components
        inline eConstThermo
        (
            const equationOfState& st,
            const scalar cv,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from Istream
        eConstThermo(Istream&);

        //- Construct as named copy
        inline eConstThermo(const word&, const eConstThermo&);

        //- Construct and return a clone
        inline autoPtr<eConstThermo> clone() const;

        // Selector from Istream
        inline static autoPtr<eConstThermo> New(Istream& is);


    // Member Functions

        // Fundamaental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar T) const;

            //- Enthalpy [J/kmol]
            inline scalar h(const scalar T) const;

            //- Sensible Enthalpy [J/kmol]
            inline scalar hs(const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar T) const;


    // Member operators

        inline void operator+=(const eConstThermo&);
        inline void operator-=(const eConstThermo&);


    // Friend operators

        friend eConstThermo operator+ <equationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );

        friend eConstThermo operator- <equationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );

        friend eConstThermo operator* <equationOfState>
        (
            const scalar,
            const eConstThermo&
        );

        friend eConstThermo operator== <equationOfState>
        (
            const eConstThermo&,
            const eConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const eConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eConstThermoI.H"

#ifdef NoRepository
#   include "eConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
