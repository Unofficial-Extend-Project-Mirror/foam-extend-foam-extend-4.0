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
    Foam::hConstThermo

Description
    Constant properties thermodynamics package
    templated into the equationOfState.

SourceFiles
    hConstThermoI.H
    hConstThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hConstThermo_H
#define hConstThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class hConstThermo;

template<class equationOfState>
inline hConstThermo<equationOfState> operator+
(
    const hConstThermo<equationOfState>&,
    const hConstThermo<equationOfState>&
);

template<class equationOfState>
inline hConstThermo<equationOfState> operator-
(
    const hConstThermo<equationOfState>&,
    const hConstThermo<equationOfState>&
);

template<class equationOfState>
inline hConstThermo<equationOfState> operator*
(
    const scalar,
    const hConstThermo<equationOfState>&
);

template<class equationOfState>
inline hConstThermo<equationOfState> operator==
(
    const hConstThermo<equationOfState>&,
    const hConstThermo<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const hConstThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hConstThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class hConstThermo
:
    public equationOfState
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private member functions

        //- Construct from components
        inline hConstThermo
        (
            const equationOfState& st,
            const scalar cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from Istream
        hConstThermo(Istream&);

        //- Construct as named copy
        inline hConstThermo(const word&, const hConstThermo&);

        //- Construct and return a clone
        inline autoPtr<hConstThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<hConstThermo> New(Istream& is);


    // Member Functions

        // Fundamaental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar T) const;

            //- Enthalpy [J/kmol]
            inline scalar h(const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar T) const;


    // Member operators

        inline void operator+=(const hConstThermo&);
        inline void operator-=(const hConstThermo&);


    // Friend operators

        friend hConstThermo operator+ <equationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator- <equationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );

        friend hConstThermo operator* <equationOfState>
        (
            const scalar,
            const hConstThermo&
        );

        friend hConstThermo operator== <equationOfState>
        (
            const hConstThermo&,
            const hConstThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const hConstThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hConstThermoI.H"

#ifdef NoRepository
#   include "hConstThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
