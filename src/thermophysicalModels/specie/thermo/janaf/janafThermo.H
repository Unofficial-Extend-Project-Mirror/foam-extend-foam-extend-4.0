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
    Foam::janafThermo

Description
    JANAF tables based thermodynamics package templated
    into the equationOfState.

SourceFiles
    janafThermoI.H
    janafThermo.C

\*---------------------------------------------------------------------------*/

#ifndef janafThermo_H
#define janafThermo_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class janafThermo;

template<class equationOfState>
inline janafThermo<equationOfState> operator+
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator-
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator*
(
    const scalar,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator==
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const janafThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janafThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class janafThermo
:
    public equationOfState
{

public:

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private member functions

        //- Check given temperature is within the range of the fitted coeffs
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline void checkT(scalar& T) const;

        //- Return the coefficients corresponding to the given temperature
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline const coeffArray& coeffs(scalar& T) const;


public:

    // Constructors

        //- Construct from components
        inline janafThermo
        (
            const equationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highCpCoeffs,
            const coeffArray& lowCpCoeffs
        );

        //- Construct from Istream
        janafThermo(Istream&);

        //- Construct as a named copy
        inline janafThermo(const word&, const janafThermo&);


    // Member Functions

        //- Heat capacity at constant pressure [J/(kmol K)]
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline scalar cp(scalar T) const;

        //- Enthalpy [J/kmol]
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline scalar h(scalar T) const;

        //- Sensible enthalpy [J/kmol]
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline scalar hs(scalar T) const;

        //- Chemical enthalpy [J/kmol]
        inline scalar hc() const;

        //- Entropy [J/(kmol K)]
        //  Note: bounding T within range.  HJ, 12/Oct/2010
        inline scalar s(scalar T) const;


    // Member operators

        inline void operator+=(const janafThermo&);
        inline void operator-=(const janafThermo&);


    // Friend operators

        friend janafThermo operator+ <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );

        friend janafThermo operator- <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );

        friend janafThermo operator* <equationOfState>
        (
            const scalar,
            const janafThermo&
        );

        friend janafThermo operator== <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const janafThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janafThermoI.H"

#ifdef NoRepository
#   include "janafThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
