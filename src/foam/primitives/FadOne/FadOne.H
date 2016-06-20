/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    FadOne

Description
    Forward automatic derivative of the first order.  Supports multi-variable
    derivative, where the number of independent variables is specified at
    compile-time.  Replaces scalar in calculations

SourceFiles
    FadOneI.H
    FadOne.C

\*---------------------------------------------------------------------------*/

#ifndef FadOne_H
#define FadOne_H

#include "error.H"
#include "IOstream.H"
#include "token.H"

#include "products.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<int nVars> class FadOne;

template<int nVars>
inline FadOne<nVars> operator+(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator+(const double, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator+(const FadOne<nVars>&, const double);

template<int nVars>
inline FadOne<nVars> operator-(const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator-(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator-(const double, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator-(const FadOne<nVars>&, const double);

template<int nVars>
inline FadOne<nVars> operator*(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator*(const double, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator*(const FadOne<nVars>&, const double);

template<int nVars>
inline FadOne<nVars> operator/(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator/(const double, const FadOne<nVars>&);

template<int nVars>
inline FadOne<nVars> operator/(const FadOne<nVars>&, const double);

template<int nVars>
inline bool operator<(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline bool operator<=(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline bool operator>(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline bool operator>=(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline bool operator==(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
inline bool operator!=(const FadOne<nVars>&, const FadOne<nVars>&);

template<int nVars>
Istream& operator>>(Istream&, FadOne<nVars>&);

template<int nVars>
Ostream& operator<<(Ostream&, const FadOne<nVars>&);


/*---------------------------------------------------------------------------*\
                           Class FadOne Declaration
\*---------------------------------------------------------------------------*/

template<int nVars>
class FadOne
{
    // Private data

        //- Length of data list
        static const label Length = nVars + 1;

        //- Value and derivatives
        double data_[Length];


    // Private Member Functions


public:

    //- Component type
    typedef FadOne<nVars> cmptType;


    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank od Scalar is 0
            nComponents = 1  // Number of components in FadOne is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const FadOne<nVars> zero;
        static const FadOne<nVars> one;


    // Constructors

        //- Construct null
        inline FadOne();

        //- Construct from value
        explicit inline FadOne(const double& v);

        //- Construct from Istream
        explicit inline FadOne(Istream&);

        //- Construct as copy
        inline FadOne(const FadOne<nVars>&);


    // Destructor - default


    // Member Functions

        //- Return value
        inline const double& value() const;
        inline double& value();

        //- Return derivative
        inline const double& deriv(const label i) const;
        inline double& deriv(const label i);


    // Member Operators

        inline void operator=(const double&);
        inline void operator=(const FadOne<nVars>&);

        inline void operator+=(const double&);
        inline void operator+=(const FadOne<nVars>&);
        inline void operator-=(const double&);
        inline void operator-=(const FadOne<nVars>&);
        inline void operator*=(const double&);
        inline void operator*=(const FadOne<nVars>&);
        inline void operator/=(const double&);
        inline void operator/=(const FadOne<nVars>&);


    // Friend Functions


    // Friend operators

        friend FadOne<nVars> operator+ <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator+ <nVars>
        (
            const double,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator+ <nVars>
        (
            const FadOne<nVars>&,
            const double
        );

        friend FadOne<nVars> operator- <nVars>(const FadOne<nVars>&);

        friend FadOne<nVars> operator- <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator- <nVars>
        (
            const double,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator- <nVars>
        (
            const FadOne<nVars>&,
            const double
        );

        friend FadOne<nVars> operator* <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator* <nVars>
        (
            const double,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator* <nVars>
        (
            const FadOne<nVars>&,
            const double
        );

        friend FadOne<nVars> operator/ <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator/ <nVars>
        (
            const double,
            const FadOne<nVars>&
        );

        friend FadOne<nVars> operator/ <nVars>
        (
            const FadOne<nVars>&,
            const double
        );

        friend bool operator< <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend bool operator<= <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend bool operator> <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend bool operator>= <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend bool operator== <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );

        friend bool operator!= <nVars>
        (
            const FadOne<nVars>&,
            const FadOne<nVars>&
        );


    // IOstream Operators

        friend Istream& operator>> <nVars>(Istream&, FadOne<nVars>&);
        friend Ostream& operator<< <nVars>(Ostream&, const FadOne<nVars>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FadOneI.H"
#include "FadOne.C"
#include "FadOneFunctions.H"
#include "FadOneFunctionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
