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
    Foam::complex

Description
    Extension to the c++ complex library type.

SourceFiles
    complexI.H
    complex.C

\*---------------------------------------------------------------------------*/

#ifndef complex_H
#define complex_H

#include "scalar.H"
#include "bool.H"
#include "word.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class complex;

inline scalar magSqr(const complex& c);
inline complex sqr(const complex& c);
inline scalar mag(const complex& c);
inline const complex& max(const complex& c1, const complex& c2);
inline const complex& min(const complex& c1, const complex& c2);
inline complex limit(const complex& c1, const complex& c2);
inline const complex& sum(const complex& c);
inline complex operator+(const complex&, const complex&);
inline complex operator-(const complex&);
inline complex operator-(const complex&, const complex&);
inline complex operator*(const complex&, const complex&);
inline complex operator/(const complex&, const complex&);
inline complex operator*(const scalar, const complex&);
inline complex operator*(const complex&, const scalar);
inline complex operator/(const complex&, const scalar);
inline complex operator/(const scalar, const complex&);
Istream& operator>>(Istream& is, complex&);
Ostream& operator<<(Ostream& os, const complex& C);


/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/

class complex
{
    // private data

        //- Real and imaginary parts of the complex number
        scalar re, im;

public:

    //- Component type
    typedef complex cmptType;


    // Static data members

        static const char* const typeName;

        static const complex zero;
        static const complex one;


    // Constructors

        //- Construct null
        inline complex();

        //- Construct given real and imaginary parts
        inline complex(const scalar Re, const scalar Im);

        //- Construct from Istream
        complex(Istream&);


    // Member functions

           // Access

               inline scalar Re() const;
               inline scalar Im() const;

           // Edit

               inline scalar& Re();
               inline scalar& Im();

           // Operators

               inline complex conjugate() const;


    // Member operators

        inline void operator=(const complex&);
        inline void operator+=(const complex&);
        inline void operator-=(const complex&);
        inline void operator*=(const complex&);
        inline void operator/=(const complex&);

        inline void operator=(const scalar);
        inline void operator+=(const scalar);
        inline void operator-=(const scalar);
        inline void operator*=(const scalar);
        inline void operator/=(const scalar);

        inline complex operator!() const;

        inline bool operator==(const complex&) const;
        inline bool operator!=(const complex&) const;


    // Friend functions

        friend scalar magSqr(const complex& c);
        friend complex sqr(const complex& c);
        friend scalar mag(const complex& c);
        friend const complex& max(const complex& c1, const complex& c2);
        friend const complex& min(const complex& c1, const complex& c2);

        friend complex limit(const complex& c1, const complex& c2);

        friend const complex& sum(const complex& c);


    // Friend operators

        friend complex operator+(const complex&, const complex&);
        friend complex operator-(const complex&);
        friend complex operator-(const complex&, const complex&);
        friend complex operator*(const complex&, const complex&);
        friend complex operator/(const complex&, const complex&);

        friend complex operator*(const scalar, const complex&);
        friend complex operator*(const complex&, const scalar);
        friend complex operator/(const complex&, const scalar);
        friend complex operator/(const scalar, const complex&);


    // IOstream operators

        friend Istream& operator>>(Istream& is, complex&);
        friend Ostream& operator<<(Ostream& os, const complex& C);

};


// * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * * //

//- Return a string representation of a complex
word name(const complex&);


//- Specify data associated with complex type is contiguous
template<>
inline bool contiguous<complex>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
