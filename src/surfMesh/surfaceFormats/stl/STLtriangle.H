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
    Foam::STLtriangle

Description
    A triangle representation for STL files.

SourceFiles
    STLtriangleI.H

\*---------------------------------------------------------------------------*/

#ifndef STLtriangle_H
#define STLtriangle_H

#include "STLpoint.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class STLtriangle Declaration
\*---------------------------------------------------------------------------*/

class STLtriangle
{
    // Private data

        //- Attribute is 16-bit
        typedef unsigned short STLattrib;

        //- The face normal, many programs write zore or other junk
        STLpoint normal_;

        //- The three points defining the triangle
        STLpoint a_, b_, c_;

        //- The attribute information could for colour or solid id, etc
        STLattrib attrib_;


public:

    // Constructors

        //- Construct null
        inline STLtriangle();

        //- Construct from components
        inline STLtriangle
        (
            const STLpoint& normal,
            const STLpoint& a,
            const STLpoint& b,
            const STLpoint& c,
            unsigned short attrib
        );

        //- Construct from istream (read binary)
        inline STLtriangle(istream&);


    // Member Functions

        // Access

            inline const STLpoint& normal() const;
            inline const STLpoint& a() const;
            inline const STLpoint& b() const;
            inline const STLpoint& c() const;
            inline unsigned short attrib() const;


        // Read

            //- Read from istream (binary)
            inline void read(istream&);


        // Write

            //- Write to ostream (binary)
            inline void write(ostream&);


    // Ostream operator

        inline friend Ostream& operator<<(Ostream&, const STLtriangle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "STLtriangleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
