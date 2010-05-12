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
    Foam::jplot

Description

SourceFiles
    jplot.C

\*---------------------------------------------------------------------------*/

#ifndef jplot_H
#define jplot_H

#include "writer.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class jplot Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jplot
:
    public writer<Type>
{

    // Private Member Functions

        //- Write header
        Ostream& writeHeader(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("jplot");


    // Constructors

        //- Construct null
        jplot();

    // Destructor

        virtual ~jplot();


    // Member Functions

        // Access

        // Write

            virtual fileName getFileName
            (
                const coordSet&,
                const wordList&
            ) const;

            void write
            (
                const coordSet&,
                const wordList&,
                const List<const Field<Type>*>&,
                Ostream& os
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "jplot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
