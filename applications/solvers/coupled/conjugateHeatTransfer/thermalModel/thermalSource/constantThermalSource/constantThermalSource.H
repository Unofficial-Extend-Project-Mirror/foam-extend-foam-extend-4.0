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
    constantThermalSource

Description
    Constant thermal properties

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

SourceFiles
    constantThermalSource.C

\*---------------------------------------------------------------------------*/

#ifndef constantThermalSource_H
#define constantThermalSource_H

#include "thermalSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantThermalSource Declaration
\*---------------------------------------------------------------------------*/

class constantThermalSource
:
    public thermalSource
{
    // Private data

        //- Strengh of the source
        dimensionedScalar S_;

        //- list of cell zones
        const wordList zones_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantThermalSource(const constantThermalSource&);

        //- Disallow default bitwise assignment
        void operator=(const constantThermalSource&);


public:

    //- Runtime type information
    TypeName("constant");

    // Static data members


    // Constructors

        //- Construct from dictionary
        constantThermalSource
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~constantThermalSource();


    // Member Functions

        //- Add source term
        virtual void addSource(volScalarField& source) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
