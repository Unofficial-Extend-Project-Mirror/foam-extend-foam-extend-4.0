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
    Foam::IDDESDelta

Description
    IDDESDelta used by the IDDES (improved low Re Spalart-Allmaras DES model)
    The min and max delta are calculated using the double distance of the min or
    max from the face centre to the cell centre.

SourceFiles
    IDDESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef IDDESDeltaDelta_H
#define IDDESDeltaDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IDDESDelta Declaration
\*---------------------------------------------------------------------------*/

class IDDESDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_;
        scalar cw_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        IDDESDelta(const IDDESDelta&);
        void operator=(const IDDESDelta&);

        //- Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("IDDESDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        IDDESDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    // Destructor

        ~IDDESDelta()
        {}


    // Member Functions

        //- Read the LESdelta dictionary
        void read(const dictionary&);

        // Correct values
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
