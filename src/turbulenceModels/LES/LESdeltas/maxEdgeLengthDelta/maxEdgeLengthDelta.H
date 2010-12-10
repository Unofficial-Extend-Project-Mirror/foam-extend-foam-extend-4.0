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
    Foam::maxEdgeLengthDelta

Description
    Maximum cell edge length delta used in LES models.

SourceFiles
    maxEdgeLengthDelta.C

\*---------------------------------------------------------------------------*/

#ifndef maxEdgeLengthDelta_H
#define maxEdgeLengthDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class maxEdgeLengthDelta Declaration
\*---------------------------------------------------------------------------*/

class maxEdgeLengthDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        maxEdgeLengthDelta(const maxEdgeLengthDelta&);

        //- Disallow default bitwise assignment
        void operator=(const maxEdgeLengthDelta&);


        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("maxEdgeLength");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        maxEdgeLengthDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    // Destructor

        virtual ~maxEdgeLengthDelta()
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
