/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::passiveProcessorFvPatch

Description
    Passive processor fvPatch.  See passiveProcessorPolyPatch

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    passiveProcessorFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef passiveProcessorFvPatch_H
#define passiveProcessorFvPatch_H

#include "coupledFvPatch.H"
#include "processorLduInterface.H"
#include "passiveProcessorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class passiveProcessorFvPatch Declaration
\*---------------------------------------------------------------------------*/

class passiveProcessorFvPatch
:
    public fvPatch
{
public:

    //- Runtime type information
    TypeName(passiveProcessorPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        passiveProcessorFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(patch, bm)
        {}


    //- Destructor
    virtual ~passiveProcessorFvPatch()
    {}


    // Member functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
