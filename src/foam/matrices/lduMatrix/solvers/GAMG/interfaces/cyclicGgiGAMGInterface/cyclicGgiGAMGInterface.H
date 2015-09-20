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
    Foam::cyclicGgiGAMGInterface

Description
    GAMG agglomerated cyclicGgi interface.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    cyclicGgiGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGgiGAMGInterface_H
#define cyclicGgiGAMGInterface_H

#include "GAMGInterface.H"
#include "ggiGAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cyclicGgiGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicGgiGAMGInterface
:
    public ggiGAMGInterface
{
public:

    //- Runtime type information
    TypeName("cyclicGgi");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicGgiGAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    // Destructor

        virtual ~cyclicGgiGAMGInterface();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
