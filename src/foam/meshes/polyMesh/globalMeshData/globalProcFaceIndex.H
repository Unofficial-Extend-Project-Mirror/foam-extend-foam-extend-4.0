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
    Foam::globalProcFaceIndex

Description
    The class creates a unique global face index for each processor face (pair)
    in the mesh.  Master and slave processor face carry the same index.

    The data is presented as a labelListList, which is set only for
    processor patches.

    Face offsets counts a number of unique faces on each processor,
    excluding slave processor patch faces, which are given master face index.
    If needed, global face index from all faces can be derived from this data

    Currently, faces are ordered with internal faces first, followed by patch
    faces in patch order, excluding slave processor patches.
    If needed, this can be changed to group processor faces with internal faces
    to facilitate parallel I/O.  HJ, 4/May/2018

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    globalProcFaceIndex.C

\*---------------------------------------------------------------------------*/

#ifndef globalProcFaceIndex_H
#define globalProcFaceIndex_H

#include "polyMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class globalProcFaceIndex Declaration
\*---------------------------------------------------------------------------*/

class globalProcFaceIndex
{
    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Number of unique faces per processor
        labelList nUniqueFaces_;

        //- Processor face index offset
        labelList procFaceOffset_;

        //- Global face label for all faces of current mesh
        //  Sized to number of live faces in the mesh
        labelList globalLabel_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalProcFaceIndex(const globalProcFaceIndex&);

        //- Disallow default bitwise assignment
        void operator=(const globalProcFaceIndex&);

        //- Calculate face index
        void calcFaceIndex();


public:

    // Constructors

        //- Construct from mesh
        globalProcFaceIndex(const polyMesh&);


    //- Destructor - default


    // Member Functions

        //- Return number of unique faces per processor
        inline const labelList& nUniqueFaces() const
        {
            return nUniqueFaces_;
        }

        //- Return face index offset per processor
        inline const labelList& procFaceOffset() const
        {
            return procFaceOffset_;
        }

        //- Return global face label for all faces of current mesh
        //  Sized to number of live faces in the mesh
        const labelList& globalLabel() const
        {
            return globalLabel_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
