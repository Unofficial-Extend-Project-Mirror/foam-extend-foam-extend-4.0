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
    Foam::globalProcPointIndex

Description
    The class creates a unique global point index for each processor points
    (pair) in the mesh. Master and slave processor points carry the same index.

    Point offsets counts a number of unique points on each processor,
    excluding slave processor patch points and global points.

    Currently, faces are ordered with internal faces first, followed by patch
    faces in patch order, excluding slave processor patches.
    If needed, this can be changed to group processor faces with internal faces
    to facilitate parallel I/O.  HJ, 4/May/2018

Author
    Henrik Rusche, Wikki GmbH

SourceFiles
    globalProcPointIndex.C

\*---------------------------------------------------------------------------*/

#ifndef globalProcPointIndex_H
#define globalProcPointIndex_H

#include "polyMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class globalProcPointIndex Declaration
\*---------------------------------------------------------------------------*/

class globalProcPointIndex
{
    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Processor point index offset
        labelList procPointOffset_;

        //- point index offset for patches in local mesh
        labelList patchPointOffset_;

        //- Global poins label for all points of current mesh
        //  Sized to number of live points in the mesh
        labelList globalLabel_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalProcPointIndex(const globalProcPointIndex&);

        //- Disallow default bitwise assignment
        void operator=(const globalProcPointIndex&);

        //- Calculate point index
        void calcPointIndex();


public:

    // Constructors

        //- Construct from mesh
        globalProcPointIndex(const polyMesh&);


    //- Destructor - default


    // Member Functions

        //- Return point index offset per processor
        inline const labelList& procPointOffset() const
        {
            return procPointOffset_;
        }

        inline const labelList& patchPointOffset() const
        {
            return patchPointOffset_;
        }

        //- Return global point label for all points of current mesh
        //  Sized to number of live points in the mesh
        const labelList globalLabel() const
        {
            return globalLabel_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
