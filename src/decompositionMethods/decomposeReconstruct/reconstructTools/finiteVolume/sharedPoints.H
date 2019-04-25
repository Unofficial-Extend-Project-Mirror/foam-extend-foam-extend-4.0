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
    Foam::sharedPoints

Description
    Calculates points shared by more than two processor patches from a list
    of meshes.

    Works purely on topology.
    Needs:
      - domain to be one single domain (i.e. all faces can be reached through
        face-cell walk).
      - patch face ordering to be ok
      - f[0] ordering on patch faces to be ok.

    The class will create a unique list of multiply shared points
    (on more than 2 mesh pieces).  On each mesh piece, it will provide a
    a list of local point labels that are globally shared (sharedPointLabels)
    and for each of such points the index in the list of globally shared
    points (sharedPointAddr).  The class also provides the number of globally
    shared points across all mesh pieces (nGlobalPoints)

Note
    Currently operating with a PtrList<fvMesh>, whereas the operation actually
    works on PtrList<polyMesh> as well, due to the calling side

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    sharedPoints.C

\*---------------------------------------------------------------------------*/

#ifndef sharedPoints_H
#define sharedPoints_H

#include "PtrList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class sharedPoints Declaration
\*---------------------------------------------------------------------------*/

class sharedPoints
{
    // Private typedefs

        //- Define procPointList as holding a list of meshPoint/processor labels
        typedef FixedList<label, 2> procPoint;
        typedef List<procPoint> procPointList;


    // Private data

        //- List of meshes
        const PtrList<fvMesh>& meshes_;

        //- Shared points used by this processor (= global point number)
        labelListList sharedPointAddr_;

        //- My mesh points corresponding to the shared points
        labelListList sharedPointLabels_;

        //- Total number of shared points.
        label nGlobalPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sharedPoints(const sharedPoints&);

        //- Disallow default bitwise assignment
        void operator=(const sharedPoints&);


        //- Calculate processor patch pairs
        //  For each processor patch, find equivalent patch on other mesh
        labelListList procPatchPairs() const;

        //- Sync mark using processor patches
        void syncMark
        (
            labelListList& markedPoints,
            const labelListList& patchPairs,
            const label fromMesh
        ) const;

        //- Calculate shared points
        void calcSharedPoints();

        //- Calculate shared points
        void calcSharedPoints(const labelListList& globalPointIndex);


public:

    // Constructors

        //- Construct from the list of meshes
        sharedPoints
        (
            const PtrList<fvMesh>& meshes
        );

        //- Construct from the list of meshes and global point addressing
        sharedPoints
        (
            const PtrList<fvMesh>& meshes,
            const labelListList& globalPointIndex
        );


    //- Destructor
    ~sharedPoints()
    {}


    // Member Functions

        // Access

            //- Shared points used by this processor (= global point number)
            inline const labelListList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- Local mesh points corresponding to the shared points
            inline const labelListList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- Number of globally shared points
            inline label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
