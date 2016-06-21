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
    Foam::hexCellLooper

Description
    Implementation of cellLooper.

    This one walks hexes in a topological way:
      - cross edge to other face
      - cross face by walking edge-point-edge across to reach the other side.
    (edges are always cut through the middle)

    For anything else (tet, prism, .. poly) it will use geomCellLooper
    (which does a purely geometric cut using a plane through cell centre)

SourceFiles
    hexCellLooper.C

\*---------------------------------------------------------------------------*/

#ifndef hexCellLooper_H
#define hexCellLooper_H

#include "geomCellLooper.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cellModel;

/*---------------------------------------------------------------------------*\
                           Class hexCellLooper Declaration
\*---------------------------------------------------------------------------*/

class hexCellLooper
:
    public geomCellLooper
{

protected:

    // Protected data

        //- Reference to hex cell shape.
        const cellModel& hex_;

private:

    // Private Member Functions

        //- walk across faces of hex. Update loop/loopWeights with edges cut.
        //  Return true if successful walk. (should be always!)
        bool walkHex
        (
            const label cellI,
            const label startFaceI,
            const label startEdgeI,

            labelList& loop,
            scalarField& loopWeights
        ) const;

        //- Convert loop into face and points
        void makeFace
        (
            const labelList& loop,
            const scalarField& loopWeights,

            labelList& faceVerts,
            pointField& facePoints
        ) const;


        //- Disallow default bitwise copy construct
        hexCellLooper(const hexCellLooper&);

        //- Disallow default bitwise assignment
        void operator=(const hexCellLooper&);


public:

    //- Runtime type information
    TypeName("hexCellLooper");


    // Constructors

        //- Construct from components
        hexCellLooper(const polyMesh& mesh);


    // Destructor

        virtual ~hexCellLooper();


    // Member Functions

        //- Create cut along circumference of cellI. Gets current mesh cuts.
        //  Cut along circumference is expressed as loop of cuts plus weights
        //  for cuts along edges (only valid for edge cuts).
        //  Return true if successful cut.
        virtual bool cut
        (
            const vector& refDir,
            const label cellI,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;


        //- Same but now also base point of cut provided (instead of always
        //  cell centre)
        virtual bool cut
        (
            const plane& cutPlane,
            const label cellI,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
