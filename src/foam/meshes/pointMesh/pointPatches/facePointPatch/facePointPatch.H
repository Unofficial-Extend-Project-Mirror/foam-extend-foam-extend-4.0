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
    Foam::facePointPatch

Description
    A pointPatch based on a polyPatch

SourceFiles
    facePointPatch.C
    facePointPatchM.C
    newPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef facePointPatch_H
#define facePointPatch_H

#include "pointPatch.H"
#include "polyPatch.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class facePointPatch Declaration
\*---------------------------------------------------------------------------*/

class facePointPatch
:
    public pointPatch
{
    // Private data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


    // Private Member Functions

        //- Disallow default bitwise copy construct
        facePointPatch(const facePointPatch&);

        //- Disallow default bitwise assignment
        void operator=(const facePointPatch&);


public:



    //- Runtime type information
    TypeName(polyPatch::typeName_());


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            facePointPatch,
            polyPatch,
            (const polyPatch& patch, const pointBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch
        facePointPatch
        (
            const polyPatch&,
            const pointBoundaryMesh&
        );


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<facePointPatch> New
        (
            const polyPatch&,
            const pointBoundaryMesh&
        );


    // Destructor

        virtual ~facePointPatch()
        {}


    // Member Functions

        //- Return the polyPatch
        const polyPatch& patch() const
        {
            return polyPatch_;
        }

        //- Return name
        virtual const word& name() const
        {
            return polyPatch_.name();
        }

        //- Return size
        virtual label size() const
        {
            return polyPatch_.nPoints();
        }

        //- Return number of faces
        virtual label nFaces() const
        {
            return polyPatch_.size();
        }

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const
        {
            return polyPatch_.index();
        }


        // Access functions for demand driven data

            //- Return mesh points
            virtual const labelList& meshPoints() const
            {
                return polyPatch_.meshPoints();
            }

            //- Return pointField of points in patch
            virtual const pointField& localPoints() const
            {
                return polyPatch_.localPoints();
            }

            //- Return point unit normals
            virtual const vectorField& pointNormals() const
            {
                return polyPatch_.pointNormals();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
