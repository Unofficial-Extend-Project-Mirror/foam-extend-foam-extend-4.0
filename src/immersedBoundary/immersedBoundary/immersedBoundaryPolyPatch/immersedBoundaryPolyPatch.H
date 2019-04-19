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
    Foam::immersedBoundaryPolyPatch

Description
    Immersed boundary patch

Author
    Zeljko Tukovic
    Reorganisation by Hrvoje Jasak
    Complete rewrite of methodology, Hrvoje Jasak

Contributors
    Inno Gatin, FMENA, Zagreb. All rights reserved.

SourceFiles
    immersedBoundaryPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryPolyPatch_H
#define immersedBoundaryPolyPatch_H

#include "polyPatch.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "standAlonePatch.H"
#include "Switch.H"
#include "volFields.H"
#include "immersedBoundaryCorrectedMeshFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryPolyPatch
:
    public polyPatch
{
    // Private data

        //- Triangular surface representing immersed boundary.
        //  Name of tri surface will be identical to the name of the patch
        triSurfaceMesh ibMesh_;

        //- Internal or external flow calculation
        Switch internalFlow_;

        //- Is this a wall walls
        Switch isWall_;

        //- Moving immersed boundary
        bool movingIb_;

        //- Time index for last update of mesh or moving boundary
        mutable label ibUpdateTimeIndex_;


    // Static data

        //- Search span factor
        //  Factor used to multiply the span of the cell to get the tri-seach
        //  span
        static const debug::tolerancesSwitch spanFactor_;


    // Demand-driven data

        //- Triangular surface search algorithm
        mutable triSurfaceSearch* triSurfSearchPtr_;


        // Immersed boundary patch data

            //- PrimitivePatch of cell-intersected immersed boundary faces
            //  Faces are calculated as intesection between the volume mesh
            //  and the ibMesh_
            mutable standAlonePatch* ibPatchPtr_;

            //- Cells intersected by the ibMesh_
            // This also corresponds to faceCells next to the IB patch
            mutable labelList* ibCellsPtr_;

            //- Centres of live part of intersected cells
            mutable vectorField* ibCellCentresPtr_;

            //- Volumes of live part of intersected cells
            mutable scalarField* ibCellVolumesPtr_;

            //- Intersected faces of the mesh
            mutable labelList* ibFacesPtr_;

            //- Centres of live part of intesected faces
            mutable vectorField* ibFaceCentresPtr_;

            //- Areas of live part of intesected faces
            mutable vectorField* ibFaceAreasPtr_;

            //- Nearest triangle for all immersed boundary patch faces
            mutable labelList* nearestTriPtr_;


        // Immersed boundary blanking

            //- Dead cells list
            mutable labelList* deadCellsPtr_;

            //- Dead faces list
            mutable labelList* deadFacesPtr_;


        // Corrected mesh geometry

            // MeshObject derived object that holds corrected global field data
            const immersedBoundaryCorrectedMeshFields& correctedFields_;

            //- Corrected face areas for the immersed boundary patch
            mutable vectorField* correctedIbPatchFaceAreasPtr_;

            //- Topo change index, used to update the corrected fields
            //  MeshObject consistently.
            //  NOTE: time index cannot be used since multiple topo changes per
            //  time-step are possible (IG 28/OCT/2018)
            mutable label topoChangeIndex_;


        // Immersed boundary motion

            //- Old points
            mutable pointField* oldIbPointsPtr_;


    // Private Member Functions

        //- Check cell intersection
        vector cellSpan(const label cellID) const;


        // Make demand-driven data

            //- Make triSurface search algorithm
            void calcTriSurfSearch() const;

            //- Calculate the immersed boundary patch geometry
            void calcImmersedBoundary() const;

            //- Calculate the corrected mesh geometry
            void calcCorrectedGeometry() const;


protected:

    // Protected Member Functions

        //- Correct patch after moving points
        virtual void movePoints(const pointField& p);


public:

    //- Runtime type information
    TypeName("immersedBoundary");


    // Constructors

        //- Construct from components
        immersedBoundaryPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        immersedBoundaryPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as copy
        immersedBoundaryPolyPatch(const immersedBoundaryPolyPatch&);

        //- Construct as copy, resetting the boundary mesh
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~immersedBoundaryPolyPatch();


    // Member Functions

        // Access

            //- Return immersed boundary surface mesh
            const triSurfaceMesh& ibMesh() const
            {
                return ibMesh_;
            }

            //- Return true if solving for flow inside the immersed boundary
            bool internalFlow() const
            {
                return internalFlow_;
            }

            // Is this a wall patch?
            virtual bool isWall() const
            {
                return isWall_;
            }

            //- Return triSurface search object
            const triSurfaceSearch& triSurfSearch() const;

            //- Return true if immersed boundary is moving
            bool movingIb() const
            {
                return movingIb_;
            }


        // Immersed boundary patch data

            //- Is the immersed boundary patch active?
            //  (is the intersection patch calculated?)
            bool active() const
            {
                return ibPatchPtr_ != nullptr;
            }

            //- Return primitivePatch of cell-intersected faces
            //  Faces are calculated as intesection between the volume mesh
            //  and the ibMesh
            const standAlonePatch& ibPatch() const;

            //- Return cells intersected by the ibMesh
            //  This also corresponds to faceCells next to the IB patch
            const labelList& ibCells() const;

            //- Return centres of live part of intersected cells
            const vectorField& ibCellCentres() const;

            //- Return volumes of live part of intersected cells
            const scalarField& ibCellVolumes() const;

            //- Return intersected faces of the mesh
            const labelList& ibFaces() const;

            //- Return centres of live part of intesected faces
            const vectorField& ibFaceCentres() const;

            //- Return areas of live part of intesected faces
            const vectorField& ibFaceAreas() const;

            //- Return nearest triangle for all immersed boundary patch faces
            const labelList& nearestTri() const;


        // Immersed boundary blanking

            //- Return dead cells
            const labelList& deadCells() const;

            //- Return dead faces
            const labelList& deadFaces() const;


        // Corrected mesh geometry

            //- Corrected cell centres
            const vectorField& correctedCellCentres() const;

            //- Corrected face centres
            const vectorField& correctedFaceCentres() const;

            //- Corrected cell volumes
            const scalarField& correctedCellVolumes() const;

            //- Corrected face areas
            const vectorField& correctedFaceAreas() const;


        // Corrected immersed patch geometry

            //- Corrected ib patch face areas
            const vectorField& correctedIbPatchFaceAreas() const;


        // Immersed boundary motion

            //- Correct patches after moving points
            void moveTriSurfacePoints(const pointField& p);

            //- Return old IB points
            const pointField& oldIbPoints() const;

            //- Motion distance vector for triangular immersed boundary faces
            tmp<vectorField> triMotionDistance() const;

            //- Motion distance vector for active immersed boundary faces
            tmp<vectorField> motionDistance() const;


        // Clear

            //- Clear geometry
            virtual void clearGeom();

            //- Clear addressing
            virtual void clearAddressing();

            //- Clear all demand-driven data
            void clearOut() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
