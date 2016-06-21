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
    Foam::fineMgMeshLevel

Description
    Fine multigrid mesh level.  Pach fvMesh information into mgMeshLevel
    interface.  Allow for access to fvMesh at the fine level.

SourceFiles
    fineMgMeshLevel.C

\*---------------------------------------------------------------------------*/

#ifndef fineMgMeshLevel_H
#define fineMgMeshLevel_H

#include "mgMeshLevel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fineMgMeshLevel Declaration
\*---------------------------------------------------------------------------*/

class fineMgMeshLevel
:
    public mgMeshLevel
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fineMgMeshLevel(const fineMgMeshLevel&);

        //- Disallow default bitwise assignment
        void operator=(const fineMgMeshLevel&);


public:

    // Static data

        //- Runtime type information
        TypeName("fineMgMeshLevel");


    // Constructors

        //- Construct from mesh
        fineMgMeshLevel(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    //- Destructor

        virtual ~fineMgMeshLevel()
        {}


    // Member Functions

        //- Is this the finest level?
        virtual bool finest() const
        {
            return true;
        }

        //- Return access to fine level
        virtual const mgMeshLevel& fineLevel() const
        {
            return *this;
        }

        //- Return access to mesh
        virtual const fvMesh& mesh() const
        {
            return mesh_;
        }


        // Sizes

            //- Return number of cells
            virtual label nCells() const
            {
                return mesh_.nCells();
            }

            //- Return number of internal faces
            virtual label nInternalFaces() const
            {
                return mesh_.nInternalFaces();
            }

            //- Return number of patches
            virtual label nPatches() const
            {
                return mesh_.boundary().size();
            }


        // Addressing

            //- Return owner addressing
            virtual const unallocLabelList& owner() const
            {
                return mesh_.owner();
            }

            //- Return neighbour addressing
            const unallocLabelList& neighbour() const
            {
                return mesh_.neighbour();
            }

            //- Return face-cell addressing for all patches
            virtual const unallocLabelList& faceCells
            (
                const label patchNo
            ) const
            {
                return mesh_.boundary()[patchNo].faceCells();
            }


        // Geometrical data

            //- Return number of valid geometric dimensions in the mesh
            virtual label nGeometricD() const
            {
                return mesh_.nGeometricD();
            }

            //- Return cell centres
            virtual const vectorField& cellCentres() const
            {
                return mesh_.C().internalField();
            }

            //- Return face centres
            virtual const vectorField& faceCentres() const
            {
                return mesh_.Cf().internalField();
            }

            //- Return cell volumes
            virtual const scalarField& cellVolumes() const
            {
                return mesh_.V().field();
            }

            //- Return face area vectors
            virtual const vectorField& faceAreas() const
            {
                return mesh_.Sf().internalField();
            }

            //- Return face area magnitudes
            virtual const scalarField& magFaceAreas() const
            {
                return mesh_.magSf().internalField();
            }

            //- Return face area vectors for patches
            virtual const vectorField& patchFaceAreas
            (
                const label patchNo
            ) const
            {
                return mesh_.Sf().boundaryField()[patchNo];
            }

            //- Return face area magnitudes for patches
            virtual const scalarField& magPatchFaceAreas
            (
                const label patchNo
            ) const
            {
                return mesh_.magSf().boundaryField()[patchNo];
            }

            //- Return cell to face vectors for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaR
            (
                const label patchNo
            ) const
            {
                return mesh_.Cf().boundaryField()[patchNo]
                    - mesh_.C().boundaryField()[patchNo].patchInternalField();
            }

            //- Return neighbour cell to face vectors for coupled patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaRNeighbour
            (
                const label patchNo
            ) const
            {
                // Check.  HJ, 4/Jun/2013
                return patchDeltaR(patchNo)
                    - mesh_.boundary()[patchNo].delta();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
