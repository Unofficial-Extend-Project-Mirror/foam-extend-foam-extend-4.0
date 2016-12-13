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
    Foam::coarseMgMeshLevel

Description
    Coarse level of multigrid hierarchy, holding addressing infromation
    for a clustered mesh

SourceFiles
    coarseMgMeshLevel.C

\*---------------------------------------------------------------------------*/

#ifndef coarseMgMeshLevel_H
#define coarseMgMeshLevel_H

#include "mgMeshLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coarseMgMeshLevel Declaration
\*---------------------------------------------------------------------------*/

class coarseMgMeshLevel
:
    public mgMeshLevel
{
    // Private data

        //- Reference to fine level
        const mgMeshLevel& fineLevel_;


        // Sizes

            //- Number of cells
            label nCells_;

            //- Number of internal faces
            label nInternalFaces_;


        // Addressing

            //- Owner addressing
            labelList owner_;

            //- Neighbour addressing
            labelList neighbour_;

            //- Face-cell addressing for patches
            labelListList faceCells_;


        // Geometrical data

            //- Cell centres
            vectorField cellCentres_;

            //- Face centres
            vectorField faceCentres_;

            //- Cell volumes
            scalarField cellVolumes_;

            //- Face area vectors
            vectorField faceAreas_;

            //- Face area magnitudes
            scalarField magFaceAreas_;


        // Restriction data

            //- Internal face restriction data
            //  Size = nFineFaces, giving coarse face index for each fine face
            //  Negative coefficient indicates the face is inside a coarse cell
            labelList faceRestrictAddr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coarseMgMeshLevel(const coarseMgMeshLevel&);

        //- Disallow default bitwise assignment
        void operator=(const coarseMgMeshLevel&);


        //- Calculate coarse level data from fine level and child array
        void calcLevel();

        //- Restrict face-based data (from fine to coarse)
         template<class Type>
         void restrictFace
         (
             UList<Type>& coarseData,
             const UList<Type>& fineData
         ) const
         {
             if
             (
                 coarseData.size() != nInternalFaces()
              || fineData.size() != fineLevel_.nInternalFaces()
             )
             {
                 FatalErrorIn
                 (
                     "void coarseMgMeshLevel::restrictFace\n"
                     "(\n"
                     "    UList<Type>& coarseData,\n"
                     "    const UList<Type>& fineData\n"
                     ") const\n"
                 )   << "Incorrect data sizes ("
                     << coarseData.size() << " " << fineData.size()
                     << ").  Should be ("
                     << nInternalFaces() << " "
                     << fineLevel_.nInternalFaces() << ")"
                     << abort(FatalError);
             }

             // Get addressing
             const labelList& c = faceRestrictAddr_;

             // Restriction = summation.  Set coarse data to zero
             coarseData = pTraits<Type>::zero;

             forAll (c, i)
             {
                 // Check address: only positive indices agglomerate
                 // into coarse faces
                 if (c[i] >= 0)
                 {
                     coarseData[c[i]] += fineData[i];
                 }
             }
         }


public:

        //- Runtime type information
        TypeName("fineMgMeshLevel");


    // Constructors

        //- Construct from fine level
        coarseMgMeshLevel(const mgMeshLevel& fineLevel);


    //- Destructor

        virtual ~coarseMgMeshLevel()
        {}


    // Member Functions

        //- Is this the finest level?
        virtual bool finest() const
        {
            return false;
        }

        //- Return access to fine level
        virtual const mgMeshLevel& fineLevel() const
        {
            return fineLevel_;
        }

        //- Return access to mesh
        virtual const fvMesh& mesh() const
        {
            return fineLevel().mesh();
        }


        // Sizes

            //- Return number of cells
            virtual label nCells() const
            {
                return nCells_;
            }

            //- Return number of internal faces
            virtual label nInternalFaces() const
            {
                return nInternalFaces_;
            }

            //- Return number of patches
            virtual label nPatches() const
            {
                return fineLevel_.nPatches();
            }


        // Addressing

            //- Return owner addressing
            virtual const unallocLabelList& owner() const
            {
                return owner_;
            }

            //- Return neighbour addressing
            const unallocLabelList& neighbour() const
            {
                return neighbour_;
            }

            //- Return face-cell addressing for all patches
            virtual const unallocLabelList& faceCells
            (
                const label patchNo
            ) const
            {
                return faceCells_[patchNo];
            }


        // Geometrical data

            //- Return number of valid geometric dimensions in the mesh
            virtual label nGeometricD() const
            {
                return fineLevel_.nGeometricD();
            }

            //- Return cell centres
            virtual const vectorField& cellCentres() const
            {
                return cellCentres_;
            }

            //- Return face centres
            virtual const vectorField& faceCentres() const
            {
                return faceCentres_;
            }

            //- Return cell volumes
            virtual const scalarField& cellVolumes() const
            {
                return cellVolumes_;
            }

            //- Return face area vectors
            virtual const vectorField& faceAreas() const
            {
                return faceAreas_;
            }

            //- Return face area magnitudes
            virtual const scalarField& magFaceAreas() const
            {
                return magFaceAreas_;
            }

            //- Return face area vectors for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual const vectorField& patchFaceAreas
            (
                const label patchNo
            ) const
            {
                return fineLevel().patchFaceAreas(patchNo);
            }

            //- Return face area magnitudes for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual const scalarField& magPatchFaceAreas
            (
                const label patchNo
            ) const
            {
                return fineLevel().magPatchFaceAreas(patchNo);
            }

            //- Return cell to face vectors for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaR
            (
                const label patchNo
            ) const
            {
                return fineLevel().patchDeltaR(patchNo);
            }

            //- Return neighbour cell to face vectors for coupled patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaRNeighbour
            (
                const label patchNo
            ) const
            {
                return fineLevel().patchDeltaRNeighbour(patchNo);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
