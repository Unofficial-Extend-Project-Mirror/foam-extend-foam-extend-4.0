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
    Foam::mgMeshLevel

Description
    Virtual base class for a level of multigrid hierarchy.

    Cells and internal faces are coarsened by clustering
    There is no clustering on the boundary: original faces are used and
    face-cell patch addressing on coarse levels returns addressing into
    coarse cell clusters.

SourceFiles
    mgMeshLevel.C

\*---------------------------------------------------------------------------*/

#ifndef mgMeshLevel_H
#define mgMeshLevel_H

#include "typeInfo.H"
#include "fvMesh.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class mgMeshLevel Declaration
\*---------------------------------------------------------------------------*/

class mgMeshLevel
{
    // Private data

        //- Number of coarse cells
        mutable label nCoarseCells_;

        //- Child array: coarse cluster index for every fine cell
        mutable labelList child_;


    // Private Member Functions

        //- Construct the CSR format addressing
        void makeCompactAddressingAndWeights
        (
            labelList& cellCells,
            labelList& cellCellOffsets,
            scalarField& faceWeights,
            scalarField& boundaryAreas
        ) const;

        //- Create child array
        void makeChild() const;


public:

    // Static Data

        //- Runtime type information
        TypeName("mgMeshLevel");


        // Coarsening parameters

            //- Minimum cluster size
            static const debug::optimisationSwitch mgMinClusterSize_;

            //- Maximum cluster size
            static const debug::optimisationSwitch mgMaxClusterSize_;


    // Constructor

        //- Construct null
        mgMeshLevel()
        :
            nCoarseCells_(-1),
            child_()
        {}


    //- Destructor

        virtual ~mgMeshLevel()
        {}


    // Member Functions

        //- Is this the finest level?
        virtual bool finest() const = 0;

        //- Return reference to fine level
        virtual const mgMeshLevel& fineLevel() const = 0;

        //- Return reference to finest level
        const mgMeshLevel& finestLevel() const;

        //- Return access to mesh
        virtual const fvMesh& mesh() const = 0;


        // Sizes

            //- Return number of cells
            virtual label nCells() const = 0;

            //- Return number of internal faces
            virtual label nInternalFaces() const = 0;

            //- Return number of patches
            virtual label nPatches() const = 0;


        // Addressing

            //- Return owner addressing
            virtual const unallocLabelList& owner() const = 0;

            //- Return neighbour addressing
            virtual const unallocLabelList& neighbour() const = 0;

            //- Return face-cell addressing for all patches
            virtual const unallocLabelList& faceCells
            (
                const label patchNo
            ) const = 0;


        // Geometrical data

            //- Return number of valid geometric dimensions in the mesh
            virtual label nGeometricD() const = 0;

            //- Return cell centres
            virtual const vectorField& cellCentres() const = 0;

            //- Return face centres
            virtual const vectorField& faceCentres() const = 0;

            //- Return cell volumes
            virtual const scalarField& cellVolumes() const = 0;

            //- Return face area vectors
            virtual const vectorField& faceAreas() const = 0;

            //- Return face area magnitudes
            virtual const scalarField& magFaceAreas() const = 0;

            //- Return face area vectors for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual const vectorField& patchFaceAreas
            (
                const label patchNo
            ) const = 0;

            //- Return face area magnitudes for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual const scalarField& magPatchFaceAreas
            (
                const label patchNo
            ) const = 0;

            //- Return cell to face vectors for patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaR
            (
                const label patchNo
            ) const = 0;

            //- Return neighbour cell to face vectors for coupled patches
            //  Note: no clustering of boundary faces allowed:
            //        using fine-level geometry
            virtual tmp<vectorField> patchDeltaRNeighbour
            (
                const label patchNo
            ) const = 0;


        // Coarsening

            //- Return number of cells
            label nCoarseCells() const;

            //- Return child array
            const labelList& child() const;

            //- Create next level from current level
            autoPtr<mgMeshLevel> makeNextLevel() const;


         // Restriction and prolongation

             //- Restrict (from fine to coarse).  Call from the level
             //  where data originates
             template<class Type>
             void restrict
             (
                 UList<Type>& coarseData,
                 const UList<Type>& fineData
             ) const
             {
                 if
                 (
                     coarseData.size() != nCoarseCells()
                  || fineData.size() != nCells()
                 )
                 {
                     FatalErrorIn
                     (
                         "void mgMeshLevel::restrict\n"
                         "(\n"
                         "    UList<Type>& coarseData,\n"
                         "    const UList<Type>& fineData\n"
                         ") const\n"
                     )   << "Incorrect data sizes ("
                         << coarseData.size() << " " << fineData.size()
                         << ").  Should be ("
                         << nCoarseCells() << " " << nCells() << ")"
                         << abort(FatalError);
                 }

                 // Get addressing
                 const labelList& c = child();

                 // Restriction = summation.  Set coarse data to zero
                 coarseData = pTraits<Type>::zero;

                 forAll (c, i)
                 {
                     coarseData[c[i]] += fineData[i];
                 }
             }

             //- Prolong (from  coarse to fine).  Call from the level
             //  where data originates
             template<class Type>
             void prolong
             (
                 const UList<Type>& coarseData,
                 UList<Type>& fineData
             ) const
             {
                 if (finest())
                 {
                     FatalErrorIn
                     (
                         "void mgMeshLevel::prolong\n"
                         "(\n"
                         "    UList<Type>& coarseData,\n"
                         "    const UList<Type>& fineData\n"
                         ") const\n"
                     )   << "Requested prolong from finest level"
                         << abort(FatalError);
                 }

                 if
                 (
                     coarseData.size() != nCells()
                  || fineData.size() != fineLevel().nCells()
                 )
                 {
                     FatalErrorIn
                     (
                         "void mgMeshLevel::prolong\n"
                         "(\n"
                         "    UList<Type>& coarseData,\n"
                         "    const UList<Type>& fineData\n"
                         ") const\n"
                     )   << "Incorrect data sizes ("
                         << coarseData.size() << " " << fineData.size()
                         << ").  Should be ("
                         << nCoarseCells() << " " << nCells() << ")"
                         << abort(FatalError);
                 }

                 // Get addressing
                 const labelList& c = fineLevel().child();

                 // Prolongation = assignment.  Initialisation not needed
                 forAll (c, i)
                 {
                     fineData[i] = coarseData[c[i]];
                 }
             }

             //- Prolong to the finest level (from  coarse to fine)
             template<class Type>
             tmp<Field<Type> > prolongToFinest
             (
                 const UList<Type>& coarseData
             ) const
             {
                 if (finest())
                 {
                     // Already on finest level
                     return tmp<Field<Type> >(new Field<Type>(coarseData));
                 }
                 else if (fineLevel().finest())
                 {
                     tmp<Field<Type> > tresult
                     (
                         new Field<Type>(fineLevel().nCells())
                     );

                     prolong(coarseData, tresult());

                     // Reached finest level, return
                     return tresult;
                 }
                 else
                 {
                     Field<Type> intermediate(fineLevel().nCells());

                     prolong(coarseData, intermediate);

                     // Recursive return
                     return fineLevel().prolongToFinest(intermediate);
                 }
             }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
