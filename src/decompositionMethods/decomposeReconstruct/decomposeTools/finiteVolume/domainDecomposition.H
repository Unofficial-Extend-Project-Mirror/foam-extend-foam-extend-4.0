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
    Foam::domainDecomposition

Description
    Automatic domain decomposition class for FOAM meshes

SourceFiles
    domainDecomposition.C
    distributeCells.C
    decomposeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecomposition_H
#define domainDecomposition_H

#include "fvMesh.H"
#include "labelList.H"
#include "SLList.H"
#include "PtrList.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Mesh decomposition control dictionary
        dictionary decompositionDict_;

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Processor label for each cell
        labelList cellToProc_;

        //- Labels of points for each processor
        labelListList procPointAddressing_;

        //- Labels of faces for each processor
        //  Note: Face turning index is stored as the sign on addressing
        //  Only the processor boundary faces are affected: if the sign of the
        //  index is negative, the processor face is the reverse of the
        //  original face. In order to do this properly, all face
        //   indices will be incremented by 1 and the decremented as
        //  necessary t avoid the problem of face number zero having no
        //  sign.  HJ, 5/Dec/2001
        labelListList procFaceAddressing_;

        //- Number of internal faces for each processor
        labelList nInternalProcFaces_;

        //- Number of live faces for each processor
        labelList nLiveProcFaces_;

        //- Labels of cells for each processor
        labelListList procCellAddressing_;

        //- Original patch index for every processor patch
        labelListList procBoundaryAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- List of globally shared point labels
        labelList globallySharedPoints_;

        //- Are there cyclic-parallel faces
        bool cyclicParallel_;


    // Private Member Functions

        //- Create cell-to processor list using domain decomposition tools
        void distributeCells();


public:

    // Declare name of the class and its debug switch
    ClassName("domainDecomposition");

    // Constructors

        //- Construct from mesh and dictionary
        domainDecomposition
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~domainDecomposition();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Number of processor in decomposition
        label nProcs() const
        {
            return nProcs_;
        }

        //- Is the decomposition data to be distributed for each processor
        bool distributed() const
        {
            return distributed_;
        }

        //- Return cell-processor decomposition labels
        const labelList& cellToProc() const
        {
            return cellToProc_;
        }

        //- Decompose mesh. Optionally remove zero-sized patches.
        void decomposeMesh(const bool filterEmptyPatches);


        // Decomposed mesh and addressing

            //- Decompose mesh for a given processor
            autoPtr<fvMesh> processorMesh
            (
                const label procI,
                const Time& procDb
            ) const;

            //- Return processor point addressing
            const labelListList& procPointAddressing() const
            {
                return procPointAddressing_;
            }

            //- Return processor face addressing
            const labelListList& procFaceAddressing() const
            {
                return procFaceAddressing_;
            }

            //- Return processor cell addressing
            const labelListList& procCellAddressing() const
            {
                return procCellAddressing_;
            }

            //- Return processor boundary addressing
            const labelListList& procBoundaryAddressing() const
            {
                return procBoundaryAddressing_;
            }


        //- Write decomposition
        bool writeDecomposition();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
