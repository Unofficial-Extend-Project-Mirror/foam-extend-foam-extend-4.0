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
    Foam::decompositionMethod

Description
    Abstract base class for decomposition

SourceFiles
    decompositionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionMethod_H
#define decompositionMethod_H

#include "polyMesh.H"
#include "pointField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class decompositionMethod Declaration
\*---------------------------------------------------------------------------*/

class decompositionMethod
{
protected:

    // Protected data

        //- Decomposition dictionary
        const dictionary& decompositionDict_;

        //- Number of processors in decomposition
        label nProcessors_;


        //- Helper to convert connectivity supplied as cellCells into
        //  simple CSR (Metis, scotch) storage
        static void calcCSR
        (
            const labelListList& globalCellCells,
            List<int>& adjncy,
            List<int>& xadj
        );

        //- Helper: convert local connectivity from the mesh
        //  into CSR (Metis, scotch) storage
        //  Treats cyclics as coupled, but not processor patches
        static void calcCSR
        (
            const polyMesh& mesh,
            List<int>& adjncy,
            List<int>& xadj
        );

        //- Helper: convert mesh connectivity into distributed CSR
        //  Very dubious coding.  HJ, 1/Mar/2011
        static void calcDistributedCSR
        (
            const polyMesh&,
            List<int>& adjncy,
            List<int>& xadj
        );

        //- Helper: determine (non-parallel) cellCells from mesh
        // with agglomeration
        static void calcCellCells
        (
            const polyMesh& mesh,
            const labelList& fineToCoarse,
            const label nCoarse,
            labelListList& cellCells
        );

        static void fixCyclics
        (
            const polyMesh& mesh,
            labelList& decomp
        );


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        decompositionMethod(const decompositionMethod&);

        //- Disallow default bitwise assignment
        void operator=(const decompositionMethod&);

        //- load missing libraries to avoid compile-time linking to external dependencies
        static void loadExternalLibraries();

public:

    //- Runtime type information
    TypeName("decompositionMethod");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionMethod,
            dictionary,
            (
                const dictionary& decompositionDict
            ),
            (decompositionDict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionMethod,
            dictionaryMesh,
            (
                const dictionary& decompositionDict,
                const polyMesh& mesh
            ),
            (decompositionDict, mesh)
        );


    // Selectors

        //- Return a pointer to the selected decomposition method
        static autoPtr<decompositionMethod> New
        (
            const dictionary& decompositionDict
        );

        //- Return a reference to the selected decomposition method
        static autoPtr<decompositionMethod> New
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Constructors

        //- Construct given the decomposition dictionary
        decompositionMethod(const dictionary& decompositionDict)
        :
            decompositionDict_(decompositionDict),
            nProcessors_
            (
                readLabel(decompositionDict.lookup("numberOfSubdomains"))
            )
        {}


    // Destructor

        virtual ~decompositionMethod()
        {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const = 0;

        //- Decompose cells
        //  If needed, use connectivity directly from the mesh
        //  Calls decompose (below) with uniform weights
        virtual labelList decompose(const pointField&);

        //- Decompose cells with weights
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        ) = 0;


        //- Decompose cell clusters
        //  Calls decompose (below) with uniform weights
        virtual labelList decompose
        (
            const labelList& fineToCoarse,
            const pointField& coarsePoints
        );

        //- Decompose cell clusters with weights on clusters
        virtual labelList decompose
        (
            const labelList& fineToCoarse,
            const pointField& coarsePoints,
            const scalarField& coarseWeights
        );


        //- Decompose cells with weights with explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
