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

\*---------------------------------------------------------------------------*/

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"
#include "scotchDecomp.H"

extern "C"
{
#define OMPI_SKIP_MPICXX
#   include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& cWeights,
    List<int>& finalDecomp
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("k-way");

    int numCells = xadj.size()-1;

    // decomposition options. 0 = use defaults
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // processor weights initialised with no size, only used if specified in
    // a file - use the data type that metis expects here
    Field<real_t> processorWeights;

    // cell weights (so on the vertices of the dual)
    List<int> cellWeights;

    // face weights (so on the edges of the dual)
    List<int> faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
                << exit(FatalError);
        }
        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }


    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");
        word weightsFile;

        if (metisCoeffs.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Method " << method
                    << " in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis method     " << method
                << nl << endl;
        }

        List<int> mOptions;
        if (metisCoeffs.readIfPresent("options", mOptions))
        {
            if (mOptions.size() != METIS_NOPTIONS)
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            forAll(mOptions, i)
            {
                options[i] = mOptions[i];
            }

            Info<< "metisDecomp : Using Metis options     " << mOptions
                << nl << endl;
        }

        if (metisCoeffs.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != nProcessors_)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nProcessors_
                    << exit(FatalError);
            }
        }

        if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "metisDecomp : Using cell-based weights." << endl;

            IOList<int> cellIOWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
            cellWeights.transfer(cellIOWeights);

            if (cellWeights.size() != xadj.size()-1)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " does not equal number of cells " << xadj.size()-1
                    << exit(FatalError);
            }
        }
    }

    int nProcs = nProcessors_;

    // output: cell -> processor addressing
    finalDecomp.setSize(numCells);

    // output: number of cut edges
    int edgeCut = 0;

    // Vertex weight info
    int* vwgtPtr = NULL;
    int* adjwgtPtr = NULL;

    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
    }
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
    }

    int one = 1;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            NULL,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            NULL,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            NULL,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            NULL,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose(const pointField&,const scalarField&)"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    List<int> adjncy;
    List<int> xadj;
    calcCSR
    (
        mesh_,
        adjncy,
        xadj
    );

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


Foam::labelList Foam::metisDecomp::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    if (fineToCoarse.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose"
            "(const labelList&, const pointField&, const scalarField&)"
        )   << "Size of cell-to-coarse map " << fineToCoarse.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    List<int> adjncy;
    List<int> xadj;
    {
        // Get cellCells on coarse mesh.
        labelListList cellCells;
        calcCellCells
        (
            mesh_,
            fineToCoarse,
            coarsePoints.size(),
            cellCells
        );

        calcCSR(cellCells, adjncy, xadj);
    }

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, coarseWeights, finalDecomp);


    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[fineToCoarse[i]];
    }

    fixCyclics(mesh_, fineDistribution);

    return fineDistribution;
}


Foam::labelList Foam::metisDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights
)
{
    if (cc.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose"
            "(const pointField&, const labelListList&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cc.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> adjncy;
    List<int> xadj;
    calcCSR(globalCellCells, adjncy, xadj);

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, cWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


// ************************************************************************* //
