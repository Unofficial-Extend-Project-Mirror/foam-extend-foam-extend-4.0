/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InClass
    decompositionMethod

\*---------------------------------------------------------------------------*/

#include "decompositionMethod.H"
#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionMethod, 0);
    defineRunTimeSelectionTable(decompositionMethod, dictionary);
    defineRunTimeSelectionTable(decompositionMethod, dictionaryMesh);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethod::calcCSR
(
    const labelListList& cellCells,
    List<int>& adjncy,
    List<int>& xadj
)
{
    // Count number of internal faces
    label nConnections = 0;

    forAll(cellCells, coarseI)
    {
        nConnections += cellCells[coarseI].size();
    }

    // Create the adjncy array as twice the size of the total number of
    // internal faces
    adjncy.setSize(nConnections);

    xadj.setSize(cellCells.size() + 1);


    // Fill in xadj
    // ~~~~~~~~~~~~
    label freeAdj = 0;

    forAll(cellCells, coarseI)
    {
        xadj[coarseI] = freeAdj;

        const labelList& cCells = cellCells[coarseI];

        forAll(cCells, i)
        {
            adjncy[freeAdj++] = cCells[i];
        }
    }
    xadj[cellCells.size()] = freeAdj;
}



void Foam::decompositionMethod::calcCSR
(
    const polyMesh& mesh,
    List<int>& adjncy,
    List<int>& xadj
)
{
    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    xadj.setSize(mesh.nCells() + 1);

    // Initialise the number of internal faces of the cells to twice the
    // number of internal faces
    label nInternalFaces = 2*mesh.nInternalFaces();

    // Check the boundary for coupled patches and add to the number of
    // internal faces
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<cyclicPolyPatch>(pbm[patchi]))
        {
            nInternalFaces += pbm[patchi].size();
        }
    }

    // Create the adjncy array the size of the total number of internal and
    // coupled faces
    adjncy.setSize(nInternalFaces);

    // Fill in xadj
    // ~~~~~~~~~~~~
    label freeAdj = 0;

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        const labelList& cFaces = mesh.cells()[cellI];

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if
            (
                mesh.isInternalFace(faceI)
             || isA<cyclicPolyPatch>(pbm[pbm.whichPatch(faceI)])
            )
            {
                freeAdj++;
            }
        }
    }
    xadj[mesh.nCells()] = freeAdj;


    // Fill in adjncy
    // ~~~~~~~~~~~~~~

    labelList nFacesPerCell(mesh.nCells(), 0);

    // Internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = mesh.faceOwner()[faceI];
        label nei = mesh.faceNeighbour()[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
    }

    // Coupled faces. Only cyclics done.
    forAll(pbm, patchi)
    {
        if (isA<cyclicPolyPatch>(pbm[patchi]))
        {
            const unallocLabelList& faceCells = pbm[patchi].faceCells();

            label sizeby2 = faceCells.size()/2;

            for (label facei=0; facei<sizeby2; facei++)
            {
                label own = faceCells[facei];
                label nei = faceCells[facei + sizeby2];

                adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
                adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
            }
        }
    }
}


void Foam::decompositionMethod::calcDistributedCSR
(
    const polyMesh& mesh,
    List<int>& adjncy,
    List<int>& xadj
)
{
    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalIndex globalCells(mesh.nCells());


    //
    // Make Metis Distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains cellCells (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    //


    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Get renumbered owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<int> globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh.nInternalFaces();

            forAll(pp, i)
            {
                globalNeighbour[bFaceI++] = globalCells.toGlobal
                (
                    faceOwner[faceI++]
                );
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour, false);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per cell
    List<int> nFacesPerCell(mesh.nCells(), 0);

    // Number of coupled faces
    label nCoupledFaces = 0;

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        nFacesPerCell[faceOwner[faceI]]++;
        nFacesPerCell[faceNeighbour[faceI]]++;
    }
    // Handle coupled faces
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces++;
                nFacesPerCell[faceOwner[faceI++]]++;
            }
        }
    }


    // Fill in xadj
    // ~~~~~~~~~~~~

    xadj.setSize(mesh.nCells() + 1);

    int freeAdj = 0;

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        freeAdj += nFacesPerCell[cellI];
    }
    xadj[mesh.nCells()] = freeAdj;



    // Fill in adjncy
    // ~~~~~~~~~~~~~~

    adjncy.setSize(2*mesh.nInternalFaces() + nCoupledFaces);

    nFacesPerCell = 0;

    // For internal faces is just offsetted owner and neighbour
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = globalCells.toGlobal(nei);
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = globalCells.toGlobal(own);
    }
    // For boundary faces is offsetted coupled neighbour
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                label own = faceOwner[faceI];
                adjncy[xadj[own] + nFacesPerCell[own]++] =
                    globalNeighbour[bFaceI];

                faceI++;
                bFaceI++;
            }
        }
    }
}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const label nCoarse,
    labelListList& cellCells
)
{
    if (fineToCoarse.size() != mesh.nCells())
    {
        FatalErrorIn
        (
            "decompositionMethod::calcCellCells"
            "(const labelList&, labelListList&) const"
        )   << "Only valid for mesh agglomeration." << exit(FatalError);
    }

    List<DynamicList<label> > dynCellCells(nCoarse);

    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = fineToCoarse[mesh.faceOwner()[faceI]];
        label nei = fineToCoarse[mesh.faceNeighbour()[faceI]];

        if (own != nei)
        {
            if (findIndex(dynCellCells[own], nei) == -1)
            {
                dynCellCells[own].append(nei);
            }
            if (findIndex(dynCellCells[nei], own) == -1)
            {
                dynCellCells[nei].append(own);
            }
        }
    }

    cellCells.setSize(dynCellCells.size());
    forAll(dynCellCells, coarseI)
    {
        cellCells[coarseI].transfer(dynCellCells[coarseI]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompositionDict
)
{
    word decompositionMethodTypeName(decompositionDict.lookup("method"));

    Info<< "Selecting decompositionMethod "
        << decompositionMethodTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(decompositionMethodTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "decompositionMethod::New"
            "(const dictionary& decompositionDict)"
        )   << "Unknown decompositionMethod "
            << decompositionMethodTypeName << endl << endl
            << "Valid decompositionMethods are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompositionDict));
}


Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
{
    word decompositionMethodTypeName(decompositionDict.lookup("method"));

    Info<< "Selecting decompositionMethod "
        << decompositionMethodTypeName << endl;

    dictionaryMeshConstructorTable::iterator cstrIter =
        dictionaryMeshConstructorTablePtr_->find(decompositionMethodTypeName);

    if (cstrIter == dictionaryMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "decompositionMethod::New"
            "(const dictionary& decompositionDict, "
            "const polyMesh& mesh)"
        )   << "Unknown decompositionMethod "
            << decompositionMethodTypeName << endl << endl
            << "Valid decompositionMethods are : " << endl
            << dictionaryMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompositionDict, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethod::decompose
(
    const pointField& points
)
{
    scalarField weights(points.size(), 1);

    return decompose(points, weights);
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints
)
{
    // Decompose based on agglomerated points
    labelList coarseDistribution(decompose(coarsePoints));

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    // Decompose based on agglomerated points
    labelList coarseDistribution(decompose(coarsePoints, coarseWeights));

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


// ************************************************************************* //
