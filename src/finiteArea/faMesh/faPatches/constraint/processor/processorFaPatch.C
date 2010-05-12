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

Description

\*---------------------------------------------------------------------------*/

#include "processorFaPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "IPstream.H"
#include "OPstream.H"
#include "transformField.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorFaPatch, 0);
addToRunTimeSelectionTable(faPatch, processorFaPatch, dictionary);


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorFaPatch::~processorFaPatch()
{
    deleteDemandDrivenData(neighbPointsPtr_);
    deleteDemandDrivenData(nonGlobalPatchPointsPtr_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void processorFaPatch::makeNonGlobalPatchPoints() const
{
    // If it is not runing parallel or there are no global points
    // create a 1->1 map

    // Can not use faGlobalMeshData at this point yet
    
    if
    (
        !Pstream::parRun()
        || !boundaryMesh().mesh()().globalData().nGlobalPoints()
//         || !boundaryMesh().mesh().globalData().nGlobalPoints()
    )
    {
        nonGlobalPatchPointsPtr_ = new labelList(nPoints());
        labelList& ngpp = *nonGlobalPatchPointsPtr_;
        forAll(ngpp, i)
        {
            ngpp[i] = i;
        }
    }
    else
    {

        // Get reference to shared points
        const labelList& sharedPoints =
            boundaryMesh().mesh()().globalData().sharedPointLabels();

        nonGlobalPatchPointsPtr_ = new labelList(nPoints());
        labelList& ngpp = *nonGlobalPatchPointsPtr_;

        const labelList& faMeshPatchPoints = pointLabels();

        const labelList& meshPoints = 
            boundaryMesh().mesh().patch().meshPoints();

        label noFiltPoints = 0;

        forAll (faMeshPatchPoints, pointI)
        {
            label curP = meshPoints[faMeshPatchPoints[pointI]];

            bool found = false;

            forAll (sharedPoints, sharedI)
            {
                if (sharedPoints[sharedI] == curP)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                ngpp[noFiltPoints] = pointI;
                noFiltPoints++;
            }
        }

        ngpp.setSize(noFiltPoints);


//         // Get reference to shared points
//         const labelList& sharedPoints =
//             boundaryMesh().mesh().globalData().sharedPointLabels();

//         nonGlobalPatchPointsPtr_ = new labelList(nPoints());
//         labelList& ngpp = *nonGlobalPatchPointsPtr_;

//         const labelList& patchPoints = pointLabels();

//         label noFiltPoints = 0;

//         forAll (patchPoints, pointI)
//         {
//             label curP = patchPoints[pointI];

//             bool found = false;

//             forAll (sharedPoints, pI)
//             {
//                 if (sharedPoints[pI] == curP)
//                 {
//                     found = true;
//                     break;
//                 }
//             }

//             if (!found)
//             {
//                 ngpp[noFiltPoints] = pointI;
//                 noFiltPoints++;
//             }
//         }

//         ngpp.setSize(noFiltPoints);
    }
}


void processorFaPatch::initGeometry()
{
    if (Pstream::parRun())
    {
        OPstream toNeighbProc
        (
            Pstream::blocking,
            neighbProcNo(), 
            3*(sizeof(label) + size()*sizeof(vector))
        );

        toNeighbProc
            << edgeCentres()
            << edgeLengths()
            << edgeFaceCentres();
    }
}


void processorFaPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc
            (
                Pstream::blocking,
                neighbProcNo(), 
                3*(sizeof(label) + size()*sizeof(vector))
            );
            fromNeighbProc
                >> neighbEdgeCentres_
                >> neighbEdgeLengths_
                >> neighbEdgeFaceCentres_;
        }

        const scalarField& magEl = magEdgeLengths();

        forAll(magEl, edgei)
        {
            scalar nmagEl = mag(neighbEdgeLengths_[edgei]);
            scalar avEl = (magEl[edgei] + nmagEl)/2.0;

            if (mag(magEl[edgei] - nmagEl)/avEl > 1e-6)
            {
                FatalErrorIn
                (
                    "processorFvPatch::makeWeights(scalarField& w) const"
                )   << "edge " << edgei
                    << " length does not match neighbour by "
                    << 100*mag(magEl[edgei] - nmagEl)/avEl
                    << "% -- possible edge ordering problem"
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            edgeCentres(),
            neighbEdgeCentres_,
            edgeNormals(),
            neighbEdgeLengths_/mag(neighbEdgeLengths_)
        );
    }
}


void processorFaPatch::initMovePoints(const pointField& p)
{
    faPatch::movePoints(p);
    initGeometry();
}


void processorFaPatch::movePoints(const pointField&)
{
    calcGeometry();
}


void processorFaPatch::initUpdateMesh()
{
    // For completeness
    faPatch::initUpdateMesh();

    deleteDemandDrivenData(neighbPointsPtr_);

    if (Pstream::parRun())
    {
        // Express all points as patch edge and index in edge.
        labelList patchEdge(nPoints());
        labelList indexInEdge(nPoints());

        const edgeList::subList patchEdges = 
            patchSlice(boundaryMesh().mesh().edges());

        const labelListList& ptEdges = pointEdges();

        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
        {
            label edgeI = ptEdges[patchPointI][0];

            patchEdge[patchPointI] = edgeI;

            const edge& e = patchEdges[edgeI];

            indexInEdge[patchPointI] = 
                findIndex
                (
                    e,
                    pointLabels()[patchPointI]
                );
        }

        OPstream toNeighbProc
        (
            Pstream::blocking,
            neighbProcNo(),
            2*sizeof(label) + 2*nPoints()*sizeof(label)
        );

        toNeighbProc
            << patchEdge
            << indexInEdge;
    }
}


void processorFaPatch::updateMesh()
{
    // For completeness
    faPatch::updateMesh();

    if (Pstream::parRun())
    {
        labelList nbrPatchEdge(nPoints());
        labelList nbrIndexInEdge(nPoints());

        {
            // Note cannot predict exact size since edgeList not (yet) sent as
            // binary entity but as List of edges.
            IPstream fromNeighbProc
            (
                Pstream::blocking,
                neighbProcNo()
            );

            fromNeighbProc
                >> nbrPatchEdge
                >> nbrIndexInEdge;
        }

        if (nbrPatchEdge.size() == nPoints())
        {
            // Convert neighbour edges and indices into face back into
            // my edges and points.
            neighbPointsPtr_ = new labelList(nPoints());
            labelList& neighbPoints = *neighbPointsPtr_;

            const edgeList::subList patchEdges = 
                patchSlice(boundaryMesh().mesh().edges());

            forAll(nbrPatchEdge, nbrPointI)
            {
                // Find edge and index in edge on this side.
                const edge& e = patchEdges[nbrPatchEdge[nbrPointI]];

                label index =  1 - nbrIndexInEdge[nbrPointI];

                label patchPointI = findIndex(pointLabels(), e[index]);

                neighbPoints[patchPointI] = nbrPointI;
            }
        }
        else
        {
            // Differing number of points. Probably patch includes
            // part of a cyclic.
            neighbPointsPtr_ = NULL;
        }
    }
}


const labelList& processorFaPatch::neighbPoints() const
{
    if (!neighbPointsPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorIn("processorFaPatch::neighbPoints() const")
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points  on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }
 
   return *neighbPointsPtr_;
}


// Make patch weighting factors
void processorFaPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        // The face normals point in the opposite direction on the other side
        scalarField neighbEdgeCentresCn
        (
            (
                neighbEdgeLengths()
               /mag(neighbEdgeLengths())
            )
          & (
              neighbEdgeCentres()
            - neighbEdgeFaceCentres())
        );

        w = neighbEdgeCentresCn/
            (
                (edgeNormals() & faPatch::delta())
              + neighbEdgeCentresCn
            );
    }
    else
    {
        w = 1.0;
    }
}


// Make patch edge - neighbour face distances
void processorFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (Pstream::parRun())
    {
        dc = (1.0 - weights())/(edgeNormals() & faPatch::delta());
    }
    else
    {
        dc = 1.0/(edgeNormals() & faPatch::delta());
    }
}


// Return delta (P to N) vectors across coupled patch
tmp<vectorField> processorFaPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                faPatch::delta()
              - (
                    neighbEdgeCentres()
                  - neighbEdgeFaceCentres()
                );
        }
        else
        {
            return
                faPatch::delta()
              - transform
                (
                    forwardT(),
                    (
                        neighbEdgeCentres()
                      - neighbEdgeFaceCentres()
                    )
                );
        }
    }
    else
    {
        return faPatch::delta();
    }
}


const labelList& processorFaPatch::nonGlobalPatchPoints() const
{
    if (!nonGlobalPatchPointsPtr_)
    {
        makeNonGlobalPatchPoints();
    }

    return *nonGlobalPatchPointsPtr_;
}

tmp<labelField> processorFaPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


void processorFaPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    send(commsType, interfaceData);
}


tmp<labelField> processorFaPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


void processorFaPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


tmp<labelField> processorFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


// Write
void processorFaPatch::write(Ostream& os) const
{
    faPatch::write(os);
    os.writeKeyword("myProcNo") << myProcNo_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbProcNo") << neighbProcNo_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
