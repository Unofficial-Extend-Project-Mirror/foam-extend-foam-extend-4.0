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

\*----------------------------------------------------------------------------*/

#include "referredCellList.H"
#include "moleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

referredCellList::referredCellList(moleculeCloud& molCloud)
:
    List<referredCell >(),
    molCloud_(molCloud)
{}


referredCellList::referredCellList
(
    moleculeCloud& molCloud,
    const List<referredCell>& referredCells,
    const List<label>& realCells
)
:
    List< referredCell >(referredCells),
    molCloud_(molCloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

referredCellList::~referredCellList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void referredCellList::setRealCellsInRange()
{
    Info<< nl << "Finding real cells in range of referred cells" << endl;

    forAll(*this, rC)
    {
        const polyMesh& mesh(molCloud_.mesh());

        referredCell& refCell = (*this)[rC];

        DynamicList<label> realCellsFoundInRange;

        const vectorList& refCellPoints = refCell.vertexPositions();

        forAll(molCloud_.realFacesWithinRCutMaxOfAnyReferringPatch(), rCF)
        {
            const label f
            (
                molCloud_.realFacesWithinRCutMaxOfAnyReferringPatch()[rCF]
            );

            if (molCloud_.testPointFaceDistance(refCellPoints,f))
            {
                const label cellO(mesh.faceOwner()[f]);

                if (findIndex(realCellsFoundInRange, cellO) == -1)
                {
                    realCellsFoundInRange.append(cellO);
                }

                if (mesh.isInternalFace(f))
                {
                    // boundary faces will not have neighbour information

                    const label cellN(mesh.faceNeighbour()[f]);

                    if (findIndex(realCellsFoundInRange, cellN) == -1)
                    {
                        realCellsFoundInRange.append(cellN);
                    }
                }
            }
        }

        forAll(molCloud_.realPointsWithinRCutMaxOfAnyReferringPatch(), rCP)
        {
            const label p
            (
                molCloud_.realPointsWithinRCutMaxOfAnyReferringPatch()[rCP]
            );

            if (molCloud_.testPointFaceDistance(p,refCell))
            {
                const labelList& pCells(mesh.pointCells()[p]);

                forAll(pCells, pC)
                {
                    const label cellI(pCells[pC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }


        const edgeList& refCellEdges = refCell.edges();

        forAll(molCloud_.realEdgesWithinRCutMaxOfAnyReferringPatch(), rCE)
        {
            const label edgeIIndex
            (
                molCloud_.realEdgesWithinRCutMaxOfAnyReferringPatch()[rCE]
            );

            const edge& eI(mesh.edges()[edgeIIndex]);

            forAll(refCellEdges, rCE)
            {
                const edge& eJ(refCellEdges[rCE]);

                if
                (
                    molCloud_.testEdgeEdgeDistance
                    (
                        eI,
                        refCellPoints[eJ.start()],
                        refCellPoints[eJ.end()]
                    )
                )
                {
                    const labelList& eICells(mesh.edgeCells()[edgeIIndex]);

                    forAll(eICells, eIC)
                    {
                        const label cellI(eICells[eIC]);

                        if (findIndex(realCellsFoundInRange, cellI) == -1)
                        {
                            realCellsFoundInRange.append(cellI);
                        }
                    }
                }
            }
        }

//         scalar rCutMaxSqr = molCloud_.rCutMax()*molCloud_.rCutMax();
//
//         forAll (molCloud_.mesh().points(), pointIIndex)
//         {
//             const point& ptI
//             (
//                 molCloud_.mesh().points()[pointIIndex]
//             );
//
//             forAll(refCellPoints, rCP)
//             {
//                 if (magSqr(ptI - refCellPoints[rCP]) <= rCutMaxSqr)
//                 {
//                     const labelList& ptICells
//                     (
//                         molCloud_.mesh().pointCells()[pointIIndex]
//                     );
//
//                     forAll(ptICells, pIC)
//                     {
//                         const label cellI(ptICells[pIC]);
//
//                         if(findIndex(realCellsFoundInRange, cellI) == -1)
//                         {
//                             realCellsFoundInRange.append(cellI);
//                         }
//                     }
//                 }
//             }
//         }

        refCell.realCells() = realCellsFoundInRange.shrink();
    }
}


void referredCellList::referMolecules()
{
    // Create referred molecules for sending using cell occupancy and
    // cellSendingReferralLists

    const List<DynamicList<molecule*> >& cellOccupancy
    (
        molCloud_.cellOccupancy()
    );

    forAll(molCloud_.cellSendingReferralLists(), cSRL)
    {
        const sendingReferralList& sRL
        (
            molCloud_.cellSendingReferralLists()[cSRL]
        );

        List<DynamicList<referredMolecule> > molsToReferOut(sRL.size());

        forAll(sRL, sRLI)
        {
            List<molecule*> realMols = cellOccupancy[sRL[sRLI]];

            forAll (realMols, rM)
            {
                molecule* mol = realMols[rM];

                molsToReferOut[sRLI].append
                (
                    referredMolecule
                    (
                        mol->id(),
                        mol->position()
                    )
                );
            }

            molsToReferOut[sRLI].shrink();
        }

        // Send lists of referred molecules to other processors

        if (sRL.destinationProc() != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                OPstream toInteractingProc
                (
                    Pstream::blocking,
                    sRL.destinationProc()
                );

                toInteractingProc << molsToReferOut;
            }
        }
        else
        {
            // Refer molecules directly for referred cells on the same
            // processor.

            const receivingReferralList& rRL
            (
                molCloud_.cellReceivingReferralLists()[cSRL]
            );

            forAll(rRL, rRLI)
            {
                forAll(rRL[rRLI], rC)
                {
                    referredCell& refCellToRefMolsTo = (*this)[rRL[rRLI][rC]];

                    refCellToRefMolsTo.referInMols(molsToReferOut[rRLI]);
                }
            }
        }
    }

    // Receive referred molecule lists to and distribute to referredCells
    // according tocellReceivingReferralLists, referredCells deal with the
    // transformations themselves

    forAll(molCloud_.cellReceivingReferralLists(), cRRL)
    {
        const receivingReferralList& rRL
        (
            molCloud_.cellReceivingReferralLists()[cRRL]
        );

        List<List<referredMolecule> > molsToReferIn(rRL.size());

        if (rRL.sourceProc() != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                IPstream fromInteractingProc
                (
                    Pstream::blocking,
                    rRL.sourceProc()
                );

                fromInteractingProc >> molsToReferIn;
            }

            forAll(rRL, rRLI)
            {
                forAll(rRL[rRLI], rC)
                {
                    referredCell& refCellToRefMolsTo = (*this)[rRL[rRLI][rC]];

                    refCellToRefMolsTo.referInMols(molsToReferIn[rRLI]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
