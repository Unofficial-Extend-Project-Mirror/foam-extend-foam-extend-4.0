/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "immersedBoundaryFvPatch.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeIbSamplingWeights() const
{
    if (debug)
    {
        Info<< "immersedBoundaryFvPatch::makeIbSamplingWeights() : "
            << "making sampling point weights"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibSamplingWeightsPtr_ || ibSamplingProcWeightsPtr_)
    {
        FatalErrorIn("void immersedBoundaryFvPatch::makeIbSamplingWeights()")
            << "sampling point weights already exist"
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();

    // Initialise the weights
    ibSamplingWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& cellWeights = *ibSamplingWeightsPtr_;

    forAll (cellWeights, cellI)
    {
        cellWeights[cellI].setSize(ibcc[cellI].size(), 0);
    }

    ibSamplingProcWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& cellProcWeights = *ibSamplingProcWeightsPtr_;

    forAll (cellProcWeights, cellI)
    {
        cellProcWeights[cellI].setSize(ibcProcC[cellI].size(), 0);
    }

    // Get sampling point locations
    const vectorField& samplingPoints = ibSamplingPoints();
    const scalarField& gammaIn = gamma().internalField();
    const vectorField& CIn = mesh_.C().internalField();

    const FieldField<Field, scalar>& gammaProc = ibProcGamma();
    const FieldField<Field, vector>& CProc = ibProcCentres();

    // Go through all cellCells and calculate inverse distance for
    // all live points
    forAll (samplingPoints, cellI)
    {
        const vector& curP = samplingPoints[cellI];

        scalar sumW = 0;

        // Local weights
        scalarList& curCW = cellWeights[cellI];

        const labelList& curCells = ibcc[cellI];

        forAll (curCells, ccI)
        {
            // Only pick live cells
            if (gammaIn[curCells[ccI]] > SMALL)
            {
                curCW[ccI] = 1/mag(CIn[curCells[ccI]] - curP);
                sumW += curCW[ccI];
            }
            else
            {
                curCW[ccI] = 0;
            }
        }

        // Processor weights
        const List<labelPair>& interpProcCells = ibcProcC[cellI];

        scalarList& curProcCW = cellProcWeights[cellI];

        forAll (interpProcCells, cProcI)
        {
            if
            (
                gammaProc
                [
                    interpProcCells[cProcI].first()
                ]
                [
                    interpProcCells[cProcI].second()
                ] > SMALL
            )
            {
                curProcCW[cProcI] =
                    1/mag
                    (
                        CProc
                        [
                            interpProcCells[cProcI].first()
                        ]
                        [
                            interpProcCells[cProcI].second()
                        ] - curP
                    );

                sumW += curProcCW[cProcI];
            }
            else
            {
                curProcCW[cProcI] = 0;
            }
        }

        // Divide through by the sum
        if (sumW < SMALL)
        {
            InfoIn
            (
                "void immersedBoundaryFvPatch::makeIbSamplingWeights()"
            )   << "Insufficient live neighbourhood for IB cell "
                << ibc[cellI] << "." << nl
                << "Please adjust radiusFactor, angleFactor or "
                << "immersedBoundaryMaxCellCellRows "
                << "in immersedBoundaryFvPatch."
                << endl;

            // Reset sum and weights and use all points
             sumW = 0;
             curCW = 0;

             forAll (curCells, ccI)
             {
                 // Use all cells
                 curCW[ccI] = 1/mag(CIn[curCells[ccI]] - curP);
                 sumW += curCW[ccI];
             }
        }

        forAll (curCells, ccI)
        {
            curCW[ccI] /= sumW;
        }

        forAll (curProcCW, cProcI)
        {
            curProcCW[cProcI] /= sumW;
        }
    }
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibSamplingWeights() const
{
    if (!ibSamplingWeightsPtr_)
    {
        makeIbSamplingWeights();
    }

    return *ibSamplingWeightsPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibSamplingProcWeights() const
{
    if (!ibSamplingProcWeightsPtr_)
    {
        makeIbSamplingWeights();
    }

    return *ibSamplingProcWeightsPtr_;
}


// ************************************************************************* //
