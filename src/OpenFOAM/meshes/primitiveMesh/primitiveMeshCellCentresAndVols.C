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
    Efficient cell-centre calculation using face-addressing, face-centres and
    face-areas.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"

#include "tetPointRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellCentresAndVols() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Calculating cell centres and cell volumes"
            << endl;
    }

    // It is an error to attempt to recalculate cellCentres
    // if the pointer is already set
    if (cellCentresPtr_ || cellVolumesPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellCentresAndVols() const")
            << "Cell centres or cell volumes already calculated"
            << abort(FatalError);
    }

    // set the accumulated cell centre to zero vector
    cellCentresPtr_ = new vectorField(nCells());
    vectorField& cellCtrs = *cellCentresPtr_;

    // Initialise cell volumes to 0
    cellVolumesPtr_ = new scalarField(nCells());
    scalarField& cellVols = *cellVolumesPtr_;

    // Make centres and volumes
    makeCellCentresAndVols(faceCentres(), faceAreas(), cellCtrs, cellVols);

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellCentresAndVols() : "
            << "Finished calculating cell centres and cell volumes"
            << endl;
    }
}


void Foam::primitiveMesh::makeCellCentresAndVols
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    vectorField& cellCtrs,
    scalarField& cellVols
) const
{
    // Clear the fields for accumulation
    cellCtrs = vector::zero;
    cellVols = 0.0;

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // first estimate the approximate cell centre as the average of
    // face centres

    vectorField cEst(nCells(), vector::zero);
    labelField nCellFaces(nCells(), 0);

    forAll (own, facei)
    {
        cEst[own[facei]] += fCtrs[facei];
        nCellFaces[own[facei]] += 1;
    }

    forAll (nei, facei)
    {
        cEst[nei[facei]] += fCtrs[facei];
        nCellFaces[nei[facei]] += 1;
    }

    forAll(cEst, celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    const faceList& allFaces = faces();
    const pointField& allPoints = points();

    forAll(own, faceI)
    {
        const face& f = allFaces[faceI];

        if (f.size() == 3)
        {
            tetPointRef tpr
            (
                allPoints[f[2]],
                allPoints[f[1]],
                allPoints[f[0]],
                cEst[own[faceI]]
            );

            scalar tetVol = tpr.mag();

            // Accumulate volume-weighted tet centre
            cellCtrs[own[faceI]] += tetVol*tpr.centre();

            // Accumulate tet volume
            cellVols[own[faceI]] += tetVol;
        }
        else
        {
            forAll(f, pI)
            {
                tetPointRef tpr
                (
                    allPoints[f[pI]],
                    allPoints[f.prevLabel(pI)],
                    fCtrs[faceI],
                    cEst[own[faceI]]
                );

                scalar tetVol = tpr.mag();

                // Accumulate volume-weighted tet centre
                cellCtrs[own[faceI]] += tetVol*tpr.centre();

                // Accumulate tet volume
                cellVols[own[faceI]] += tetVol;
            }
        }
    }

    forAll(nei, faceI)
    {
        const face& f = allFaces[faceI];

        if (f.size() == 3)
        {
            tetPointRef tpr
            (
                allPoints[f[0]],
                allPoints[f[1]],
                allPoints[f[2]],
                cEst[nei[faceI]]
            );

            scalar tetVol = tpr.mag();

            // Accumulate volume-weighted tet centre
            cellCtrs[nei[faceI]] += tetVol*tpr.centre();

            // Accumulate tet volume
            cellVols[nei[faceI]] += tetVol;
        }
        else
        {
            forAll(f, pI)
            {
                tetPointRef tpr
                (
                    allPoints[f[pI]],
                    allPoints[f.nextLabel(pI)],
                    fCtrs[faceI],
                    cEst[nei[faceI]]
                );

                scalar tetVol = tpr.mag();

                // Accumulate volume-weighted tet centre
                cellCtrs[nei[faceI]] += tetVol*tpr.centre();

                // Accumulate tet volume
                cellVols[nei[faceI]] += tetVol;
            }
        }
    }

    cellCtrs /= cellVols + VSMALL;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::cellCentres() const
{
    if (!cellCentresPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellCentresPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        calcCellCentresAndVols();
    }

    return *cellVolumesPtr_;
}


// ************************************************************************* //
