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

void Foam::immersedBoundaryFvPatch::makeInvDirichletMatrices() const
{
    if (debug)
    {
        Info<< "immersedBoundaryFvPatch::makeInvDirichletMatrices() : "
            << "making immersed boundary inverse matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invDirichletMatricesPtr_)
    {
        FatalErrorIn
        (
            "void immersedBoundaryFvPatch::makeInvDirichletMatrices()"
        )   << "immersed boundary inverse least squares matrices already exist"
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();
    const vectorField& ibp = ibPoints();

    invDirichletMatricesPtr_ =
        new PtrList<scalarRectangularMatrix>(ibc.size());
    PtrList<scalarRectangularMatrix>& idm = *invDirichletMatricesPtr_;

    const vectorField& C = mesh_.C().internalField();

    scalarField conditionNumber(ibc.size(), 0.0);

    // Initialize maxRowSum for debug
    scalar maxRowSum = 0.0;

    const FieldField<Field, vector>& procC = ibProcCentres();

    label nCoeffs = 5;

    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    forAll (idm, cellI)
    {
        const labelList& interpCells = ibcc[cellI];
        const List<labelPair>& interpProcCells = ibcProcC[cellI];

        vectorField allPoints
        (
            interpCells.size()
          + interpProcCells.size(),
            vector::zero
        );

        if (allPoints.size() < nCoeffs)
        {
            FatalErrorIn
            (
                "void immersedBoundaryFvPatch::makeInvDirichletMatrices()"
            )   << "allPoints.size() < " << nCoeffs << " : "
                << allPoints.size() << abort(FatalError);
        }

        label pointID = 0;

        // Cells
        for (label i = 0; i < interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Processor cells
        for (label i = 0; i < interpProcCells.size(); i++)
        {
            allPoints[pointID++] =
                procC
                [
                    interpProcCells[i].first()
                ]
                [
                    interpProcCells[i].second()
                ];
        }

        // Weights calculation

        vector origin = C[ibc[cellI]];

        scalarField curDist = mag(allPoints - origin);

        // Calculate weights
        scalarField W =
            0.5*
            (
                1
              + cos(mathematicalConstant::pi*curDist/(1.1*max(curDist)))
            );

        idm.set
        (
            cellI,
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size(),
                0.0
            )
        );
        scalarRectangularMatrix& curMatrix = idm[cellI];

        scalarRectangularMatrix M
        (
            allPoints.size(),
            nCoeffs,
            0.0
        );

        origin = ibp[cellI];

        for(label i = 0; i < allPoints.size(); i++)
        {
            scalar X = allPoints[i].x() - origin.x();
            scalar Y = allPoints[i].y() - origin.y();

            label coeff = 0;
            M[i][coeff++] = X;
            M[i][coeff++] = Y;
            M[i][coeff++] = X*Y;
            M[i][coeff++] = sqr(X);
            M[i][coeff++] = sqr(Y);

            if (mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[i].z() - origin.z();
                M[i][coeff++] = Z;
                M[i][coeff++] = X*Z;
                M[i][coeff++] = Y*Z;
                M[i][coeff++] = sqr(Z);
            }
        }

        for (label i = 0; i < M.n(); i++)
        {
            for (label j = 0; j < M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        for (label i = 0; i < lsM.n(); i++)
        {
            for (label j = 0; j < lsM.m(); j++)
            {
                for (label k=0; k<M.n(); k++)
                {
                    lsM[i][j] += M[k][i]*M[k][j];
                }
            }
        }

        if (debug)
        {
            // Calculate matrix norm
            maxRowSum = 0.0;

            for (label i = 0; i < lsM.n(); i++)
            {
                scalar curRowSum = 0.0;

                for (label j = 0; j < lsM.m(); j++)
                {
                    curRowSum += lsM[i][j];
                }

                if (curRowSum > maxRowSum)
                {
                    maxRowSum = curRowSum;
                }
            }

            conditionNumber[cellI] = maxRowSum;
        }

        // Calculate inverse
        scalarSquareMatrix invLsM = lsM.LUinvert();

        for (label i = 0; i < lsM.n(); i++)
        {
            for (label j = 0; j < M.n(); j++)
            {
                for (label k = 0; k < lsM.n(); k++)
                {
                    curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                }
            }
        }

        if (debug)
        {
            // Calculate condition number
            maxRowSum = 0.0;

            for (label i = 0; i < lsM.n(); i++)
            {
                scalar curRowSum = 0.0;

                for (label j = 0; j < lsM.m(); j++)
                {
                    curRowSum += invLsM[i][j];
                }

                if (curRowSum > maxRowSum)
                {
                    maxRowSum = curRowSum;
                }
            }

            conditionNumber[cellI] *= maxRowSum;

            InfoIn
            (
                "void immersedBoundaryFvPatch::"
                "makeInvDirichletMatrices() const"
            )   << "Max Dirichlet matrix condition number: "
                << gMax(conditionNumber) << endl;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeInvNeumannMatrices() const
{
    if (debug)
    {
        Info<< "immersedBoundaryFvPatch::makeInvNeumannMatrices() : "
            << "making immersed boundary inverse least sqares matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invNeumannMatricesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeInvNeumannMatrices()")
            << "immersed boundary inverse least squares matrices already exist"
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();
    const vectorField& ibp = ibPoints();

    // Note: the algorithm is originally written with inward-facing normals
    // and subsequently changed: IB surface normals point outwards
    // HJ, 21/May/2012
    const vectorField& ibn = ibNormals();

    invNeumannMatricesPtr_ =
        new PtrList<scalarRectangularMatrix>(ibc.size());
    PtrList<scalarRectangularMatrix>& inm = *invNeumannMatricesPtr_;

    const vectorField& C = mesh_.C().internalField();

    scalarField conditionNumber(ibc.size(), 0);

    // Initialize maxRowSum for debug
    scalar maxRowSum = 0.0;

    const FieldField<Field, vector>& procC = ibProcCentres();

    label nCoeffs = 6;

    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    forAll (inm, cellI)
    {
        const labelList& interpCells = ibcc[cellI];
        const List<labelPair>& interpProcCells = ibcProcC[cellI];

        vectorField allPoints
        (
            interpCells.size() + 1
          + interpProcCells.size(),
            vector::zero
        );

        label pointID = 0;

        // Cells
        for (label i = 0; i < interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // IB point
        allPoints[pointID++] = ibp[cellI];

        // Processor cells
        for (label i = 0; i < interpProcCells.size(); i++)
        {
            allPoints[pointID++] =
                procC
                [
                    interpProcCells[i].first()
                ]
                [
                    interpProcCells[i].second()
                ];
        }

        // Weights calculation

        vector origin = C[ibc[cellI]];

        scalarField curR = mag(allPoints - origin);

        // Calculate weights
        scalarField W = 1 - curR/(1.1*max(curR));


        inm.set
        (
            cellI,
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size(),
                0.0
            )
        );
        scalarRectangularMatrix& curMatrix = inm[cellI];

        scalarRectangularMatrix M
        (
            allPoints.size(),
            nCoeffs,
            0.0
        );

        pointID = 0;
        origin = ibp[cellI];
        for (label i = 0; i < interpCells.size(); i++)
        {
            scalar X = allPoints[pointID].x() - origin.x();
            scalar Y = allPoints[pointID].y() - origin.y();

            M[pointID][0] = 1.0;
            M[pointID][1] = X;
            M[pointID][2] = Y;
            M[pointID][3] = X*Y;
            M[pointID][4] = sqr(X);
            M[pointID][5] = sqr(Y);

            if (mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[pointID].z() - origin.z();

                M[pointID][6] = Z;
                M[pointID][7] = X*Z;
                M[pointID][8] = Y*Z;
                M[pointID][9] = sqr(Z);
            }

            pointID++;
        }

        scalar X = allPoints[pointID].x() - origin.x();
        scalar Y = allPoints[pointID].y() - origin.y();

        M[pointID][0] = 0;
        M[pointID][1] = -ibn[cellI].x();
        M[pointID][2] = -ibn[cellI].y();
        M[pointID][3] =
        (
           -ibn[cellI].x()*Y
          - ibn[cellI].y()*X
        );
        M[pointID][4] = -2*ibn[cellI].x()*X;
        M[pointID][5] = -2*ibn[cellI].y()*Y;

        if (mesh_.nGeometricD() == 3)
        {
            scalar Z = allPoints[pointID].z() - origin.z();

            M[pointID][6] = -ibn[cellI].z();
            M[pointID][7] =
            (
               -ibn[cellI].x()*Z
              - ibn[cellI].z()*X
            );
            M[pointID][8] =
            (
               -ibn[cellI].y()*Z
              - ibn[cellI].z()*Y
            );
            M[pointID][9] = -2*ibn[cellI].z()*Z;
        }

        pointID++;

        for(label i = 0; i < interpProcCells.size(); i++)
        {
            scalar X = allPoints[pointID].x() - origin.x();
            scalar Y = allPoints[pointID].y() - origin.y();

            M[pointID][0] = 1.0;
            M[pointID][1] = X;
            M[pointID][2] = Y;
            M[pointID][3] = X*Y;
            M[pointID][4] = sqr(X);
            M[pointID][5] = sqr(Y);

            if (mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[pointID].z() - origin.z();

                M[pointID][6] = Z;
                M[pointID][7] = X*Z;
                M[pointID][8] = Y*Z;
                M[pointID][9] = sqr(Z);
            }

            pointID++;
        }

        for (label i = 0; i < M.n(); i++)
        {
            for (label j = 0; j < M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        for (label i = 0; i < lsM.n(); i++)
        {
            for (label j = 0; j < lsM.m(); j++)
            {
                for (label k = 0; k < M.n(); k++)
                {
                    lsM[i][j] += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        if (debug)
        {
            maxRowSum = 0.0;

            for (label i = 0; i < lsM.n(); i++)
            {
                scalar curRowSum = 0.0;

                for (label j = 0; j < lsM.m(); j++)
                {
                    curRowSum += lsM[i][j];
                }

                if (curRowSum > maxRowSum)
                {
                    maxRowSum = curRowSum;
                }
            }

            conditionNumber[cellI] = maxRowSum;
        }

        // Calculate inverse
        scalarSquareMatrix invLsM = lsM.LUinvert();

        for (label i = 0; i < lsM.n(); i++)
        {
            for (label j = 0; j < M.n(); j++)
            {
                for (label k = 0; k < lsM.n(); k++)
                {
                    curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                }
            }
        }

        // Calculate condition number
        if (debug)
        {
            maxRowSum = 0.0;

            for (label i = 0; i < lsM.n(); i++)
            {
                scalar curRowSum = 0.0;

                for (label j = 0; j < lsM.m(); j++)
                {
                    curRowSum += invLsM[i][j];
                }

                if (curRowSum > maxRowSum)
                {
                    maxRowSum = curRowSum;
                }
            }

            conditionNumber[cellI] *= maxRowSum;

            if (conditionNumber[cellI] > 1e6)
            {
                labelList CC = interpCells;
                sort(CC);

                Info<< "Condition = " << conditionNumber[cellI] << nl
                    << "cell cells: " << CC
                    << "M = " << lsM
                    << endl;
            }

            InfoIn
            (
                "void immersedBoundaryFvPatch::makeInvNeumannMatrices() const"
            )   << "Max Neumann matrix condition number: "
            << gMax(conditionNumber) << endl;
        }
    }
}


const Foam::PtrList<Foam::scalarRectangularMatrix>&
Foam::immersedBoundaryFvPatch::invDirichletMatrices() const
{
    if (!invDirichletMatricesPtr_)
    {
        makeInvDirichletMatrices();
    }

    return *invDirichletMatricesPtr_;
}


const Foam::PtrList<Foam::scalarRectangularMatrix>&
Foam::immersedBoundaryFvPatch::invNeumannMatrices() const
{
    if (!invNeumannMatricesPtr_)
    {
        makeInvNeumannMatrices();
    }

    return *invNeumannMatricesPtr_;
}


// ************************************************************************* //
