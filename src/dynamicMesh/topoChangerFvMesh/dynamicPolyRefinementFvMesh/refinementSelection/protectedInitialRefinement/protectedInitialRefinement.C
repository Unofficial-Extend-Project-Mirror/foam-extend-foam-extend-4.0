/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "protectedInitialRefinement.H"
#include "labelIOField.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(protectedInitialRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    protectedInitialRefinement,
    dictionary
);

}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::protectedInitialRefinement::cellLevelAsField() const
{
    // Get cell level
    const labelIOField& cLevel = mesh().lookupObject<labelIOField>("cellLevel");

    // Create cell level field as volScalarField
    tmp<volScalarField> tCellLevelField
    (
        new volScalarField
        (
            IOobject
            (
                "cellLevel",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& cellLevelField = tCellLevelField();
    scalarField& cellLevelFieldIn = cellLevelField.internalField();

    // Set the field
    forAll (cLevel, cellI)
    {
        cellLevelFieldIn[cellI] = cLevel[cellI];
    }
    cellLevelField.correctBoundaryConditions();

    return tCellLevelField;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::protectedInitialRefinement::protectedInitialRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    initialCellLevel_
    (
        IOobject
        (
            "initialCellLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        cellLevelAsField()()
    ),
    // Note: conversion from label to scalar
    minProtectedLevel_
    (
        coeffDict().lookupOrDefault<label>("minProtectedLevel", 0)
    )
{
    Info<< "Will not allow initial refinement below level: "
        << minProtectedLevel_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::protectedInitialRefinement::~protectedInitialRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::protectedInitialRefinement::refinementCellCandidates() const
{
    // All cells are candidates for refinement
    labelList refinementCandidates(mesh().nCells());
    forAll (refinementCandidates, cellI)
    {
        refinementCandidates[cellI] = cellI;
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " (all) cells as refinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::protectedInitialRefinement::unrefinementPointCandidates() const
{
    // Create mark-up field for points that are found on initially refined cells
    // Note: this can be done more intelligently by either cacheing until load
    // balancing actually occurs or holding pointLevel as data member and
    // updating it on load balancing. The overall overhead compared to AMR and
    // DLB is really small, so there's no need to worry about this yet.
    // VV, 10/July/2018
    boolList pointsOnRefinedCells(mesh().nPoints(), false);

    // Get initial cell level and mesh data
    const scalarField& initialCellLevelIn = initialCellLevel_.internalField();
    const labelListList& meshCellPoints = mesh().cellPoints();

    // Loop through all cells and count number of protected points
    label nProtectedPoints = 0;

    // Get current cell level
    const labelIOField& curCellLevel =
        mesh().lookupObject<labelIOField>("cellLevel");

    forAll (initialCellLevelIn, cellI)
    {
        const scalar& cl = initialCellLevelIn[cellI];
        const scalar& curcl = curCellLevel[cellI];

        if
        (
            (cl > SMALL)
         && (curcl - minProtectedLevel_ < 0.5)
        )
        {
            // Cell has been refined during meshing and the original level is
            // equal to the current level, mark all of its points
            const labelList& curCellPoints = meshCellPoints[cellI];

            forAll (curCellPoints, i)
            {
                // Get marker for this point
                bool& ptOnRefCell = pointsOnRefinedCells[curCellPoints[i]];

                if (!ptOnRefCell)
                {
                    // This points has not been marked yet, mark it and
                    // increment the counter for protected points
                    ptOnRefCell = true;
                    ++nProtectedPoints;
                }
            }
        }
    }

    // Create the list for unrefinement candidates
    labelList unrefinementCandidates(mesh().nPoints() - nProtectedPoints);
    label nUnrefPoints = 0;

    forAll (pointsOnRefinedCells, pointI)
    {
        if (!pointsOnRefinedCells[pointI])
        {
            // This point is an unrefinement candidate, set it and increment
            unrefinementCandidates[nUnrefPoints++] = pointI;
        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //
