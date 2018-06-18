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

#include "fieldBoundsRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(fieldBoundsRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    fieldBoundsRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldBoundsRefinement::fieldBoundsRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    fieldName_(coeffDict().lookup("fieldName")),
    lowerBound_(readScalar(coeffDict().lookup("lowerBound"))),
    upperBound_(readScalar(coeffDict().lookup("upperBound"))),
    cellPointCellSmoothing_
    (
        coeffDict().lookupOrDefault<Switch>("cellPointCellSmoothing", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::fieldBoundsRefinement::~fieldBoundsRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::fieldBoundsRefinement::refinementCellCandidates() const
{
    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);

    // Create temporary for the field (sharing the reference to volume field)
    tmp<volScalarField> tvf(vField);

    // Use cell-point-cell interpolation to smooth out the field
    if (cellPointCellSmoothing_)
    {
        // Create volume to point interpolation object
        const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

        // Interpolate from cell centres to points
        pointScalarField pField(vpi.interpolate(vField));

        // Create point to volume interpolation object
        const pointMesh& pMesh = pointMesh::New(mesh());
        const pointVolInterpolation pvi(pMesh, mesh());

        // Interpolate from points back to cell centres
        tmp<volScalarField> tInterpolatedVf = pvi.interpolate(pField);

        // Assign to the temporary object
        tvf = tInterpolatedVf;
    }

    // Get const reference to (possibly smoothed volume field)
    const volScalarField& vf = tvf();

    // Create storage for collection of cells. Assume that one in five cells
    // will be refined to prevent excessive resizing.
    dynamicLabelList refinementCandidates(Foam::max(100, mesh().nCells()/5));

    // Loop through internal field and collect cells to refine
    const scalarField& vfIn = vf.internalField();

    forAll (vfIn, cellI)
    {
        // Get current cell value
        const scalar& cellValue = vfIn[cellI];

        if (cellValue > lowerBound_ && cellValue < upperBound_)
        {
            // Cell value is within the bounds, append cell for potential
            // refinement
            refinementCandidates.append(cellI);
        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " cells as refinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::fieldBoundsRefinement::unrefinementPointCandidates() const
{
    // Get the field
    const volScalarField& vField =
        mesh().lookupObject<volScalarField>(fieldName_);

    // Create volume to point interpolation object
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh());

    // Interpolate the volume field from cell centres to points and get the
    // internal point field
    pointScalarField pField(vpi.interpolate(vField));
    const scalarField& pFieldIn = pField.internalField();

    // Create storage for collection of candidates. Assume that one in ten
    // mesh points will be unrefined to prevent excessive resizing
    dynamicLabelList unrefinementCandidates
    (
        Foam::max(100, mesh().nPoints()/10)
    );

    // Loop through all split points and select candidates to unrefine
    forAll (pField, pointI)
    {
        // Get point value
        const scalar pointValue = pFieldIn[pointI];

        if (pointValue > upperBound_ || pointValue < lowerBound_)
        {
            // Point value is outside of bounds, append point for potential
            // unrefinement
            unrefinementCandidates.append(pointI);
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
