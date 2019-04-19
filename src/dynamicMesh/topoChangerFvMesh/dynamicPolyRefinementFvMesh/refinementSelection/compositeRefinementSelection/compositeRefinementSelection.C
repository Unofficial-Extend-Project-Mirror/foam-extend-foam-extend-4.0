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

#include "compositeRefinementSelection.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(compositeRefinementSelection, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    compositeRefinementSelection,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compositeRefinementSelection::compositeRefinementSelection
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    baseRefinementSelections_()
{
    // Read basic refinement selections
    PtrList<entry> baseRefSelectionEntries
    (
        coeffDict().lookup("baseRefinementSelections")
    );

    baseRefinementSelections_.setSize(baseRefSelectionEntries.size());

    forAll (baseRefinementSelections_, brsI)
    {
        baseRefinementSelections_.set
        (
            brsI,
            refinementSelection::New
            (
                mesh,
                baseRefSelectionEntries[brsI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::compositeRefinementSelection::~compositeRefinementSelection()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::compositeRefinementSelection::refinementCellCandidates() const
{
    // Final refinement cell candidates are defined as the intersection of all
    // sets (obtained with different refinement selection algorithms)
    // In order to define the intersection in a straightforward and efficient
    // way, we will create a labelField for all cells counting the number of
    // times this cell has been selected for refinement. If all selection
    // algorithms have selected the cell, this cell is marked as a final
    // refinement candidate.

    // Create field counting the number of selections
    labelField nSelections(mesh().nCells(), 0);

    // Loop through all base refinement selections
    forAll (baseRefinementSelections_, brsI)
    {
        // Get refinement candidates from this base selection algorithm. Note:
        // list is transferred
        const labelList curRefCandidates
        (
            baseRefinementSelections_[brsI].refinementCellCandidates()
        );

        // Increment the number of selections for selected cells
        forAll (curRefCandidates, i)
        {
            ++nSelections[curRefCandidates[i]];
        }
    }

    // Create storage for collection of final cell candidates. Assume that
    // one fifth of the cells will be marked to prevent excessive resizing
    dynamicLabelList refinementCandidates(mesh().nCells()/5);

    // Get number of active selection algorithms
    const label nBaseSelections = baseRefinementSelections_.size();

    // Loop through all cells and collect final refinement candidates
    forAll (nSelections, cellI)
    {
        if (nSelections[cellI] == nBaseSelections)
        {
            // Cell has been marked by all selection algorithms, append it
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
Foam::compositeRefinementSelection::unrefinementPointCandidates() const
{
    // Final unrefinement point candidates are defined as the intersection of
    // all sets (obtained with different refinement selection algorithms) In
    // order to define the intersection in a straightforward and efficient way,
    // we will create a labelField for all points counting the number of times
    // this point has been selected for unrefinement. If all selection
    // algorithms have selected the point, this point is marked as a final
    // unrefinement split point candidate.

    // Create a field counting the number of selections
    labelField nSelections(mesh().nPoints(), 0);

    // Loop through all base refinement selections
    forAll (baseRefinementSelections_, brsI)
    {
        // Get unrefinement candidates from this base selection algorithm. Note:
        // list is transferred
        const labelList curUnrefCandidates
        (
            baseRefinementSelections_[brsI].unrefinementPointCandidates()
        );

        // Increment the number of selections for selected points
        forAll (curUnrefCandidates, i)
        {
            ++nSelections[curUnrefCandidates[i]];
        }
    }

    // Create storage for collection of final point candidates. Assume that one
    // tenth of the points will be marked to prevent excessive resizing
    dynamicLabelList unrefinementCandidates(mesh().nPoints()/10);

    // Get number of active selection algorithms
    const label nBaseSelections = baseRefinementSelections_.size();

    // Loop through all points and collect final unrefinement candidates
    forAll (nSelections, pointI)
    {
        if (nSelections[pointI] == nBaseSelections)
        {
            // Point has been marked by all selection algorithms, append it
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
