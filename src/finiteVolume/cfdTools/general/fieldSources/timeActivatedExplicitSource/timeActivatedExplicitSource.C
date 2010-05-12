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

\*---------------------------------------------------------------------------*/

#include "timeActivatedExplicitSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedExplicitSource::timeActivatedExplicitSource
(
    const word& sourceName,
    const fvMesh& mesh
)
:
    dict_
    (
        IOobject
        (
            sourceName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    runTime_(mesh.time()),
    cellSource_(dict_.lookup("cellSource")),
    timeStart_(dimensionedScalar(dict_.lookup("timeStart")).value()),
    duration_(dimensionedScalar(dict_.lookup("duration")).value()),
    onValue_(dict_.lookup("onValue")),
    offValue_(dict_.lookup("offValue")),
    currentValue_(dimensionedScalar("zero", onValue_.dimensions(), 0.0)),
    cellSelector_
    (
        topoSetSource::New
        (
            cellSource_,
            mesh,
            dict_.subDict(cellSource_ + "Coeffs")
        )
    ),
    selectedCellSet_
    (
        mesh,
        "timeActivatedExplicitSourceCellSet",
        mesh.nCells()/10 + 1  // Reasonable size estimate.
    )
{
    // Check dimensions of on/off values are consistent
    if (onValue_.dimensions() != offValue_.dimensions())
    {
        FatalErrorIn
        (
            "Foam::timeActivatedExplicitSource::timeActivatedExplicitSource"
        )<< "Dimensions of on and off values must be equal" << nl
         << "onValue = " << onValue_ << nl
         << "offValue = " << offValue_ << exit(FatalError);
    }

    // Create the cell set
    cellSelector_->applyToSet
    (
        topoSetSource::NEW,
        selectedCellSet_
    );

    // Give some feedback
    Info<< "timeVaryingExplitSource(" << sourceName << ")" << nl
        << "Selected " << returnReduce(selectedCellSet_.size(), sumOp<label>())
        << " cells." << endl;

    // Initialise the value
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::timeActivatedExplicitSource::timeStart() const
{
    return timeStart_;
}


Foam::scalar Foam::timeActivatedExplicitSource::duration() const
{
    return duration_;
}


const Foam::dimensionedScalar&
Foam::timeActivatedExplicitSource::currentValue() const
{
    return currentValue_;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::timeActivatedExplicitSource::Su() const
{
    tmp<DimensionedField<scalar, volMesh> > tSource
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "timeActivatedExplicitSource",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", onValue_.dimensions(), 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& sourceField = tSource();

    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        sourceField[iter.key()] = currentValue_.value();
    }

    return tSource;
}


void Foam::timeActivatedExplicitSource::update()
{
    if
    (
        (runTime_.time().value() >= timeStart_)
     && (runTime_.time().value() <= timeStart_ + duration_)
    )
    {
        currentValue_ = onValue_;
    }
    else
    {
        currentValue_ = offValue_;
    }
}


// ************************************************************************* //
