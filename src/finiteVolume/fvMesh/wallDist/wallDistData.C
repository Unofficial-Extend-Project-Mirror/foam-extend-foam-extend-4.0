/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "wallDistData.H"
#include "patchDataWave.H"
#include "wallPolyPatch.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class TransferType>
Foam::wallDistData<TransferType>::wallDistData
(
    const Foam::fvMesh& mesh,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const bool correctWalls
)
:
    volScalarField
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    cellDistFuncs(mesh),
    field_(field),
    correctWalls_(correctWalls),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class TransferType>
Foam::wallDistData<TransferType>::~wallDistData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct for mesh geom/topo changes
template<class TransferType>
void Foam::wallDistData<TransferType>::correct()
{
    Info<< "wallDistData.correct() called" << endl;

    const polyMesh& mesh = cellDistFuncs::mesh();

    //
    // Fill data on wall patches with initial values
    //

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    labelHashSet wallPatchIDs(bMesh.size());
    forAll(bMesh, patchI)
    {
        if (bMesh[patchI].isWall())
        {
            wallPatchIDs.insert(patchI);
        }
    }

    // Get patchids of walls
    // labelHashSet wallPatchIDs(getPatchIDs<wallPolyPatch>());

    // Collect pointers to data on patches
    UPtrList<Field<Type> > patchData(mesh.boundaryMesh().size());

    forAll(field_.boundaryField(), patchI)
    {
        patchData.set(patchI, &field_.boundaryField()[patchI]);
    }

    // Do mesh wave
    patchDataWave<TransferType> wave
    (
        mesh,
        wallPatchIDs,
        patchData,
        correctWalls_
    );

    // Transfer cell values from wave into *this and field_
    transfer(wave.distance());

    field_.transfer(wave.cellData());

    // Transfer values on patches into boundaryField of *this and field_
    forAll(boundaryField(), patchI)
    {
        scalarField& waveFld = wave.patchDistance()[patchI];

        if (!isA<emptyFvPatchScalarField>(boundaryField()[patchI]))
        {
            boundaryField()[patchI].transfer(waveFld);

            Field<Type>& wavePatchData = wave.patchData()[patchI];

            field_.boundaryField()[patchI].transfer(wavePatchData);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();
}


// ************************************************************************* //
