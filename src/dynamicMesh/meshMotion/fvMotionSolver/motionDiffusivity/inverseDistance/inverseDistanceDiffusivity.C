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

#include "inverseDistanceDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "patchWave.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::inverseDistanceDiffusivity
(
    const fvMotionSolver& mSolver,
    Istream& mdData
)
:
    uniformDiffusivity(mSolver, mdData),
    patchNames_(mdData)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::~inverseDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::inverseDistanceDiffusivity::y() const
{
    const polyMesh& mesh = mSolver().mesh();

    labelHashSet patchSet(mesh.boundaryMesh().patchSet(patchNames_));

    if (patchSet.size())
    {
        return tmp<scalarField>
        (
            new scalarField(patchWave(mesh, patchSet, false).distance())
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(mesh.nCells(), 1.0));
    }
}


void Foam::inverseDistanceDiffusivity::correct()
{
    const fvMesh& mesh = mSolver().mesh();

    volScalarField y_
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );
    y_.internalField() = y();
    y_.correctBoundaryConditions();

    faceDiffusivity_ = 1.0/fvc::interpolate(y_);
}


// ************************************************************************* //
