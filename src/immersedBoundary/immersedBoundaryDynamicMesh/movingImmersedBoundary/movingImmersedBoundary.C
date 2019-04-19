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

\*---------------------------------------------------------------------------*/

#include "movingImmersedBoundary.H"
#include "immersedBoundaryPolyPatch.H"
#include "mixedIbFvPatchFields.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingImmersedBoundary::movingImmersedBoundary
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    sbmfPtr_(solidBodyMotionFunction::New(dict, mesh.time())),
    refIbSurface_
    (
        IOobject
        (
            name  + ".ftr",
            mesh.time().constant(),      // instance
            "triSurface",                // local
            mesh,                        // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingImmersedBoundary::~movingImmersedBoundary()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingImmersedBoundary::movePoints() const
{
    // Get ibMesh from patch
    const label patchID = mesh().boundaryMesh().findPatchID(name());

    if (patchID < 0)
    {
        FatalErrorIn
        (
            "void movingImmersedBoundary::movePoints() const"
        )   << "Patch " << name() << " not found.  Available patch names: "
            << mesh().boundaryMesh().names()
            << abort(FatalError);
    }

    const immersedBoundaryPolyPatch& cibPatch =
        refCast<const immersedBoundaryPolyPatch>
        (
            mesh().boundaryMesh()[patchID]
        );

    // Get non-const reference to patch
    immersedBoundaryPolyPatch& ibPatch =
        const_cast<immersedBoundaryPolyPatch&>(cibPatch);

    // Move points
    ibPatch.moveTriSurfacePoints
    (
        transform(sbmfPtr_->transformation(), refIbSurface_.points())
    );
}


// ************************************************************************* //
