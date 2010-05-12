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

#include "displacementSBRStressFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementSBRStressFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        displacementSBRStressFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::displacementSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    fvMotionSolver(mesh),
    points0_
    (
        pointIOField
        (
            IOobject
            (
                "points",
                time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh_
    ),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement_.dimensions(),
            vector::zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(*this, lookup("diffusivity"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementSBRStressFvMotionSolver::
~displacementSBRStressFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementSBRStressFvMotionSolver::curPoints() const
{
    vpi_.interpolate(cellDisplacement_, pointDisplacement_);

    tmp<pointField> tcurPoints
    (
        points0_ + pointDisplacement_.internalField()
    );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::displacementSBRStressFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();
    pointDisplacement_.boundaryField().updateCoeffs();

    surfaceScalarField Df = diffusivityPtr_->operator()();

    volTensorField gradCd = fvc::grad(cellDisplacement_);

    Foam::solve
    (
        fvm::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               (
                   cellDisplacement_.mesh().Sf()
                 & fvc::interpolate(gradCd.T() - gradCd)
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
            )
        )

        /*
      - fvc::laplacian
        (
            2*Df,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            Df
           *(
               (
                   cellDisplacement_.mesh().Sf()
                 & fvc::interpolate(gradCd + gradCd.T())
               )

               // Solid-body rotation "lambda" term
             - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
           )
        )
        */
    );
}


void Foam::displacementSBRStressFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    fvMotionSolver::updateMesh(mpm);

    // Map points0_
    // Map points0_. Bit special since we somehow have to come up with
    // a sensible points0 position for introduced points.
    // Find out scaling between points0 and current points

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        mpm.hasMotionPoints()
      ? mpm.preMotionPoints()
      : fvMesh_.points()
    );

    // Note: boundBox does reduce
    const boundBox bb0(points0_, true);
    const vector span0(bb0.max()-bb0.min());
    const boundBox bb(points, true);
    const vector span(bb.max()-bb.min());

    vector scaleFactors(cmptDivide(span0, span));

    pointField newPoints0(mpm.pointMap().size());

    forAll(newPoints0, pointI)
    {
        label oldPointI = mpm.pointMap()[pointI];

        if (oldPointI >= 0)
        {
            label masterPointI = mpm.reversePointMap()[oldPointI];

            if (masterPointI == pointI)
            {
                newPoints0[pointI] = points0_[oldPointI];
            }
            else
            {
                // New point. Assume motion is scaling.
                newPoints0[pointI] =
                    points0_[oldPointI]
                  + cmptMultiply
                    (
                        scaleFactors,
                        points[pointI]-points[masterPointI]
                    );
            }
        }
        else
        {
            FatalErrorIn
            (
                "displacementSBRStressFvMotionSolver::updateMesh"
                "(const mapPolyMesh& mpm)"
            )   << "Cannot work out coordinates of introduced vertices."
                << " New vertex " << pointI << " at coordinate "
                << points[pointI] << exit(FatalError);
        }
    }
    points0_.transfer(newPoints0);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(NULL);
    diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
