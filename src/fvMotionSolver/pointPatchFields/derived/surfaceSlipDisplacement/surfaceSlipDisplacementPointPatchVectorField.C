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

#include "surfaceSlipDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "fvMesh.H"
#include "displacementLaplacianFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
NamedEnum<surfaceSlipDisplacementPointPatchVectorField::followMode, 3>::
names[] =
{
    "nearest",
    "pointNormal",
    "fixedNormal"
};

const NamedEnum<surfaceSlipDisplacementPointPatchVectorField::followMode, 3>
    surfaceSlipDisplacementPointPatchVectorField::followModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(p, iF),
    projectMode_(NEAREST),
    projectDir_(vector::zero),
    wedgePlane_(-1)
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchVectorField(p, iF),
    surfacesDict_(dict.subDict("geometry")),
    projectMode_(followModeNames_.read(dict.lookup("followMode"))),
    projectDir_(dict.lookup("projectDirection")),
    wedgePlane_(readLabel(dict.lookup("wedgePlane"))),
    frozenPointsZone_(dict.lookup("frozenPointsZone"))
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper&
)
:
    pointPatchVectorField(p, iF),
    surfacesDict_(ppf.surfacesDict()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf
)
:
    pointPatchVectorField(ppf),
    surfacesDict_(ppf.surfacesDict()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(ppf, iF),
    surfacesDict_(ppf.surfacesDict()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const searchableSurfaces& surfaceSlipDisplacementPointPatchVectorField::
surfaces() const
{
    if (!surfacesPtr_.valid())
    {
        surfacesPtr_.reset
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",                              // dummy name
                    db().time().constant(),             // directory
                    "triSurface",                       // instance
                    db().time(),                        // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfacesDict_
            )
        );
    }
    return surfacesPtr_();
}


void surfaceSlipDisplacementPointPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    const polyMesh& mesh = patch().boundaryMesh().mesh()();

    //const scalar deltaT = mesh.time().deltaT().value();

    // Construct large enough vector in direction of projectDir so
    // we're guaranteed to hit something.

    const scalar projectLen = mag(mesh.bounds().max()-mesh.bounds().min());

    // For case of fixed projection vector:
    vector projectVec;
    if (projectMode_ == FIXEDNORMAL)
    {
        vector n = projectDir_/mag(projectDir_);
        projectVec = projectLen*n;
    }

    //- Per point projection vector:

    const pointField& localPoints = patch().localPoints();
    const labelList& meshPoints = patch().meshPoints();

    vectorField displacement(this->patchInternalField());


    // Get fixed points (bit of a hack)
    const pointZone* zonePtr = NULL;

    if (frozenPointsZone_.size() > 0)
    {
        const pointZoneMesh& pZones = mesh.pointZones();

        zonePtr = &pZones[pZones.findZoneID(frozenPointsZone_)];

        Pout<< "surfaceSlipDisplacementPointPatchVectorField : Fixing all "
            << zonePtr->size() << " points in pointZone " << zonePtr->name()
            << endl;
    }

    // Get the starting locations from the motionSolver
    const displacementLaplacianFvMotionSolver& motionSolver =
        mesh.lookupObject<displacementLaplacianFvMotionSolver>
        (
            "dynamicMeshDict"
        );
    const pointField& points0 = motionSolver.points0();


//XXXXXX


    pointField start(meshPoints.size());
    forAll(start, i)
    {
        start[i] = points0[meshPoints[i]] + displacement[i];
    }

    if (projectMode_ == NEAREST)
    {
        List<pointIndexHit> nearest;
        labelList hitSurfaces;
        surfaces().findNearest
        (
            start,
            scalarField(start.size(), sqr(projectLen)),
            hitSurfaces,
            nearest
        );

        forAll(nearest, i)
        {
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                displacement[i] =
                    nearest[i].hitPoint()
                  - points0[meshPoints[i]];
            }
            else
            {
                Pout<< "    point:" << meshPoints[i]
                    << " coord:" << localPoints[i]
                    << "  did not find any surface within " << projectLen
                    << endl;
            }
        }
    }
    else
    {
        // Do tests on all points. Combine later on.

        // 1. Check if already on surface
        List<pointIndexHit> nearest;
        {
            labelList nearestSurface;
            surfaces().findNearest
            (
                start,
                scalarField(start.size(), sqr(SMALL)),
                nearestSurface,
                nearest
            );
        }

        // 2. intersection. (combined later on with information from nearest
        // above)
        vectorField projectVecs(start.size(), projectVec);

        if (projectMode_ == POINTNORMAL)
        {
            projectVecs = projectLen*patch().pointNormals();
        }

        // Knock out any wedge component
        scalarField offset(start.size(), 0.0);
        if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
        {
            forAll(offset, i)
            {
                offset[i] = start[i][wedgePlane_];
                start[i][wedgePlane_] = 0;
                projectVecs[i][wedgePlane_] = 0;
            }
        }

        List<pointIndexHit> rightHit;
        {
            labelList rightSurf;
            surfaces().findAnyIntersection
            (
                start,
                start+projectVecs,
                rightSurf,
                rightHit
            );
        }
        
        List<pointIndexHit> leftHit;
        {
            labelList leftSurf;
            surfaces().findAnyIntersection
            (
                start,
                start-projectVecs,
                leftSurf,
                leftHit
            );
        }

        // 3. Choose either -fixed, nearest, right, left.
        forAll(displacement, i)
        {
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                // Found nearest.
                displacement[i] =
                    nearest[i].hitPoint()
                  - points0[meshPoints[i]];
            }
            else
            {
                pointIndexHit interPt;

                if (rightHit[i].hit())
                {
                    if (leftHit[i].hit())
                    {
                        if
                        (
                            magSqr(rightHit[i].hitPoint()-start[i])
                          < magSqr(leftHit[i].hitPoint()-start[i])
                        )
                        {
                            interPt = rightHit[i];
                        }
                        else
                        {
                            interPt = leftHit[i];
                        }
                    }
                    else
                    {
                        interPt = rightHit[i];
                    }
                }
                else
                {
                    if (leftHit[i].hit())
                    {
                        interPt = leftHit[i];
                    }
                }


                if (interPt.hit())
                {
                    if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
                    {
                        interPt.rawPoint()[wedgePlane_] += offset[i];
                    }
                    displacement[i] = interPt.rawPoint()-points0[meshPoints[i]];
                }
                else
                {
                    Pout<< "    point:" << meshPoints[i]
                        << " coord:" << localPoints[i]
                        << "  did not find any intersection between ray from "
                        << start[i]-projectVecs[i]
                        << " to " << start[i]+projectVecs[i]
                        << endl;
                }
            }
        }
    }

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->internalField());

    //setInInternalField(iF, motionU);
    setInInternalField(iF, displacement);

    pointPatchVectorField::evaluate(commsType);
}


void surfaceSlipDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchVectorField::write(os);
    os.writeKeyword("geometry") << surfacesDict_
        << token::END_STATEMENT << nl;
    os.writeKeyword("followMode") << followModeNames_[projectMode_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("projectDirection") << projectDir_
        << token::END_STATEMENT << nl;
    os.writeKeyword("wedgePlane") << wedgePlane_
        << token::END_STATEMENT << nl;
    os.writeKeyword("frozenPointsZone") << frozenPointsZone_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    surfaceSlipDisplacementPointPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
