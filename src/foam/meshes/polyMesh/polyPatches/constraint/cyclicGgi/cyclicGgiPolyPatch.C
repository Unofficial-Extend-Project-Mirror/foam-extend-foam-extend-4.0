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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "RodriguesRotation.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "standAlonePatch.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGgiPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicGgiPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicGgiPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cyclicGgiPolyPatch::checkDefinition() const
{
    // A little bit of sanity check The rotation angle/axis is
    // specified in both the master and slave patch of the
    // cyclicGgi.  This is a pain, but the other alternatives
    // would be:
    //
    //   - Specify in only of the two patches boundary
    //   definition : - which one to chose?  - which default
    //   value to chose for the non-initialized value - Use a
    //   specific dictionary for this... Nope, too cumbersome.
    //
    // So, we impose that the boundary definition of both
    // patches must specify the same information If not, well,
    // we stop the simulation and ask for a fix.

    if (!active())
    {
        // No need to check anything, the shadow is not initialized properly.
        // This will happen with blockMesh when defining cyclicGGI patches.
        // Return quietly
        return;
    }

    if
    (
        (mag(rotationAngle()) - mag(cyclicShadow().rotationAngle())) > SMALL
     || cmptSum(rotationAxis() - cyclicShadow().rotationAxis()) > SMALL
    )
    {
        FatalErrorIn("void cyclicGgiPolyPatch::check() const")
            << "    Rotation angle for patch name           : "
            << name() << " is: " << rotationAngle()
            << " axis: " << rotationAxis() << nl
            << "    Rotation angle for shadow patch name: "
            << shadowName() << " is: "
            << cyclicShadow().rotationAngle() << " axis: "
            << cyclicShadow().rotationAxis() << nl
            << "    Both values need to be opposite in "
            << "the boundary file. "
            << abort(FatalError);
    }

    if
    (
        (mag(separationOffset() + cyclicShadow().separationOffset())) > SMALL
    )
    {
        FatalErrorIn("void cyclicGgiPolyPatch::check() const")
            << "Separation offset for patch name           : "
            << name() << " is: " << separationOffset()
            << "    Separation offset for shadow patch name: "
            << shadowName() << " is: "
            << cyclicShadow().separationOffset() << " axis: "
            << "    Both values need to be opposite in "
            << "the boundary file. "
            << abort(FatalError);
    }

    if (debug > 1 && master())
    {
        Info<< "Writing transformed slave patch as VTK." << nl
            << "Master: " << name()
            << " Slave: " << shadowName()
            << " Angle (master to slave): " << rotationAngle() << " deg"
            << " Axis: " << rotationAxis()
            << " Separation: " << separationOffset() << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            pointField transformedPoints = cyclicShadow().localPoints();

            tensor rot = RodriguesRotation(rotationAxis_,  -rotationAngle_);

            transform(transformedPoints, rot, transformedPoints);

            // Add separation offset to transformed points.  HJ, 24/Nov/2009
            transformedPoints += cyclicShadow().separationOffset();

            standAlonePatch::writeVTK
            (
                fvPath/fileName("cyclicGgi" + name() + cyclicShadow().name()),
                cyclicShadow().localFaces(),
                transformedPoints
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicGgiPolyPatch::cyclicGgiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    ggiPolyPatch(name, size, start, index, bm),
    separationOffset_(vector::zero),
    rotationAxis_(vector(0.0, 0.0, 1.0)),
    rotationAngle_(0.0)
{}


// Construct from components
Foam::cyclicGgiPolyPatch::cyclicGgiPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowName,
    const word& zoneName,
    const bool bridgeOverlap,
    const vector& separationOffset,
    const vector& rotationAxis,
    const scalar rotationAngle
)
:
    ggiPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        shadowName,
        zoneName,
        bridgeOverlap
    ),
    separationOffset_(separationOffset),
    rotationAxis_(rotationAxis),
    rotationAngle_(rotationAngle)
{}


// Construct from dictionary
Foam::cyclicGgiPolyPatch::cyclicGgiPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    ggiPolyPatch(name, dict, index, bm),
    separationOffset_(dict.lookup("separationOffset")),
    rotationAxis_(dict.lookup("rotationAxis")),
    rotationAngle_(readScalar(dict.lookup("rotationAngle")))
{}


// Construct as copy, resetting the boundary mesh
Foam::cyclicGgiPolyPatch::cyclicGgiPolyPatch
(
    const cyclicGgiPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    ggiPolyPatch(pp, bm),
    separationOffset_(pp.separationOffset_),
    rotationAxis_(pp.rotationAxis_),
    rotationAngle_(pp.rotationAngle_)
{}


// Construct as copy, resetting the face list and boundary mesh data
Foam::cyclicGgiPolyPatch::cyclicGgiPolyPatch
(
    const cyclicGgiPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    ggiPolyPatch(pp, bm, index, newSize, newStart),
    separationOffset_(pp.separationOffset_),
    rotationAxis_(pp.rotationAxis_),
    rotationAngle_(pp.rotationAngle_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicGgiPolyPatch::~cyclicGgiPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicGgiPolyPatch& Foam::cyclicGgiPolyPatch::cyclicShadow() const
{
    return refCast<const cyclicGgiPolyPatch>(boundaryMesh()[shadowIndex()]);
}


void Foam::cyclicGgiPolyPatch::calcTransforms() const
{
    if (active() && debug)
    {
        // Check definition of the cyclic pair
        checkDefinition();
    }

    // For computing the rotation tensors forwardT and reverseT, we
    // can see from Robert Magnan's post on the Forum dated
    // 27/May/2008 that we cannot use the actual implementation of
    // calcTransformTensors and rotationTensor.
    //
    // It is also not possible to use Robert solution because for
    // non-conformal cyclic meshes, we cannot usualy find a pair of
    // matching faces for computing the tensors. So we are using
    // user-supplied values instead.
    //
    // Since these tensors are defined as private in the
    // coupledPolyPatch class definition, we will override these
    // values here and compute the tensors ourselves using the
    // Rodrigues Rotation formula.

    // We compute the rotationTensor from rotationAxis_ and
    // rotationAngle_ We compute the separation vector from
    // separationOffset_

    // All transforms are constant: size = 1.  HJ, 18/Feb/2009

    if (mag(rotationAngle_) > SMALL)
    {
        // Rotation tensor computed from rotationAxis_ and rotationAngle_
        // Note: cyclics already have opposing signs for the rotation
        //       so there is no need for a special practice.  HJ, 30/Jun/2009
        forwardT_ = tensorField
        (
            1,
            RodriguesRotation(rotationAxis_,  -rotationAngle_)
        );

        reverseT_ = tensorField
        (
            1,
            RodriguesRotation(rotationAxis_, rotationAngle_)
        );
    }
    else
    {
        forwardT_.setSize(0);
        reverseT_.setSize(0);
    }

    // Handling the separation offset separatly
    if (mag(separationOffset_) > SMALL)
    {
        separation_ = vectorField(1, separationOffset_);
    }
    else
    {
        separation_.setSize(0);
    }

    if (debug > 1 && master())
    {
        if (patchToPatch().uncoveredMasterFaces().size() > 0)
        {
            // Write uncovered master faces
            Info<< "Writing uncovered master faces for patch "
                << name() << " as VTK." << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            indirectPrimitivePatch::writeVTK
            (
                fvPath/fileName("uncoveredCyclicGgiFaces" + name()),
                IndirectList<face>
                (
                    localFaces(),
                    patchToPatch().uncoveredMasterFaces()
                ),
                localPoints()
            );
        }

        if (patchToPatch().uncoveredSlaveFaces().size() > 0)
        {
            // Write uncovered master faces
            Info<< "Writing uncovered shadow faces for patch "
                << shadowName() << " as VTK." << endl;

            const polyMesh& mesh = boundaryMesh().mesh();

            fileName fvPath(mesh.time().path()/"VTK");
            mkDir(fvPath);

            indirectPrimitivePatch::writeVTK
            (
                fvPath/fileName("uncoveredCyclicGgiFaces" + shadowName()),
                IndirectList<face>
                (
                    shadow().localFaces(),
                    patchToPatch().uncoveredSlaveFaces()
                ),
                shadow().localPoints()
            );
        }

        // Check for bridge overlap
        if (!bridgeOverlap())
        {
            if
            (
                patchToPatch().uncoveredMasterFaces().size() > 0
             || patchToPatch().uncoveredSlaveFaces().size() > 0
            )
            {
                FatalErrorIn("label cyclicGgiPolyPatch::shadowIndex() const")
                    << "cyclic ggi patch " << name() << " with shadow "
                    << shadowName() << " has "
                    << patchToPatch().uncoveredMasterFaces().size()
                    << " uncovered master faces and "
                    << patchToPatch().uncoveredSlaveFaces().size()
                    << " uncovered slave faces.  Bridging is switched off. "
                    << abort(FatalError);
            }
        }
    }
}


// Write
void Foam::cyclicGgiPolyPatch::write(Ostream& os) const
{
    ggiPolyPatch::write(os);
    os.writeKeyword("rotationAxis") << rotationAxis_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAngle") << rotationAngle_
        << token::END_STATEMENT << nl;
    os.writeKeyword("separationOffset") << separationOffset_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
