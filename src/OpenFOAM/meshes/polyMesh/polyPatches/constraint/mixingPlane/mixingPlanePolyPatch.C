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

Author
    Martin Beaudoin, Hydro-Quebec, 2009. All rights reserved.

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlanePolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "Time.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixingPlanePolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcPatchToPatch() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        // Create dummy interpolation profile
        pointField iProfile;

        if
        (
            assemblyType_
         == MixingPlaneInterpolationName::USER_DEFINED
        )
        {
            Info<< "Reading interpolation profile from file: "
                << userProfileFile_ << endl;

            iProfile = vectorIOField
            (
                IOobject
                (
                    userProfileFile_,
                    boundaryMesh().mesh().time().constant(),
                    boundaryMesh().mesh().time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            if (iProfile.empty())
            {
                FatalErrorIn
                (
                    "void mixingPlanePolyPatch::calcPatchToPatch() const"
                )   << "Empty user-defined mixing plane profile for patch "
                    << name() << " read from file " << userProfileFile_
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            Info<< "Creating mixingPlaneInterpolation for patch "
                << name() << " with shadow " << shadowName() << nl
                << "assemblyType = "
                << MixingPlaneInterpolationName::assemblyNames_[assemblyType_]
                << " " << assemblyType_
                << " orientationType = "
                << MixingPlaneInterpolationName::orientationNames_
                   [orientationType_]
                << " " << orientationType_
                << endl;

        }

        patchToPatchPtr_ =
            new mixingPlaneInterpolation
            (
                *this,
                shadow(),
                csPtr_(),
                assemblyType_,
                orientationType_,
                iProfile
            );
    }
    else
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcPatchToPatch() const")
            << "Attempting to create MixingPlaneInterpolation on a shadow"
            << abort(FatalError);
    }

    if (debug > 1 && master())
    {
        Info<< "Writing transformed mixing plane patches as VTK." << nl
            << "Master: " << name()
            << " Slave: " << shadowName()
            << endl;

        const polyMesh& mesh = boundaryMesh().mesh();

        fileName fvPath(mesh.time().path()/"VTK");
        mkDir(fvPath);

        patchToPatchPtr_->mixingPlanePatch().writeVTK
        (
            fvPath/fileName
            (
                "mixingPlaneRibbon_" + name() + "_" + shadow().name()
            )
        );

        patchToPatchPtr_->transformedMasterPatch().writeVTK
        (
            fvPath/fileName
            (
                "mixingPlaneMaster_" + name() + "_" + shadow().name()
            )
        );

        patchToPatchPtr_->transformedShadowPatch().writeVTK
        (
            fvPath/fileName
            (
                "mixingPlaneShadow_" + name() + "_" + shadow().name()
            )
        );
    }
}


void Foam::mixingPlanePolyPatch::calcReconFaceCellCentres() const
{
    // Create neighbouring face centres using interpolation
    if (master())
    {
        reconFaceCellCentresPtr_ = new vectorField
        (
            interpolate
            (
                shadow().faceCellCentres()
              - shadow().faceCentres()
            )
          + faceCentres()
        );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlanePolyPatch::calcReconFaceCellCentres() const"
        )   << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::mixingPlanePolyPatch::clearOut()
{
    deleteDemandDrivenData(patchToPatchPtr_);
    deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_(fileName::null),
    csPtr_
    (
        new coordinateSystem
        (
            "mixingCS",
            vector::zero,
            vector(0, 0, 1),
            vector(1, 0, 0)
        )
    ),
    assemblyType_(mixingPlaneInterpolation::USER_DEFINED),
    orientationType_(mixingPlaneInterpolation::UNKNOWN),
    userProfileFile_(fileName::null),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& shadowName,
    const coordinateSystem& cs,
    const mixingPlaneInterpolation::assembly assemblyType,
    const mixingPlaneInterpolation::orientation orientationType,
    const fileName& userProfileFile
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadowName_(shadowName),
    csPtr_(cs.clone()),
    assemblyType_(assemblyType),
    orientationType_(orientationType),
    userProfileFile_(fileName::null),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    shadowName_(dict.lookup("shadowPatch")),
    csPtr_
    (
        new coordinateSystem
        (
            "mixingCS",
            vector::zero,
            vector(0, 0, 1),
            vector(1, 0, 0)
        )
    ),
    assemblyType_(mixingPlaneInterpolation::USER_DEFINED),
    orientationType_(mixingPlaneInterpolation::UNKNOWN),
    userProfileFile_(fileName::null),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{
    // When construting from dictionary, only master side information will be
    // read and used.  This requires special check, because polyBoundaryMesh
    // is being filled at this point.  See fix in findPatchID.
    // HJ, and MB, 28/Jan/2011

    // Check if shadow exists.  If so, we are on the slave side
    polyPatchID shadow(shadowName_, boundaryMesh());

    if (!shadow.active())
    {
        // Master side, read additional data

        csPtr_ =
            coordinateSystem::New
            (
                "mixingCS",
                dict.subDict("coordinateSystem")
            );

        assemblyType_ =
            MixingPlaneInterpolationName::assemblyNames_.read
            (
                dict.lookup("assembly")
            );

        orientationType_ =
            MixingPlaneInterpolationName::orientationNames_.read
            (
                dict.lookup("orientation")
            );

        if (assemblyType_ == MixingPlaneInterpolationName::USER_DEFINED)
        {
            if (dict.found("userProfileFile"))
            {
                dict.lookup("userProfileFile") >> userProfileFile_;
            }
            else
            {
                FatalIOErrorIn
                (
                    "mixingPlanePolyPatch::mixingPlanePolyPatch\n"
                    "(\n"
                    "    const word& name,\n"
                    "    const dictionary& dict,\n"
                    "    const label index,\n"
                    "    const polyBoundaryMesh& bm\n"
                    ")",
                    dict
                )   << "Patch: " << name << " : Missing profile entry for "
                    << "userDefined profile"
                    << abort(FatalIOError);
            }
        }
    }
}


Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadowName_(pp.shadowName_),
    csPtr_(pp.csPtr_->clone()),
    assemblyType_(pp.assemblyType_),
    orientationType_(pp.orientationType_),
    userProfileFile_(pp.userProfileFile_),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


//- Construct as copy, resetting the face list and boundary mesh data
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    shadowName_(pp.shadowName_),
    csPtr_(pp.csPtr_->clone()),
    assemblyType_(pp.assemblyType_),
    orientationType_(pp.orientationType_),
    userProfileFile_(pp.userProfileFile_),
    shadowIndex_(-1),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingPlanePolyPatch::~mixingPlanePolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::mixingPlanePolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1 && shadowName_ != Foam::word::null)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadowName_, boundaryMesh());

        if (!shadow.active())
        {
            WarningIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowName_
                << " not found.  Please check your MixingPlane definition.  "
                << "This may be fine at mesh generation stage."
                << endl;

            // Return a large label to indicate "undefined" or slave side
            return 99999;
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a mixingPlane
        if (!isA<mixingPlanePolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "Shadow of mixingPlane patch " << name()
                << " named " << shadowName() << " is not a mixingPlane.  "
                << "Type: " << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        // Check for mixingPlane onto self
        if (index() == shadowIndex_)
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "mixingPlane patch " << name()
                << " created as its own shadow"
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}


const Foam::mixingPlanePolyPatch& Foam::mixingPlanePolyPatch::shadow() const
{
    return refCast<const mixingPlanePolyPatch>(boundaryMesh()[shadowIndex()]);
}


const Foam::coordinateSystem& Foam::mixingPlanePolyPatch::cs() const
{
    if (master())
    {
        return csPtr_();
    }
    else
    {
        return shadow().cs();
    }
}


const Foam::mixingPlaneInterpolation&
Foam::mixingPlanePolyPatch::patchToPatch() const
{
    if (master())
    {
        if (!patchToPatchPtr_)
        {
            Info<< "Initializing the mixingPlane interpolator between "
                << "master/shadow patches: "
                << name() << "/" << shadowName()
                << endl;


            calcPatchToPatch();
        }

        return *patchToPatchPtr_;
    }
    else
    {
        return shadow().patchToPatch();
    }
}


const Foam::vectorField&
Foam::mixingPlanePolyPatch::reconFaceCellCentres() const
{
    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}


void Foam::mixingPlanePolyPatch::initAddressing()
{
    polyPatch::initAddressing();
}


void Foam::mixingPlanePolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::mixingPlanePolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::mixingPlanePolyPatch::calcGeometry()
{
    // Reconstruct the cell face centres
    if (patchToPatchPtr_ && master())
    {
        // Compute the neighbour face cell center
        reconFaceCellCentres();

        // Next, identify which cells are located at these locations

        // Next, compute the weighting factors in order to properly interpolate
        // the field values at those locations. We will be using an inverse
        // distance interpolation scheme.
    }

    calcTransforms();
    polyPatch::calcGeometry();
}


void Foam::mixingPlanePolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::mixingPlanePolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    clearOut();
}


void Foam::mixingPlanePolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::mixingPlanePolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    clearOut();
}


void Foam::mixingPlanePolyPatch::calcTransforms()
{
    forwardT_.setSize(0);
    reverseT_.setSize(0);
    separation_.setSize(0);
}


void Foam::mixingPlanePolyPatch::initOrder(const primitivePatch& pp) const
{}


bool Foam::mixingPlanePolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size(), -1);
    rotation.setSize(pp.size(), 0);

    // Nothing changes
    return false;
}


void Foam::mixingPlanePolyPatch::syncOrder() const
{}


void Foam::mixingPlanePolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("shadowPatch") << shadowName_
        << token::END_STATEMENT << nl;

    // Note: only master writes the data
    if (master() || shadowIndex_ == -1)
    {
        // Write coordinate system dictionary.  Check by hand.  HJ, 26/Jan/2011
        os.writeKeyword("coordinateSystem");
        csPtr_().writeDict(os, true);

        os.writeKeyword("assembly")
            << MixingPlaneInterpolationName::assemblyNames_[assemblyType_]
                << token::END_STATEMENT << nl;

        os.writeKeyword("orientation")
            << MixingPlaneInterpolationName::orientationNames_
                   [orientationType_]
            << token::END_STATEMENT << nl;

        if (userProfileFile_ != fileName::null)
        {
            os.writeKeyword("userProfileFile") << userProfileFile_
                << token::END_STATEMENT << nl;
        }
    }
}


// ************************************************************************* //
