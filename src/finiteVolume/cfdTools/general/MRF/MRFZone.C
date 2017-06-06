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

#include "MRFZone.H"
#include "cylindricalCS.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "faceSet.H"
#include "geometricOneField.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::MRFZone, 0);


const Foam::debug::tolerancesSwitch Foam::MRFZone::rotAngle_
(
    "MRFZoneRotAngle",
    10
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0: not in zone
    //  1: moving with frame
    //  2: other
    labelList faceType(mesh_.nFaces(), 0);

    // Determine faces in cell zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    if (cellZoneID_ != -1)
    {
        const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
        forAll (cellLabels, i)
        {
            zoneCell[cellLabels[i]] = true;
        }
    }

    label nZoneFaces = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (zoneCell[own[faceI]] || zoneCell[nei[faceI]])
        {
            // By default, set type "other" for faces in MRF zone
            faceType[faceI] = 1;
            nZoneFaces++;
        }
    }


    labelHashSet excludedPatches(excludedPatchLabels_);

    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.isWall() && !excludedPatches.found(patchI))
        {
            forAll (pp, i)
            {
                label faceI = pp.start() + i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 1;
                    nZoneFaces++;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll (pp, i)
            {
                label faceI = pp.start() + i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 2;
                    nZoneFaces++;
                }
            }
        }
    }

    // Now we have for faceType:
    //  0   : face not in cellZone
    //  1   : internal face or wall patch or empty patch
    //  2   : coupled patch face or excluded patch face

    // Sort into lists per patch.

    internalFaces_.setSize(mesh_.nFaces());
    label nInternal = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceType[faceI] == 1)
        {
            internalFaces_[nInternal++] = faceI;
        }
    }
    internalFaces_.setSize(nInternal);

    labelList nIncludedFaces(patches.size(), 0);
    labelList nExcludedFaces(patches.size(), 0);

    register label faceI;

    forAll (patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll (pp, patchFaceI)
        {
            faceI = pp.start() + patchFaceI;

            if (faceType[faceI] == 1)
            {
                nIncludedFaces[patchi]++;
            }
            else if (faceType[faceI] == 2)
            {
                nExcludedFaces[patchi]++;
            }
        }
    }

    includedFaces_.setSize(patches.size());
    excludedFaces_.setSize(patches.size());
    forAll (nIncludedFaces, patchI)
    {
        includedFaces_[patchI].setSize(nIncludedFaces[patchI]);
        excludedFaces_[patchI].setSize(nExcludedFaces[patchI]);
    }
    nIncludedFaces = 0;
    nExcludedFaces = 0;

    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll (pp, patchFaceI)
        {
            faceI = pp.start() + patchFaceI;

            if (faceType[faceI] == 1)
            {
                includedFaces_[patchI][nIncludedFaces[patchI]++] = patchFaceI;
            }
            else if (faceType[faceI] == 2)
            {
                excludedFaces_[patchI][nExcludedFaces[patchI]++] = patchFaceI;
            }
        }
    }

    if (debug)
    {
        faceSet internalFaces
        (
            mesh_,
            "internalFaces",
            labelHashSet(internalFaces_)
        );
        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;
        internalFaces.write();

        faceSet MRFFaces(mesh_, "includedFaces", 100);

        forAll (includedFaces_, patchI)
        {
            forAll (includedFaces_[patchI], i)
            {
                label patchFaceI = includedFaces_[patchI][i];
                MRFFaces.insert(patches[patchI].start() + patchFaceI);
            }
        }
        Pout<< "Writing " << MRFFaces.size()
            << " patch faces in MRF zone to faceSet "
            << MRFFaces.name() << endl;
        MRFFaces.write();

        faceSet excludedFaces(mesh_, "excludedFaces", 100);
        forAll (excludedFaces_, patchI)
        {
            forAll (excludedFaces_[patchI], i)
            {
                label patchFaceI = excludedFaces_[patchI][i];
                excludedFaces.insert(patches[patchI].start() + patchFaceI);
            }
        }
        Pout<< "Writing " << excludedFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << excludedFaces.name() << endl;
        excludedFaces.write();
    }
}


void Foam::MRFZone::calcMeshVelocity() const
{
    // Calculate mesh velocity from deformed mesh without executing mesh motion
    // HJ, 6/Jun/2017

    if (debug)
    {
        InfoIn("void MRFZone::calcMeshVelocity() const")
            << "Calculating mesh velocity for zone " << name_
            << endl;
    }

    if (meshVelocityPtr_)
    {
        FatalErrorIn("void MRFZone::calcMeshVelocity() const")
            << "Mesh velocity for zone " << name_
            << " already calculated"
            << abort(FatalError);
    }

    // Create the mesh velocity
    meshVelocityPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "meshVelocity" + name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity*dimArea, 0)
        );
    surfaceScalarField& meshVel = *meshVelocityPtr_;

    // Record the time of creation of mesh velocity
    meshVelTime_ = mesh_.time().value();

    // If there is no rotation, return
    if (omega_.value() < SMALL)
    {
        return;
    }

    // Create moving points mask
    scalarField movingPointsMask(mesh_.allPoints().size(), 0);

    const pointField& p = mesh_.allPoints();
    const faceList& f = mesh_.allFaces();
    const cellList& c = mesh_.cells();

    const labelList& cellAddr = mesh_.cellZones()[cellZoneID_];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    // Create a rotational coordinate system

    // Create direction by transposing axis
    vector dir;

    if (axis_.value().x() > SMALL || axis_.value().y() > SMALL)
    {
        dir = vector(axis_.value().y(), axis_.value().x(), axis_.value().z());
    }
    else if (axis_.value().z() > SMALL)
    {
        dir = vector(axis_.value().x(), axis_.value().z(), axis_.value().y());
    }
    else
    {
        FatalErrorIn("void MRFZone::calcMeshVelocity() const")
            << "Cannot transpose axis: " << axis_.value()
            << " for MRF zone " << name_
            << abort(FatalError);
    }

    cylindricalCS cs
    (
        "cs",
        origin_.value(),
        axis_.value(),
        dir,
        false            // inDegrees, for consistency
    );

    // Calculate time-step to be used in rotation from rotAngle_
    // Consider using ramped omega?  HJ, 6/Jun/2017
    const scalar deltaT =
        rotAngle_()*mathematicalConstant::pi/(180*omega_.value());

    if (debug)
    {
        InfoIn("void MRFZone::calcMeshVelocity() const")
            << "deltaT for zone " << name_ << " = " << deltaT
            << endl;
    }

    // Calculate new points
    const vectorField newP =
        cs.globalPosition
        (
            cs.localPosition(p)
          + vector(0, mag(Omega())*deltaT, 0)*movingPointsMask
        );

    // Calculate mesh velocity for all moving faces
    register label faceI, patchFaceI;

    scalarField& meshVelIn = meshVel.internalField();

    forAll (internalFaces_, i)
    {
        faceI = internalFaces_[i];
        meshVelIn[faceI] = f[faceI].sweptVol(p, newP)/deltaT;
    }

    // Included patches

    forAll (includedFaces_, patchI)
    {
        const label patchStart = mesh_.boundaryMesh()[patchI].start();

        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            meshVel.boundaryField()[patchI][patchFaceI] =
                f[patchStart + patchFaceI].sweptVol(p, newP)/deltaT;
        }
    }

    // Excluded patches
    forAll (excludedFaces_, patchI)
    {
        const label patchStart = mesh_.boundaryMesh()[patchI].start();

        forAll (excludedFaces_[patchI], i)
        {
            patchFaceI = excludedFaces_[patchI][i];

            meshVel.boundaryField()[patchI][patchFaceI] =
                f[patchStart + patchFaceI].sweptVol(p, newP)/deltaT;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone(const fvMesh& mesh, Istream& is)
:
    name_(is),
    mesh_(mesh),
    dict_(is),
    cellZoneID_(mesh_.cellZones().findZoneID(name_)),
    excludedPatchNames_
    (
        dict_.lookupOrDefault("nonRotatingPatches", wordList(0))
    ),
    origin_(dict_.lookup("origin")),
    axis_(dict_.lookup("axis")),
    omega_(dict_.lookup("omega")),
    rampTime_(dict_.lookupOrDefault<scalar>("rampTime", 0)),
    meshVelocityPtr_(NULL),
    meshVelTime_(-1)
{
    if (dict_.found("patches"))
    {
        WarningIn("MRFZone(const fvMesh&, Istream&)")
            << "Ignoring entry 'patches'\n"
            << "    By default all patches within the rotating region rotate."
            << nl << "    Optionally supply excluded patches using "
            << "'nonRotatingPatches'."
            << endl;
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (mag(axis_.value()) < SMALL)
    {
        FatalErrorIn("MRFZone(const fvMesh&, Istream&)")
            << "Axis vector has zero magnitude: " << axis_
            << ".  This is not allowed"
            << abort(FatalError);
    }

    axis_ = axis_/mag(axis_);

    excludedPatchLabels_.setSize(excludedPatchNames_.size());

    forAll (excludedPatchNames_, i)
    {
        excludedPatchLabels_[i] = patches.findPatchID(excludedPatchNames_[i]);

        if (excludedPatchLabels_[i] == -1)
        {
            FatalErrorIn
            (
                "Foam::MRFZone::MRFZone(const fvMesh&, Istream&)"
            )   << "cannot find MRF patch " << excludedPatchNames_[i]
                << exit(FatalError);
        }
    }

    bool cellZoneFound = (cellZoneID_ != -1);
    reduce(cellZoneFound, orOp<bool>());

    if (!cellZoneFound)
    {
        FatalErrorIn
        (
            "Foam::MRFZone::MRFZone(const fvMesh&, Istream&)"
        )   << "cannot find MRF cellZone " << name_
            << exit(FatalError);
    }


    Info<< "Creating MRF for cell zone " << name_ << ". rpm = "
        << 30*omega_.value()/mathematicalConstant::pi
        << endl;

    setMRFFaces();
}


const Foam::surfaceScalarField& Foam::MRFZone::meshVelocity() const
{
    // Check time
    if
    (
        rampTime_ > SMALL
     && mesh_.time().value() < rampTime_
     && mag(meshVelTime_ - mesh_.time().value()) > SMALL
    )
    {
        // MRF is ramping.  Recalculate mesh velocity
        Info<< "Clearing mesh velocity" << endl;
        deleteDemandDrivenData(meshVelocityPtr_);
    }

    // Calculate mesh velocity
    if (!meshVelocityPtr_)
    {
        calcMeshVelocity();
    }

    return *meshVelocityPtr_;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::MRFZone::~MRFZone()
{
    deleteDemandDrivenData(meshVelocityPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFZone::Omega() const
{
    if (rampTime_ < SMALL)
    {
        return omega_.value()*axis_.value();
    }
    else
    {
        // Ramping
        const scalar t = mesh_.time().value();
        const scalar ramp = sin(2*pi/(4*rampTime_)*Foam::min(rampTime_, t));
        Info<< "ramp: " << ramp << endl;
        return ramp*omega_.value()*axis_.value();
    }
}


void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector rotVel = Omega();

    forAll (cells, i)
    {
        label celli = cells[i];
        Usource[celli] -= V[celli]*(rotVel ^ U[celli]);
    }
}


void Foam::MRFZone::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector rotVel = Omega();

    forAll (cells, i)
    {
        label celli = cells[i];
        Usource[celli] -= V[celli]*rho[celli]*(rotVel ^ U[celli]);
    }
}


void Foam::MRFZone::relativeVelocity(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll (cells, i)
    {
        label celli = cells[i];
        U[celli] -= (rotVel ^ (C[celli] - origin));
    }

    // Included faces
    forAll (includedFaces_, patchi)
    {
        forAll (includedFaces_[patchi], i)
        {
            label patchFaceI = includedFaces_[patchi][i];
            U.boundaryField()[patchi][patchFaceI] = vector::zero;
        }
    }

    // Excluded faces
    forAll (excludedFaces_, patchi)
    {
        forAll (excludedFaces_[patchi], i)
        {
            label patchFaceI = excludedFaces_[patchi][i];
            U.boundaryField()[patchi][patchFaceI] -=
                (rotVel ^ (C.boundaryField()[patchi][patchFaceI] - origin));
        }
    }
}


void Foam::MRFZone::absoluteVelocity(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll (cells, i)
    {
        label celli = cells[i];
        U[celli] += (rotVel ^ (C[celli] - origin));
    }

    // Included faces
    forAll (includedFaces_, patchi)
    {
        forAll (includedFaces_[patchi], i)
        {
            label patchFaceI = includedFaces_[patchi][i];
            U.boundaryField()[patchi][patchFaceI] =
                (rotVel ^ (C.boundaryField()[patchi][patchFaceI] - origin));
        }
    }

    // Excluded faces
    forAll (excludedFaces_, patchi)
    {
        forAll (excludedFaces_[patchi], i)
        {
            label patchFaceI = excludedFaces_[patchi][i];
            U.boundaryField()[patchi][patchFaceI] +=
                (rotVel ^ (C.boundaryField()[patchi][patchFaceI] - origin));
        }
    }
}


void Foam::MRFZone::relativeFlux(surfaceScalarField& phi) const
{
    relativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::relativeFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    relativeRhoFlux(rho, phi);
}


void Foam::MRFZone::absoluteFlux(surfaceScalarField& phi) const
{
    absoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::absoluteFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    absoluteRhoFlux(rho, phi);
}


void Foam::MRFZone::meshPhi
(
    surfaceScalarField& phi
) const
{
    phi += meshVelocity();
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    register label patchFaceI;

    // Included patches
    forAll (includedFaces_, patchI)
    {
        const vectorField& patchC = mesh_.Cf().boundaryField()[patchI];

        vectorField pfld(U.boundaryField()[patchI]);

        forAll (includedFaces_[patchI], i)
        {
            patchFaceI = includedFaces_[patchI][i];

            pfld[patchFaceI] = (rotVel ^ (patchC[patchFaceI] - origin));
        }

        U.boundaryField()[patchI] == pfld;
    }
}


void Foam::MRFZone::Su
(
    const volScalarField& phi,
    const volVectorField& gradPhi,
    volScalarField& source
) const
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    forAll (cells, i)
    {
        source[cells[i]] =
            (rotVel ^ (C[cells[i]] - origin)) & gradPhi[cells[i]];
    }
}


void Foam::MRFZone::Su
(
    const volVectorField& phi,
    const volTensorField& gradPhi,
    volVectorField& source
) const
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const volVectorField& C = mesh_.C();

    const vector& origin = origin_.value();
    const vector rotVel = Omega();

    forAll (cells, i)
    {
        source[cells[i]] =
            ((rotVel ^ (C[cells[i]] - origin)) & gradPhi[cells[i]])
          - (rotVel ^ phi[cells[i]]);
    }
}


// ************************************************************************* //
