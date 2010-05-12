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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "faceSet.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces
(
    labelList& faceType,
    const labelList& excludedPatchIDs
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Knock out coupled patches
    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const polyPatch& pp = patches[patchi];

            forAll(pp, j)
            {
                label faceI = pp.start()+j;

                if (faceType[faceI] == 1)
                {
                    faceType[faceI] = 2;
                }
            }
        }
    }

    // All explicitly provided exclusions
    forAll(excludedPatchIDs, i)
    {
        const polyPatch& pp = patches[excludedPatchIDs[i]];

        forAll(pp, j)
        {
            label faceI = pp.start()+j;

            if (faceType[faceI] == 1)
            {
                faceType[faceI] = 2;
            }
        }
    }

    // Collect into lists per patch.
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

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label faceI = pp.start()+patchFacei;

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
    forAll(nIncludedFaces, patchi)
    {
        includedFaces_[patchi].setSize(nIncludedFaces[patchi]);
        excludedFaces_[patchi].setSize(nExcludedFaces[patchi]);
    }
    nIncludedFaces = 0;
    nExcludedFaces = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label faceI = pp.start()+patchFacei;

            if (faceType[faceI] == 1)
            {
                includedFaces_[patchi][nIncludedFaces[patchi]++] = patchFacei;
            }
            else if (faceType[faceI] == 2)
            {
                excludedFaces_[patchi][nExcludedFaces[patchi]++] = patchFacei;
            }
        }
    }

    //if (debug)
    //{
    //    faceSet internalFaces(mesh_, "internalFaces", internalFaces_);
    //    Pout<< "Writing internal faces in MRF zone to faceSet "
    //        << internalFaces.name() << endl;
    //    internalFaces.write();
    //}
    //{
    //    faceSet MRFFaces(mesh_, "includedFaces", 100);
    //    forAll(includedFaces_, patchi)
    //    {
    //        forAll(includedFaces_[patchi], i)
    //        {
    //            label patchFacei = includedFaces_[patchi][i];
    //            MRFFaces.insert(patches[patchi].start()+patchFacei);
    //        }
    //    }
    //    Pout<< "Writing patch faces in MRF zone to faceSet "
    //        << MRFFaces.name() << endl;
    //    MRFFaces.write();
    //}
    //{
    //    faceSet excludedFaces(mesh_, "excludedFaces", 100);
    //    forAll(excludedFaces_, patchi)
    //    {
    //        forAll(excludedFaces_[patchi], i)
    //        {
    //            label patchFacei = excludedFaces_[patchi][i];
    //            excludedFaces.insert(patches[patchi].start()+patchFacei);
    //        }
    //    }
    //    Pout<< "Writing faces in MRF zone with special handling to faceSet "
    //        << excludedFaces.name() << endl;
    //    excludedFaces.write();
    //}
}


void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0:not in zone
    //  1:moving with frame
    //  2:other
    labelList faceType(mesh_.nFaces(), 0);

    bool faceZoneFound = (faceZoneID_ != -1);
    reduce(faceZoneFound, orOp<bool>());

    if (faceZoneFound)
    {
        // Explicitly provided faces.
        if (faceZoneID_ != -1)
        {
            const labelList& zoneFaces = mesh_.faceZones()[faceZoneID_];

            forAll(zoneFaces, i)
            {
                faceType[zoneFaces[i]] = 1;
            }

            if (allPatchesMove_)
            {
                // Explicitly provided patches that do not move
                setMRFFaces(faceType, patchLabels_);
            }
            else
            {
                setMRFFaces(faceType, labelList(0));
            }
        }
    }
    else
    {
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
            forAll(cellLabels, i)
            {
                zoneCell[cellLabels[i]] = true;
            }
        }


        label nZoneFaces = 0;

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (zoneCell[own[faceI]] || zoneCell[nei[faceI]])
            {
                faceType[faceI] = 1;
                nZoneFaces++;
            }
        }
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (!isA<emptyPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    if (zoneCell[own[faceI]])
                    {
                        faceType[faceI] = 1;
                        nZoneFaces++;
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, faceType, maxEqOp<label>(), false);


        Info<< nl
            << "MRFZone " << name_ << " : did not find a faceZone; using "
            << returnReduce(nZoneFaces, sumOp<label>())
            << " faces from the cellZone instead." << endl;


        if (allPatchesMove_)
        {
            // Explicitly provided excluded patches
            setMRFFaces(faceType, patchLabels_);
        }
        else
        {
            setMRFFaces(faceType, labelList(0));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone(const fvMesh& mesh, Istream& is)
:
    mesh_(mesh),
    name_(is),
    dict_(is),
    cellZoneID_(mesh_.cellZones().findZoneID(name_)),
    faceZoneID_(mesh_.faceZones().findZoneID(name_)),
    allPatchesMove_(dict_.found("nonRotatingPatches")),
    patchNames_
    (
        allPatchesMove_
      ? dict_.lookup("nonRotatingPatches")
      : dict_.lookup("patches")
    ),
    origin_(dict_.lookup("origin")),
    axis_(dict_.lookup("axis")),
    omega_(dict_.lookup("omega")),
    Omega_("Omega", omega_*axis_)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    axis_ = axis_/mag(axis_);
    Omega_ = omega_*axis_;

    patchLabels_.setSize(patchNames_.size());

    forAll(patchNames_, i)
    {
        patchLabels_[i] = patches.findPatchID(patchNames_[i]);

        if (patchLabels_[i] == -1)
        {
            FatalErrorIn
            (
                "Foam::MRFZone::MRFZone(const fvMesh&, Istream&)"
            )   << "cannot find MRF patch " << patchNames_[i]
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

    setMRFFaces();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    if (UEqn.dimensions() != UEqn.psi().dimensions()*dimVolume/dimTime)
    {
        FatalErrorIn
        (
            "void MRFZone::addCoriolis(fvVectorMatrix& UEqn) const"
        )   << "Equation for " << UEqn.psi().name()
            << " is not defined in terms of velocity. UEqn: "
            << UEqn.dimensions() << " expected "
            << UEqn.psi().dimensions()*dimVolume/dimTime
            << abort(FatalError);
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi().internalField();
    const vector& Omega = Omega_.value();

    forAll(cells, i)
    {
        Usource[cells[i]] -= V[cells[i]]*(Omega ^ U[cells[i]]);
    }
}


void Foam::MRFZone::relativeFlux(surfaceScalarField& phi) const
{
    if (phi.dimensions() != dimVelocity*dimArea)
    {
        FatalErrorIn
        (
            "void MRFZone::relativeFlux(surfaceScalarField& phi) const"
        )   << phi.name() << " is not a volumetric flux field"
            << abort(FatalError);
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector& origin = origin_.value();
    const vector& Omega = Omega_.value();

    // Access to flux field.  Optimisation.  HJ, 12/Dec/2009
    scalarField& phiIn = phi.internalField();

    surfaceScalarField::GeometricBoundaryField& phiBoundary =
        phi.boundaryField();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        phiIn[facei] -= (Omega ^ (Cf[facei] - origin)) & Sf[facei];
    }

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            phiBoundary[patchi][patchFacei] = 0.0;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            phiBoundary[patchi][patchFacei] -=
                (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin))
              & Sf.boundaryField()[patchi][patchFacei];
        }
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

    if
    (
        UEqn.dimensions()
     != rho.dimensions()*UEqn.psi().dimensions()*dimVolume/dimTime
    )
    {
        FatalErrorIn
        (
            "void Foam::MRFZone::addCoriolis\n"
            "(\n"
            "    const volScalarField& rho,\n"
            "    fvVectorMatrix& UEqn\n"
            ") const"
        )   << "Equation for " << UEqn.psi().name()
            << " is not defined in terms of momentum.  UEqn: "
            << UEqn.dimensions() << " expected "
            << rho.dimensions()*UEqn.psi().dimensions()*dimVolume/dimTime
            << abort(FatalError);
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi().internalField();
    const scalarField& rhoIn = rho.internalField();
    const vector& Omega = Omega_.value();

    forAll(cells, i)
    {
        Usource[cells[i]] -= rhoIn[cells[i]]*V[cells[i]]*(Omega ^ U[cells[i]]);
    }
}


void Foam::MRFZone::relativeFlux
(
    const surfaceScalarField& rhof,
    surfaceScalarField& phi
) const
{
    if (phi.dimensions() != dimVelocity*dimArea*rhof.dimensions())
    {
        FatalErrorIn
        (
            "void MRFZone::relativeFlux(surfaceScalarField& phi) const"
        )   << "Dimensions of flux field" << phi.name()
            << " are not consistent with a flux field"
            << abort(FatalError);
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector& origin = origin_.value();
    const vector& Omega = Omega_.value();

    // Prefactor field
    const scalarField& rhofIn = rhof.internalField();

    const surfaceScalarField::GeometricBoundaryField& rhofBoundary =
        rhof.boundaryField();

    // Access to flux field.  Optimisation.  HJ, 12/Dec/2009
    scalarField& phiIn = phi.internalField();

    surfaceScalarField::GeometricBoundaryField& phiBoundary =
        phi.boundaryField();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        phiIn[facei] -=
            rhofIn[facei]*(Omega ^ (Cf[facei] - origin)) & Sf[facei];
    }

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            phiBoundary[patchi][patchFacei] = 0.0;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            phiBoundary[patchi][patchFacei] -=
                rhofBoundary[patchi][patchFacei]*
                (
                    (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin))
                  & Sf.boundaryField()[patchi][patchFacei]
                );
        }
    }
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    const vector& origin = origin_.value();
    const vector& Omega = Omega_.value();

    // Included patches
    forAll(includedFaces_, patchi)
    {
        const vectorField& patchC = mesh_.Cf().boundaryField()[patchi];

        vectorField pfld(U.boundaryField()[patchi]);

        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            pfld[patchFacei] = (Omega ^ (patchC[patchFacei] - origin));
        }

        U.boundaryField()[patchi] == pfld;
    }
}


// ************************************************************************* //
