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

#include "solidSolver.H"
#include "volFields.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidSolver, 0);
    defineRunTimeSelectionTable(solidSolver, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidSolver::calcGlobalFaceZones() const
{
    // Find global face zones
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidSolver::calcGlobalFaceZones() const"
        )
            << "Global face zones already fonud"
                << abort(FatalError);
    }

    SLList<label> globalFaceZonesSet;

    const faceZoneMesh& faceZones = mesh().faceZones();

    forAll(faceZones, zoneI)
    {
        const faceZone& curFaceZone = faceZones[zoneI];

        bool globalFaceZone = false;

        forAll(curFaceZone, faceI)
        {
            // if unused face exist
            if (curFaceZone[faceI] >= mesh().nFaces())
            {
//                 globalFaceZonesSet.insert(zoneI);
                globalFaceZone = true;
                break;
            }
        }

        reduce(globalFaceZone, orOp<bool>());

        if (globalFaceZone)
        {
            globalFaceZonesSet.insert(zoneI);
        }
    }

    globalFaceZonesPtr_ = new labelList(globalFaceZonesSet);
}


void Foam::solidSolver::calcGlobalToLocalFaceZonePointMap() const
{
    // Find global face zones
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidSolver::calcGlobalToLocalFaceZonePointMap() const"
        )
            << "Global to local face zones point map already exists"
                << abort(FatalError);
    }

    globalToLocalFaceZonePointMapPtr_ =
        new labelListList(globalFaceZones().size());

    labelListList& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    forAll(globalFaceZones(), zoneI)
    {
        label curZoneID = globalFaceZones()[zoneI];

        labelList curMap(mesh().faceZones()[curZoneID]().nPoints(), -1);

        vectorField fzGlobalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        //- set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints *= 0.0;
        }

        //- pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        //- now every proc has the master's list of FZ points
        //- every proc must now find the mapping from their local FZ points to
        //- the global FZ points

        const vectorField& fzLocalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        const edgeList& fzLocalEdges =
            mesh().faceZones()[curZoneID]().edges();

        const labelListList& fzPointEdges =
            mesh().faceZones()[curZoneID]().pointEdges();

        scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

        forAll(minEdgeLength, pI)
        {
            const labelList& curPointEdges = fzPointEdges[pI];

            forAll(curPointEdges, eI)
            {
                scalar Le = fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                if (Le < minEdgeLength[pI])
                {
                    minEdgeLength[pI] = Le;
                }
            }
        }

        forAll(fzGlobalPoints, globalPointI)
        {
            boolList visited(fzLocalPoints.size(), false);

            forAll(fzLocalPoints, procPointI)
            {
                if (!visited[procPointI])
                {
                    visited[procPointI] = true;

                    label nextPoint = procPointI;

                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[nextPoint]
                          - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[nextPoint])
                    {
                        curMap[globalPointI] = nextPoint;
                        break;
                    }

                    label found = false;

                    while (nextPoint != -1)
                    {
                        const labelList& nextPointEdges =
                            fzPointEdges[nextPoint];

                        scalar minDist = GREAT;
                        label index = -1;
                        forAll(nextPointEdges, edgeI)
                        {
                            label curNgbPoint =
                                fzLocalEdges[nextPointEdges[edgeI]]
                               .otherVertex(nextPoint);

                            if (!visited[curNgbPoint])
                            {
                                visited[curNgbPoint] = true;

                                scalar curDist =
                                    mag
                                    (
                                        fzLocalPoints[curNgbPoint]
                                      - fzGlobalPoints[globalPointI]
                                    );

                                if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                                {
                                    curMap[globalPointI] = curNgbPoint;
                                    found = true;
                                    break;
                                }
                                else if (curDist < minDist)
                                {
                                    minDist = curDist;
                                    index = curNgbPoint;
                                }
                            }
                        }

                        nextPoint = index;
                    }

                    if (found)
                    {
                        break;
                    }
                }
            }

//             forAll(fzLocalPoints, procPointI)
//             {
//                 scalar curDist =
//                     mag
//                     (
//                         fzLocalPoints[procPointI]
//                       - fzGlobalPoints[globalPointI]
//                     );

//                 if (curDist < 1e-4*minEdgeLength[procPointI])
//                 {
//                     curMap[globalPointI] = procPointI;
//                     break;
//                 }
//             }
        }

        forAll(curMap, globalPointI)
        {
            if (curMap[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "solidSolver::calcGlobalToLocalFaceZonePointMap()"
                )
                    << "local to global face zone point map is not correct"
                        << abort(FatalError);
            }
        }

        globalToLocalFaceZonePointMap[zoneI] = curMap;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidSolver::solidSolver
(
    const word& type,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE // must be AUTO_WRITE
        )
    ),
    mesh_(mesh),
    solidProperties_(subDict(type + "Coeffs")),
    bodyForce_
    (
        IOobject
        (
            "bodyForce",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL),
    nCorrectors_(),
    interpFaceGrad_(false),
    lsCellGrad_(false)
{
    if (solidProperties().found("interpFaceGrad"))
    {
        interpFaceGrad_ =
            Switch(solidProperties().lookup("interpFaceGrad"));

        if (interpFaceGrad_)
        {
            Info << "Using interpolation to calculate face "
                << "gradient of displacement" << endl;
        }
    }

    if (solidProperties().found("lsCellGrad"))
    {
        lsCellGrad_ =
            Switch(solidProperties().lookup("lsCellGrad"));

        if (lsCellGrad_)
        {
            Info << "Using least squares method to calculate cell "
                << "gradient of displacement" << endl;
        }
    }

    // Disabling the 3rd direction for axisymmetric cases
    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp =
                refCast<const wedgePolyPatch>
                (
                    mesh_.boundaryMesh()[patchI]
                );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        Info<< nl << "Axisymmetric case: disabling the 3rd direction"
            << nl << endl;

        // We will const_cast as it is currently the tidiest way,
        // until polyMesh is modified or gives write access to solutionD
        Vector<label>& solD = const_cast<Vector<label>&>(mesh_.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        Info << wedgeDirVec << endl;

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;
            }
            else
            {
                solD[cmpt] = 1;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidSolver::~solidSolver()
{
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidSolver::hydPressure() const
{
    notImplemented("tmp<volScalarField> Foam::solidSolver::hydPressure()");

    tmp<volScalarField> tHydPressure
    (
        new volScalarField
        (
            IOobject
            (
                "hydPressure",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0)
        )
    );

    return tHydPressure;
}

bool Foam::solidSolver::thermalStress() const
{
    return mesh().objectRegistry::found("thermalProperties");
}

const Foam::labelList& Foam::solidSolver::globalFaceZones() const
{
    if (!globalFaceZonesPtr_)
    {
        calcGlobalFaceZones();
    }

    return *globalFaceZonesPtr_;
}

const Foam::labelListList&
Foam::solidSolver::globalToLocalFaceZonePointMap() const
{
    if (!globalToLocalFaceZonePointMapPtr_)
    {
        calcGlobalToLocalFaceZonePointMap();
    }

    return *globalToLocalFaceZonePointMapPtr_;
}

bool Foam::solidSolver::read()
{
    if (regIOobject::read())
    {
        solidProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


//- Patch face acceleration
Foam::tmp<Foam::vectorField> Foam::solidSolver::patchFaceAcceleration
(
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().boundary()[patchID].size(),
            vector::zero
        )
    );

    return tAcceleration;
}


//- Patch face velocity
Foam::tmp<Foam::vectorField> Foam::solidSolver::faceZoneAcceleration
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );

    return tAcceleration;
}


//- Patch face velocity
Foam::tmp<Foam::vectorField> Foam::solidSolver::faceZoneVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tVelocity
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );

    return tVelocity;
}


// ************************************************************************* //
