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

Class
    solidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::solidContactFvPatchVectorField::calcZoneIndex() const
{
    if (zoneIndexPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcZoneIndex() const"
        )   << "zoneIndexPtr_ already set" << abort(FatalError);
    }

    zoneIndexPtr_ = new label(-1);
    label& zoneIndex = *zoneIndexPtr_;

    word zoneName = patch().name() + "FaceZone";

    faceZoneID zone(zoneName, patch().boundaryMesh().mesh().faceZones());

    if (!zone.active())
    {
        FatalErrorIn("solidContactFvPatchScalarField")
            << "Face zone name " << zoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    zoneIndex = zone.index();
}


Foam::label Foam::solidContactFvPatchVectorField::zoneIndex() const
{
    if (!zoneIndexPtr_)
    {
        calcZoneIndex();
    }

    return *zoneIndexPtr_;
}


void Foam::solidContactFvPatchVectorField::calcShadowZoneIndex() const
{
    if (shadowZoneIndexPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcShadowZoneIndex() const"
        )   << "shadowZoneIndexPtr_ already set" << abort(FatalError);
    }

    shadowZoneIndexPtr_ = new label(-1);
    label& shadowZoneIndex = *shadowZoneIndexPtr_;

    word shadowZoneName = shadowPatchName_ + "FaceZone";

    faceZoneID shadowZone
    (
        shadowZoneName, patch().boundaryMesh().mesh().faceZones()
    );

    if (!shadowZone.active())
    {
        FatalErrorIn("solidContactFvPatchScalarField")
            << "Face zone name " << shadowZoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    shadowZoneIndex = shadowZone.index();
}


Foam::label Foam::solidContactFvPatchVectorField::shadowZoneIndex() const
{
    if (!shadowZoneIndexPtr_)
    {
        calcShadowZoneIndex();
    }

    return *shadowZoneIndexPtr_;
}


void Foam::solidContactFvPatchVectorField::calcNormalModel() const
{
    if (normalModelPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "normalModelPtr_ already set" << abort(FatalError);
    }

    normalModelPtr_ =
        normalContactModel::New
        (
            word(dict().lookup("normalContactModel")),
            patch(),
            dict(),
            patch().index(), // master
            shadowPatchIndex(), // slave
            zoneIndex(), // master face zone ID
            shadowZoneIndex(), // slave face zone ID
            zone(),
            shadowZone()
        ).ptr();
}


Foam::normalContactModel&
Foam::solidContactFvPatchVectorField::normalModel()
{
    if (master())
    {
        if (!normalModelPtr_)
        {
            calcNormalModel();
        }

        return *normalModelPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.normalModel();
}


const Foam::normalContactModel&
Foam::solidContactFvPatchVectorField::normalModel() const
{
    if (master())
    {
        if (!normalModelPtr_)
        {
            calcNormalModel();
        }

        return *normalModelPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.normalModel();
}


void Foam::solidContactFvPatchVectorField::calcFrictionModel() const
{
    if (frictionModelPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcFrictionModel() const"
        )   << "frictionModelPtr_ already set" << abort(FatalError);
    }

    frictionModelPtr_ =
        frictionContactModel::New
        (
            word(dict().lookup("frictionContactModel")),
            patch(),
            dict(),
            patch().index(), // master
            shadowPatchIndex(), // slave
            zoneIndex(), // master face zone ID
            shadowZoneIndex() // slave face zone ID
        ).ptr();
}


Foam::frictionContactModel&
Foam::solidContactFvPatchVectorField::frictionModel()
{
    if (master())
    {
        if (!frictionModelPtr_)
        {
            //Info<< "before" << flush
            //    << word(dict().subDict(name+"FrictionModelDict")) << endl;
            calcFrictionModel();
        }

        return *frictionModelPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.frictionModel();
}


const Foam::frictionContactModel&
Foam::solidContactFvPatchVectorField::frictionModel() const
{
    if (master())
    {
        if (!frictionModelPtr_)
        {
            calcFrictionModel();
        }

        return *frictionModelPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.frictionModel();
}


void Foam::solidContactFvPatchVectorField::calcAllPointsDeformed() const
{
    if (allPointsDeformedPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcAllPointsDeformed() const"
        )   << "allPointsDeformedPtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    allPointsDeformedPtr_ = new vectorField(mesh.allPoints());
}


Foam::vectorField&
Foam::solidContactFvPatchVectorField::allPointsDeformed()
{
    if (master())
    {
        if (!allPointsDeformedPtr_)
        {
            calcAllPointsDeformed();
        }

        return *allPointsDeformedPtr_;
    }

    // slave
    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.allPointsDeformed();
}


const Foam::vectorField&
Foam::solidContactFvPatchVectorField::allPointsDeformed() const
{
    if (master())
    {
        if (!allPointsDeformedPtr_)
        {
            calcAllPointsDeformed();
        }

        return *allPointsDeformedPtr_;
    }

    // slave
    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.allPointsDeformed();
}


void Foam::solidContactFvPatchVectorField::updateAllPointsDeformed()
{
    if (master())
    {
        vectorField& p = allPointsDeformed();

        // lookup displacement field
        const pointVectorField& field =
            this->db().objectRegistry::lookupObject<pointVectorField>
            (
                "point" + this->dimensionedInternalField().name()
            );

        // total field is used for non-moving mesh incremental solvers
        bool incrementalNonMovingMesh = false;
        if (fieldName_ == "DD" && nonLinear_ == nonLinearGeometry::OFF)
        {
            incrementalNonMovingMesh = true;
        }

        const pointVectorField& totalField =
            this->db().objectRegistry::lookupObject<pointVectorField>
            (
                "pointD" //"pointU" // total displacement point field
            ).oldTime();

        // pointDU does not include unused points from the global face zones
        const label nUsedPoints = field.size();

        // Move zone points

        const vectorField& oldZoneP = oldZonePoints();
        const vectorField& oldShadowZoneP = oldShadowZonePoints();

        // Zone may not have been created so we will use the zones in the mesh
        //const labelList& zoneMeshPoints = zone().meshPoints();
        //const labelList& shadowZoneMeshPoints = shadowZone().meshPoints();
        const fvMesh& mesh = patch().boundaryMesh().mesh();
        const labelList& zoneMeshPoints =
            mesh.faceZones()[zoneIndex()]().meshPoints();
        const labelList& shadowZoneMeshPoints =
            mesh.faceZones()[shadowZoneIndex()]().meshPoints();

        forAll(zoneMeshPoints, pI)
        {
            const label pointID = zoneMeshPoints[pI];

            p[pointID] = oldZoneP[pI];

            if (pointID < nUsedPoints)
            {
                p[pointID] += field[pointID];

                if (incrementalNonMovingMesh)
                {
                    p[pointID] += totalField[pointID];
                }
            }
        }

        forAll(shadowZoneMeshPoints, pI)
        {
            const label pointID = shadowZoneMeshPoints[pI];

            p[pointID] = oldShadowZoneP[pI];

            if (pointID < nUsedPoints)
            {
                p[pointID] += field[pointID];

                if (incrementalNonMovingMesh)
                {
                    p[pointID] += totalField[pointID];
                }
            }
        }

        // Set values for unused points

        const vectorField& pointFieldI = field.internalField();
        const vectorField& pointTotalFieldI = field.internalField();

        const labelList& gFaceZones = globalFaceZones();

        forAll(gFaceZones, zoneI)
        {
            const label curZoneID = gFaceZones[zoneI];

            const labelList& curMap =
                globalToLocalFaceZonePointMap()[zoneI];

            const labelList& curZoneMeshPoints =
                mesh.faceZones()[curZoneID]().meshPoints();

            vectorField curGlobalZonePointDispl
                (
                    curZoneMeshPoints.size(),
                    vector::zero
                );

            //-Inter-proc points are shared by multiple procs
            // pointNumProc is the number of procs which a point lies on
            scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];

                if(curZoneMeshPoints[localPoint] < mesh.nPoints())
                {
                    label procPoint = curZoneMeshPoints[localPoint];

                    curGlobalZonePointDispl[globalPointI] =
                        pointFieldI[procPoint];

                    if (incrementalNonMovingMesh)
                    {
                        curGlobalZonePointDispl[globalPointI] +=
                            pointTotalFieldI[procPoint];
                    }

                    pointNumProcs[globalPointI] = 1;
                }
            }

            if (Pstream::parRun())
            {
                reduce(curGlobalZonePointDispl, sumOp<vectorField>());
                reduce(pointNumProcs, sumOp<scalarField>());

                //- now average the displacement between all procs
                curGlobalZonePointDispl /= pointNumProcs;
            }

            //- The curZonePointsDisplGlobal now contains the correct face zone
            //  displacement in a global master processor order, now convert
            //  them back into the local proc order

            vectorField curZonePointDispl
            (
                curZoneMeshPoints.size(),
                vector::zero
            );

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];

                curZonePointDispl[localPoint] =
                    curGlobalZonePointDispl[globalPointI];
            }

            forAll(curZonePointDispl, pointI)
            {
                // unused points
                if (curZoneMeshPoints[pointI] >= mesh.nPoints())
                {
                    // Pout<< "correcting motion for point "
                    //     << curZoneMeshPoints[pointI] << endl;
                    p[curZoneMeshPoints[pointI]] += curZonePointDispl[pointI];
                }
            }
        }
    }
    else
    {
        const volVectorField& field =
            this->db().objectRegistry::lookupObject<volVectorField>
            (
                this->dimensionedInternalField().name()
            );

        // We will const cast the shadow patch to update the master points

        solidContactFvPatchVectorField& shadowPatchField =
            const_cast<solidContactFvPatchVectorField&>
            (
                refCast<const solidContactFvPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndex_]
                )
            );

        return shadowPatchField.updateAllPointsDeformed();
    }
}


const Foam::vectorField&
Foam::solidContactFvPatchVectorField::oldZonePoints() const
{
    if (master())
    {
        if (!oldZonePointsPtr_)
        {
            calcOldZonePoints();
        }

        return *oldZonePointsPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.oldZonePoints();
}


void Foam::solidContactFvPatchVectorField::calcOldZonePoints() const
{
    if (oldZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcOldZonePoints() const"
        )   << "pointer already set"
            << abort(FatalError);
    }

    if
    (
        nonLinear_ != nonLinearGeometry::OFF
     && nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
     && nonLinear_ != nonLinearGeometry::TOTAL_LAGRANGIAN
    )
    {
        FatalErrorIn
        (
            "void vectorField& Foam::solidContactFvPatchVectorField::"
            "oldZonePoints() const"
        )   << "only implemented for nonLinear off, "
            << "updatedLagrangianKirchhoff, and totalLagrangian"
            << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    oldZonePointsPtr_ =
        new vectorField(mesh.faceZones()[zoneIndex()]().localPoints());
}


const Foam::vectorField&
Foam::solidContactFvPatchVectorField::oldShadowZonePoints() const
{
    if (master())
    {
        if (!oldShadowZonePointsPtr_)
        {
            calcOldShadowZonePoints();
        }

        return *oldShadowZonePointsPtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.oldZonePoints();
}


void Foam::solidContactFvPatchVectorField::calcOldShadowZonePoints() const
{
    if (oldShadowZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcOldShadowZonePoints() const"
        )   << "pointer already set"
            << abort(FatalError);
    }

    if
    (
        nonLinear_ != nonLinearGeometry::OFF
     && nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
     && nonLinear_ != nonLinearGeometry::TOTAL_LAGRANGIAN
    )
    {
        FatalErrorIn
        (
            "void vectorField& Foam::solidContactFvPatchVectorField::"
            "oldShadowZonePoints() const"
        )   << "only implemented for nonLinear off, "
            << "updatedLagrangianKirchhoff, and totalLagrangian"
            << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    oldShadowZonePointsPtr_ =
        new vectorField(mesh.faceZones()[shadowZoneIndex()]().localPoints());
}


void Foam::solidContactFvPatchVectorField::calcZone() const
{
    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcZone() const"
        )   << "zonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    zonePtr_ =
        new primitiveFacePatch
        (
            faceList(mesh.faceZones()[zoneIndex()].size()),
            allPointsDeformed()
        );

    primitiveFacePatch& patch = *zonePtr_;

    const faceList& f = mesh.allFaces();

    const labelList& addr = mesh.faceZones()[zoneIndex()];
    const boolList& flip = mesh.faceZones()[zoneIndex()].flipMap();

    forAll (addr, faceI)
    {
        if (flip[faceI])
        {
            patch[faceI] = f[addr[faceI]].reverseFace();
        }
        else
        {
            patch[faceI] = f[addr[faceI]];
        }
    }
}


const Foam::primitiveFacePatch&
Foam::solidContactFvPatchVectorField::zone() const
{
    if (master())
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.zone();
}


Foam::primitiveFacePatch& Foam::solidContactFvPatchVectorField::zone()
{
    if (master())
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.zone();
}


void Foam::solidContactFvPatchVectorField::calcShadowZone() const
{
    if (shadowZonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcShadowZone() const"
        )   << "shadowZonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    shadowZonePtr_ =
        new primitiveFacePatch
        (
            faceList(mesh.faceZones()[shadowZoneIndex()].size()),
            allPointsDeformed()
        );

    primitiveFacePatch& patch = *shadowZonePtr_;

    const faceList& f = mesh.allFaces();

    const labelList& addr = mesh.faceZones()[shadowZoneIndex()];
    const boolList& flip = mesh.faceZones()[shadowZoneIndex()].flipMap();

    forAll (addr, faceI)
    {
        if (flip[faceI])
        {
            patch[faceI] = f[addr[faceI]].reverseFace();
        }
        else
        {
            patch[faceI] = f[addr[faceI]];
        }
    }
}


const Foam::primitiveFacePatch&
Foam::solidContactFvPatchVectorField::shadowZone() const
{
    if (master())
    {
        if (!shadowZonePtr_)
        {
            calcShadowZone();
        }

        return *shadowZonePtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.shadowZone();
}


Foam::primitiveFacePatch& Foam::solidContactFvPatchVectorField::shadowZone()
{
    if (master())
    {
        if (!shadowZonePtr_)
        {
            calcShadowZone();
        }

        return *shadowZonePtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.shadowZone();
}


void Foam::solidContactFvPatchVectorField::calcZoneToZone() const
{
    // Create zone-to-zone interpolation
    if (zoneToZonePtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::calcZoneToZone() const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    // Check master and slave patch
    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    if (master())
    {
        if (shadowPatchField.master() == true)
        {
            FatalErrorIn("solidContactFvPatchScalarField")
                << "There are two master patches"
                << abort(FatalError);
        }
    }
    else
    {
        if (shadowPatchField.master() == false)
        {
            FatalErrorIn("solidContactFvPatchScalarField")
                << "There is no master patch"
                << abort(FatalError);
        }
    }

    if (master())
    {
        // Create interpolation for patches
        zoneToZonePtr_ =
            new ggiZoneInterpolation
            (
                zone(),
                shadowZone(),
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                0,              // Non-overlapping face tolerances
                0,              // HJ, 24/Oct/2008
                true,           // Rescale weighting factors.  Bug fix, MB.
                ggiInterpolation::AABB
                //ggiInterpolation::BB_OCTREE  // Octree search, MB.
                //extendedGgiInterpolation::BB_OCTREE
            );
    }
    else
    {
        FatalErrorIn
        (
            "void solidContactFvPatchVectorField::calcZoneToZone() const"
        )   << "Attempting to create GGIInterpolation on a slave"
            << abort(FatalError);
    }
}


const Foam::ggiZoneInterpolation&
Foam::solidContactFvPatchVectorField::zoneToZone() const
{
    if (master())
    {
        if (!zoneToZonePtr_)
        {
            if (debug)
            {
                word zoneName =
                    patch().boundaryMesh().mesh().faceZones()
                    [
                        zoneIndex()
                    ].name();

                word shadowZoneName =
                    patch().boundaryMesh().mesh()
                    .faceZones()[shadowZoneIndex()].name();

                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow zones: "
                    << zoneName << "/" << shadowZoneName
                    << endl;
            }

            calcZoneToZone();
        }

        return *zoneToZonePtr_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.zoneToZone();
}


Foam::label Foam::solidContactFvPatchVectorField::correctionFreq() const
{
    if (master())
    {
        return correctionFreq_;
    }

    // slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.correctionFreq();
}


Foam::ggiZoneInterpolation&
Foam::solidContactFvPatchVectorField::zoneToZone()
{
    if (master())
    {
        if (!zoneToZonePtr_)
        {
            word zoneName =
                patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

            word shadowZoneName =
                patch().boundaryMesh().mesh()
               .faceZones()[shadowZoneIndex()].name();

            if (debug)
            {
                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow zones: "
                    << zoneName << "/" << shadowZoneName
                    << endl;
            }

            calcZoneToZone();
        }

        return *zoneToZonePtr_;
    }

    // Slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.zoneToZone();
}


void Foam::solidContactFvPatchVectorField::calcSlaveFaceNormals() const
{
    // Create slave face normals
    if (slaveFaceNormalsPtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::calcSlaveFaceNormals() const"
        )   << "slaveFaceNormals pointer already calculated"
            << abort(FatalError);
    }

    slaveFaceNormalsPtr_ =
        new vectorField
        (
            patch().boundaryMesh().mesh().boundary()[shadowPatchIndex()].size(),
            vector::zero
        );
}


Foam::vectorField&
Foam::solidContactFvPatchVectorField::slaveFaceNormals()
{
    if (master())
    {
        if (!slaveFaceNormalsPtr_)
        {
            calcSlaveFaceNormals();
        }

        return *slaveFaceNormalsPtr_;
    }

    // Slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.slaveFaceNormals();
}


void Foam::solidContactFvPatchVectorField::calcGlobalFaceZones() const
{
    // Create slave face normals
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::"
            "calcGlobalFaceZones() const"
        )   << "pointer already calculated"
            << abort(FatalError);
    }

    globalFaceZonesPtr_ = new labelList(0);

    labelList& globalFaceZones = *globalFaceZonesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (Pstream::parRun())
    {
        SLList<label> globalFaceZonesSet;

        // Lookup globalFaceZones from decomposeParDict
        IOdictionary decompDict
            (
                IOobject
                (
                    "decomposeParDict",
                    mesh.time().time().system(),
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        if (decompDict.found("globalFaceZones"))
        {
            wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

            const faceZoneMesh& faceZones = mesh.faceZones();

            forAll(globalFaceZoneNames, nameI)
            {
                const label zoneID =
                    faceZones.findZoneID(globalFaceZoneNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn
                    (
                        "solidContactFvPatchVectorField::"
                        "solidContactFvPatchVectorField"
                    )   << "cannot find globalFaceZone:"
                        << " " << globalFaceZoneNames[nameI]
                        << abort(FatalError);
                }

                globalFaceZonesSet.insert(zoneID);
            }
        }

        globalFaceZones = labelList(globalFaceZonesSet);
    }
}


const Foam::labelList&
Foam::solidContactFvPatchVectorField::globalFaceZones() const
{
    if (master())
    {
        if (!globalFaceZonesPtr_)
        {
            calcGlobalFaceZones();
        }

        return *globalFaceZonesPtr_;
    }

    // Slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.globalFaceZones();
}


void Foam::solidContactFvPatchVectorField::
calcGlobalToLocalFaceZonePointMap() const
{
    // Create slave face normals
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::"
            "calcGlobalToLocalFaceZonePointMap() const"
        )   << "slaveFaceNormals pointer already calculated"
            << abort(FatalError);
    }

    const labelList& gFaceZones = globalFaceZones();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    fileName mapName = "globalToLocalFaceZonePointMapSolidContact";

    word timeName = mesh.time().timeName();

    bool mapExists = false;

    {
        IOobject mapHeader
        (
            mapName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            mapExists = true;
        }
        else
        {
            // Check previous time-step
            instantList times = mesh.time().times();

            if (times.size() > 1)
            {
                word prevTimeName = times[times.size() - 1].name();

                Pout<< "Reading face map from " << prevTimeName << endl;

                IOobject mapPrevHeader
                (
                    mapName,
                    prevTimeName,
                    mesh,
                    IOobject::MUST_READ
                );

                if (mapPrevHeader.headerOk())
                {
                    mapExists = true;
                    timeName = prevTimeName;
                }
            }
        }
    }

    globalToLocalFaceZonePointMapPtr_ =
        new IOList<labelList>
        (
            IOobject
            (
                mapName,
                timeName,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            gFaceZones.size()
        );

    IOList<labelList>& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    // Calculate map if it has not been read
    if (!mapExists && Pstream::parRun())
    {
        forAll(gFaceZones, zoneI)
        {
            label curZoneID = gFaceZones[zoneI];

            Info<< "Creating faceMap for globalFaceZones "
                << mesh.faceZones()[curZoneID].name()<< endl;

            labelList curMap(mesh.faceZones()[curZoneID]().nPoints(), -1);

            vectorField fzGlobalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            // Set all slave points to zero because only the master order is
            // used
            if (!Pstream::master())
            {
                fzGlobalPoints = vector::zero;
            }

            //- pass points to all procs
            reduce(fzGlobalPoints, sumOp<vectorField>());

            // Now every proc has the master's list of FZ points
            // every proc must now find the mapping from their local FZ points
            // to the global FZ points

            const vectorField& fzLocalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            const edgeList& fzLocalEdges =
                mesh.faceZones()[curZoneID]().edges();

            const labelListList& fzPointEdges =
                mesh.faceZones()[curZoneID]().pointEdges();

            scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

            forAll(minEdgeLength, pI)
            {
                const labelList& curPointEdges = fzPointEdges[pI];

                forAll(curPointEdges, eI)
                {
                    scalar Le =
                        fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                    if (Le < minEdgeLength[pI])
                    {
                        minEdgeLength[pI] = Le;
                    }
                }
            }

            forAll(fzGlobalPoints, globalPointI)
            {
                //scalar minDist = GREAT;
                bool pointFound = false;

                forAll(fzLocalPoints, procPointI)
                {
                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[procPointI]
                            - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[procPointI])
                    {
                        curMap[globalPointI] = procPointI;

                        if (pointFound)
                        {
                            FatalError
                                << "findGlobalFaceZones: point found twice!"
                                << abort(FatalError);
                        }
                        pointFound = true;
                    }
                }
            }

            forAll(curMap, globalPointI)
            {
                if (curMap[globalPointI] == -1)
                {
                    FatalErrorIn
                    (
                        "solidContactFvPatchVectorField::"
                        "calcGlobalToLocalFaceZonePointMap()"
                    )   << "local to global face zone point map is not correct"
                        << " for zone " << zoneI
                        << abort(FatalError);
                }
            }

            globalToLocalFaceZonePointMap[zoneI] = curMap;
        }
    }
    else
    {
        Info<< "globalToLocalFaceZonePointMap read from file" << endl;
    }
}


const Foam::IOList<Foam::labelList>&
Foam::solidContactFvPatchVectorField::globalToLocalFaceZonePointMap() const
{
    if (master())
    {
        if (!globalToLocalFaceZonePointMapPtr_)
        {
            calcGlobalToLocalFaceZonePointMap();
        }

        return *globalToLocalFaceZonePointMapPtr_;
    }

    // Slave

    const volVectorField& field =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.globalToLocalFaceZonePointMap();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    master_("undefined"),
    shadowPatchName_("undefined"),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    contactActive_(false),
    rigidMaster_(false),
    dict_(NULL),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    allPointsDeformedPtr_(NULL),
    oldZonePointsPtr_(NULL),
    oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    iCorr_(0),
    correctionFreq_(0),
    orthotropic_(false),
    forceCorrection_(false),
    forceCorrectionFriction_(false),
    slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL),
    contactResidual_(1.0),
    nonLinear_(nonLinearGeometry::OFF)
{}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    contactActive_(ptf.contactActive_),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    allPointsDeformedPtr_(NULL),
    oldZonePointsPtr_(NULL),
    oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    iCorr_(ptf.iCorr_),
    correctionFreq_(ptf.correctionFreq_),
    orthotropic_(ptf.orthotropic_),
    forceCorrection_(ptf.forceCorrection_),
    forceCorrectionFriction_(ptf.forceCorrectionFriction_),
    slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL),
    contactResidual_(ptf.contactResidual_),
    nonLinear_(ptf.nonLinear_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.allPointsDeformedPtr_)
    {
        allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    }

    if (ptf.oldZonePointsPtr_)
    {
        oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    }

    if (ptf.oldShadowZonePointsPtr_)
    {
        oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    if (ptf.slaveFaceNormalsPtr_)
    {
        slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    master_(dict.lookup("master")),
    shadowPatchName_(dict.lookup("shadowPatch")),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    contactActive_(true), // hard-code true
    rigidMaster_(false),
    dict_(dict),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    allPointsDeformedPtr_(NULL),
    oldZonePointsPtr_(NULL),
    oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    iCorr_(0),
    correctionFreq_(1), // hard-code to 1
    orthotropic_(false),
    forceCorrection_(false),
    forceCorrectionFriction_(false),
    slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL),
    contactResidual_(1.0),
    // We should probably look nonLinear up just before using it as it may
    // change during the simulation
    nonLinear_
    (
        nonLinearGeometry::nonLinearNames_.read
        (
            patch().boundaryMesh().mesh().solutionDict().subDict
            (
                "solidMechanics"
            ).lookup("nonLinear")
        )
    )
{
    Info<< "Creating " << solidContactFvPatchVectorField::typeName << " patch"
        << endl;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Shadow patch index
    polyPatchID shadow(shadowPatchName_, mesh.boundaryMesh());

    if (!shadow.active())
    {
        FatalErrorIn("solidContactFvPatchScalarField")
            << "Shadow patch name " << shadowPatchName_ << " not found."
            << abort(FatalError);
    }

    shadowPatchIndex_ = shadow.index();

    // Master creates contact laws
    if (master_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));
    }

    if (dict.found("refValue"))
    {
        Info<< "    Reading refValue allowing restart of case" << endl;
        this->refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = vector::zero;
    }

    if (dict.found("refGradient"))
    {
        Info<< "    Reading refGrad allowing restart of case" << endl;
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    if (dict.found("valueFraction"))
    {
        Info<< "    Reading valueFraction allowing restart of case" << endl;
        this->valueFraction() =
            symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        this->valueFraction() = symmTensor::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector> normalValue = transform(valueFraction(), refValue());

        Field<vector> gradValue =
            this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

        Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    contactActive_(ptf.contactActive_),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(ptf.normalModelPtr_),
    frictionModelPtr_(ptf.frictionModelPtr_),
    allPointsDeformedPtr_(ptf.allPointsDeformedPtr_),
    oldZonePointsPtr_(ptf.oldZonePointsPtr_),
    oldShadowZonePointsPtr_(ptf.oldShadowZonePointsPtr_),
    zonePtr_(ptf.zonePtr_),
    shadowZonePtr_(ptf.shadowZonePtr_),
    zoneToZonePtr_(ptf.zoneToZonePtr_),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    iCorr_(ptf.iCorr_),
    correctionFreq_(ptf.correctionFreq_),
    orthotropic_(ptf.orthotropic_),
    forceCorrection_(ptf.forceCorrection_),
    forceCorrectionFriction_(ptf.forceCorrectionFriction_),
    slaveFaceNormalsPtr_(ptf.slaveFaceNormalsPtr_),
    globalFaceZonesPtr_(ptf.globalFaceZonesPtr_),
    globalToLocalFaceZonePointMapPtr_(ptf.globalToLocalFaceZonePointMapPtr_),
    contactResidual_(ptf.contactResidual_),
    nonLinear_(ptf.nonLinear_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.allPointsDeformedPtr_)
    {
        allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    }

    if (ptf.oldZonePointsPtr_)
    {
        oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    }

    if (ptf.oldShadowZonePointsPtr_)
    {
        oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    if (ptf.slaveFaceNormalsPtr_)
    {
        slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    contactActive_(ptf.contactActive_),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(ptf.normalModelPtr_),
    frictionModelPtr_(ptf.frictionModelPtr_),
    // allPointsDeformedPtr_(ptf.allPointsDeformedPtr_),
    // oldZonePointsPtr_(ptf.oldZonePointsPtr_),
    // oldShadowZonePointsPtr_(ptf.oldShadowZonePointsPtr_),
    // zonePtr_(ptf.zonePtr_),
    // shadowZonePtr_(ptf.shadowZonePtr_),
    // zoneToZonePtr_(ptf.zoneToZonePtr_),
    //normalModelPtr_(NULL),
    //frictionModelPtr_(NULL),
    allPointsDeformedPtr_(NULL),
    oldZonePointsPtr_(NULL),
    oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    iCorr_(ptf.iCorr_),
    correctionFreq_(ptf.correctionFreq_),
    orthotropic_(ptf.orthotropic_),
    forceCorrection_(ptf.forceCorrection_),
    forceCorrectionFriction_(ptf.forceCorrectionFriction_),
    // slaveFaceNormalsPtr_(ptf.slaveFaceNormalsPtr_),
    slaveFaceNormalsPtr_(NULL),
    //globalFaceZonesPtr_(ptf.globalFaceZonesPtr_),
    globalFaceZonesPtr_(NULL),
    //globalToLocalFaceZonePointMapPtr_(ptf.globalToLocalFaceZonePointMapPtr_),
    globalToLocalFaceZonePointMapPtr_(NULL),
    contactResidual_(ptf.contactResidual_),
    nonLinear_(ptf.nonLinear_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.allPointsDeformedPtr_)
    {
        allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    }

    if (ptf.oldZonePointsPtr_)
    {
        oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    }

    if (ptf.oldShadowZonePointsPtr_)
    {
        oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    if (ptf.slaveFaceNormalsPtr_)
    {
        slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::~solidContactFvPatchVectorField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::solidContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);

    // Force all data to be re-created when needed
    clearOut();
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::solidContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void Foam::solidContactFvPatchVectorField::clearOut()
{
    deleteDemandDrivenData(zoneIndexPtr_);
    deleteDemandDrivenData(shadowZoneIndexPtr_);
    deleteDemandDrivenData(normalModelPtr_);
    deleteDemandDrivenData(frictionModelPtr_);
    deleteDemandDrivenData(allPointsDeformedPtr_);
    deleteDemandDrivenData(oldZonePointsPtr_);
    deleteDemandDrivenData(oldShadowZonePointsPtr_);
    deleteDemandDrivenData(zonePtr_);
    deleteDemandDrivenData(shadowZonePtr_);
    deleteDemandDrivenData(zoneToZonePtr_);
    deleteDemandDrivenData(slaveFaceNormalsPtr_);
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


void Foam::solidContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (contactActive_)
    {
        // if (!normalModelPtr_ || !frictionModelPtr_)
        // {
        //     if (!master_)
        //     {
        //         const volVectorField& Ufield =
        //             this->db().objectRegistry::lookupObject<volVectorField>
        //             (
        //                 fieldName_
        //             );

        //         const solidContactFvPatchVectorField& Upatch =
        //             refCast<const solidContactFvPatchVectorField>
        //             (
        //                 Ufield.boundaryField()[shadowPatchIndex()]
        //             );
                // Info<< "    Slave contact patch " << patch().name()
                //     << " grabbing normalModel pointer from master"
                //     << endl;
                // normalModelPtr_ = Upatch.normalModelPtr();
                // if (!normalModelPtr_)
                // {
                //     FatalError
                //         << "\nThe patch " << patch().name()
                //         << " has a NULL normalModel pointer!"
                //         << abort(FatalError);
                // }

                // Info<< "    Slave contact patch " << patch().name()
                //     << " grabbing frictionModel pointer from master"
                //     << endl;
                // frictionModelPtr_ = Upatch.frictionModelPtr();
                // if (!frictionModelPtr_)
                // {
                //     FatalError
                //         << "\nThe patch " << patch().name()
                //         << " has a NULL frictionModel pointer!"
                //         << abort(FatalError);
                // }

                // lookup master's correction frequency
                //correctionFreq_ = Upatch.correctionFreq();
        //     }
        //     else
        //     {
        //         FatalErrorIn("solidContactFvPatchVectorField::updateCoeff()")
        //             << "NULL contactLaw\ncontactLaw not created by master."
        //             << abort(FatalError);
        //     }
        // }

        // if it is a new time step then reset iCorr
        if (curTimeIndex_ != this->db().time().timeIndex())
        {
            curTimeIndex_ = this->db().time().timeIndex();
            iCorr_ = 0;

            // update old face zone points
            if
            (
                master_ &&
                (
                    nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN
                    || nonLinear_ ==
                        nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
                    || nonLinear_ == nonLinearGeometry::DEFORMED_LAGRANGIAN
              // || nonLinear_ == nonLinearGeometry::DEFORMED_LAGRANGIAN_HUGHES)
                )
            )
            {
                deleteDemandDrivenData(oldZonePointsPtr_);
                deleteDemandDrivenData(oldShadowZonePointsPtr_);
            }
        }


        // Update position of zones and interpolation weights
        updateAllPointsDeformed();
        zoneToZone().movePoints(tensorField(0), tensorField(0), vectorField(0));
        zone().clearGeom();
        shadowZone().clearGeom();


        // the time taken to correct the contact may not be negligible
        // so reduce the correctiion frequency can speed up the simulation
        if
        (
            iCorr_++ % correctionFreq() == 0
            || forceCorrection_
            || forceCorrectionFriction_
        )
        {
            forceCorrection_ = false;

            if (master())
            {
                // Correct laws
                // the normal model sets what face normal to use
                // on the slave e.g. use actual face normals, or master normals
                // interpolated to the slave, undeformed/deformed normals.
                // the friction model then uses these face normals

                normalModel().correct
                    (
                        zone(),
                        shadowZone(),
                        alg_,
                        dir_,
                        fieldName_,
                        orthotropic_,
                        nonLinear_,
                        slaveFaceNormals(),
                        zoneToZone()
                    );

                frictionModel().correct
                    (
                        normalModel().slavePressure(),
                        zone(),
                        shadowZone(),
                        alg_,
                        dir_,
                        "ggi",
                        fieldName_,
                        orthotropic_,
                        nonLinear_,
                        slaveFaceNormals()
                    );

                forceCorrectionFriction_ = false;
            }
        } // if correction frequency
    } // if contactActive

    // Fix: philipc 25-09-13 We always set boundary conditions even when contact
    // is not active as we need to
    // iterative enforce the traction condition
    if (!master_)
    {
        // set refValue, refGrad and valueFraction

        // refValue
        refValue() =
            normalModel().slaveDisp()
            + frictionModel().slaveDisp();

        // refGrad - set traction
        refGrad() =
            tractionBoundaryGradient().snGrad
            (
                frictionModel().slaveTraction()
              + normalModel().slavePressure(),    // surface traction
                scalarField(patch().size(),0.0),  // surface pressure
                fieldName_,                       // working field name
                "U",                              // total field name
                patch(),                          // polyPatch
                bool(fieldName_ == "DU")          // incremental
            );

        // valueFraction
        valueFraction() =
            normalModel().slaveValueFrac()
            + frictionModel().slaveValueFrac();

        // slave calculates contact residual
        contactResidual_ = normalModel().residual();
    }
    else
    {
        // master is always traction condition
        // interpolate traction from slave
        // Dirichlet-Neumann uses lagged traction
        // penalty uses same as for slave traction

        if (rigidMaster_)
        {
            // set to master to traction free if it is rigid
            refGrad() =
                tractionBoundaryGradient().snGrad
                (
                    vectorField(patch().size(),vector::zero), // traction
                    scalarField(patch().size(),0.0),  // surface pressure
                    fieldName_,                       // working field name
                    "U",                              // total field name
                    patch(),                          // polyPatch
                    bool(fieldName_ == "DU")          // incremental
                );
        }
        else
        {
            // Interpolate slave traction to the master

            vectorField slavePatchTraction =
                -frictionModel().slaveTractionForMaster()
                -normalModel().slavePressure();

            vectorField slaveZoneTraction =
                zoneField
                (
                    shadowZoneIndex(),
                    shadowPatchIndex(),
                    slavePatchTraction
                );

            // Face-to-face
            vectorField masterZoneTraction =
                zoneToZone().slaveToMaster(slaveZoneTraction);

            // Face-to-point -- point-to-point -- point-to-face
            // Does not seem to solve problem of "jittering" stresses during
            // sliding
            // vectorField slaveZonePointTraction =
            //     zoneFaceToPointInterpolate
            //     (
            //         shadowZoneIndex(), slaveZoneTraction
            //     );

            // // vectorField masterZonePointTraction =
            // //     zoneToZone().slaveToMasterPointInterpolate
            // //     (
            // //         slaveZonePointTraction
            // //     );

            // PatchToPatchInterpolation<primitiveFacePatch, primitiveFacePatch>
            //     zoneToZoneInvDist(shadowZone(), zone());

            // vectorField masterZonePointTraction =
            //     zoneToZoneInvDist.pointInterpolate
            //     (
            //         slaveZonePointTraction
            //     );

            // vectorField masterZoneTraction =
            //     zonePointToFaceInterpolate
            //     (
            //         zoneIndex(), masterZonePointTraction
            //     );

            vectorField masterPatchTraction =
                patchField
                (
                    patch().index(),
                    zoneIndex(),
                    masterZoneTraction
                );

            refGrad() =
                tractionBoundaryGradient().snGrad
                (
                    masterPatchTraction, // surface traction
                    scalarField(patch().size(),0.0),  // surface pressure
                    fieldName_,                       // working field name
                    "U",                              // total field name
                    patch(),                          // polyPatch
                    bool(fieldName_ == "DU")          // incremental
                );
        }
    }

    // Remove components in 3rd direction for 2-D meshes
    if (patch().boundaryMesh().mesh().nSolutionD() == 2)
    {
        const Vector<label>& solD = patch().boundaryMesh().mesh().solutionD();

        forAll(solD, dirI)
        {
            if (solD[dirI] == -1)
            {
                forAll(patch(), facei)
                {
                    refValue()[facei][dirI] = 0.0;
                    refGrad()[facei][dirI] = 0.0;
                }
            }
        }
    }

    directionMixedFvPatchVectorField::updateCoeffs();
}


void Foam::solidContactFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField normalValue = transform(valueFraction(), refValue());

    //- non-orthogonal correction vectors needed to calculate gradValue
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField gradValue =
      this->patchInternalField()
      + (k&gradField.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
      transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}

Foam::tmp<Foam::Field<Foam::vector> > Foam::solidContactFvPatchVectorField::
snGrad() const
{
    vectorField pif = this->patchInternalField();

    vectorField normalValue = transform(valueFraction(), refValue());

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    //- correction vector
    vectorField k = delta - n*(n&delta);

    vectorField gradValue =
      this->patchInternalField()
      + (k&gradField.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    vectorField transformGradValue =
      transform(I - valueFraction(), gradValue);

    return
      (
       (normalValue + transformGradValue)
       - (pif + (k&gradField.patchInternalField()))
       )*this->patch().deltaCoeffs();
}


//- Increment of dissipated energy due to friction
Foam::tmp<Foam::scalarField> Foam::solidContactFvPatchVectorField::Qc() const
{
    // Consider storing Qc instead of recalculting multiple times

    if (!master())
    {
        FatalErrorIn
            (
                "solidContact::Qc()"
            )
            << "Only master can call Qc function!"
            << abort(FatalError);
    }

    tmp<scalarField> tQc(new scalarField(patch().size(), 0.0));
    scalarField& Qc = tQc();

    // For now, we assume traction is constant over time-step
    // Todo: use trapezoidal rule
    vectorField curTraction(Qc.size(), vector::zero);

    // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
    // has converged and sigma has been updated and mesh moved
    if
    (
        this->db().objectRegistry::foundObject<volSymmTensorField>
        (
            "sigmaCauchy"
        )
    )
    {
        const symmTensorField& sigma =
            this->db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigmaCauchy"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }
    else
    {
        const symmTensorField& sigma =
            this->db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }

    // Calculate slip

    vectorField slavePatchSlip = frictionModel().slip();

    vectorField slaveZoneSlip =
        zoneField
        (
            shadowZoneIndex(),
            shadowPatchIndex(),
            slavePatchSlip
        );

    // Interpolate from slave to master
    // Face-to-face

    vectorField masterZoneSlip =
        zoneToZone().slaveToMaster(slaveZoneSlip);

    // Face-to-point -- point-to-point -- point-to-face
    // Does not seem to solve problem of "jittering" stresses during
    // sliding
    // vectorField slaveZonePointSlip =
    //     zoneFaceToPointInterpolate
    //     (
    //         shadowZoneIndex(), slaveZoneSlip
    //     );

    // // vectorField masterZonePointSlip =
    // //     zoneToZone().slaveToMasterPointInterpolate
    // //     (
    // //         slaveZonePointSlip
    // //     );

    // PatchToPatchInterpolation<primitiveFacePatch, primitiveFacePatch>
    //     zoneToZoneInvDist(shadowZone(), zone());

    // vectorField masterZonePointSlip =
    //     zoneToZoneInvDist.pointInterpolate
    //     (
    //         slaveZonePointSlip
    //     );

    // vectorField masterZoneSlip =
    //     zonePointToFaceInterpolate
    //     (
    //         zoneIndex(), masterZonePointSlip
    //     );

    vectorField masterPatchSlip =
        patchField
        (
            patch().index(),
            zoneIndex(),
            masterZoneSlip
        );

    // Bug fix 12-Sep-14: heat flux rate, we must divide by time-step
    scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    // Increment of dissipated frictional energy for this timestep
    // The dot product of the traction vectors and the slip vectors gives the
    // dissipated frictional energy per unit area; which is always positive
    Qc = mag(curTraction & (masterPatchSlip/deltaT));

    return tQc;
}


// Write
void Foam::solidContactFvPatchVectorField::write(Ostream& os) const
{
    //Info<< "writing..."<<flush;
    directionMixedFvPatchVectorField::write(os);

    os.writeKeyword("master") << master_ << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch")
        << patch().boundaryMesh().mesh().boundary()[shadowPatchIndex()].name()
        << token::END_STATEMENT << nl;
    os.writeKeyword("orthotropic")
        << orthotropic_ << token::END_STATEMENT << nl;
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("contactActive") << contactActive_
        << token::END_STATEMENT << nl;

    if (master_)
    {
        os.writeKeyword("rigidMaster") << rigidMaster_
            << token::END_STATEMENT << nl;
        if (contactActive_)
        {
            os.writeKeyword("normalContactModel")
                << normalModel().type() << token::END_STATEMENT << nl;
            normalModel().writeDict(os);

            os.writeKeyword("frictionContactModel") << frictionModel().type()
                << token::END_STATEMENT << nl;
            frictionModel().writeDict(os);

            os.writeKeyword("correctionFrequency") << correctionFreq()
                << token::END_STATEMENT << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchVectorField, solidContactFvPatchVectorField);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace Foam

// ************************************************************************* //
