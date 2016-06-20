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
    normalContactModel

\*---------------------------------------------------------------------------*/

#include "normalContactModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(normalContactModel, 0);
defineRunTimeSelectionTable(normalContactModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalContactModel::normalContactModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID,
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch
)
:
    name_(name),
    patch_(patch),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID),
    masterFaceZoneID_(masterFaceZoneID),
    slaveFaceZoneID_(slaveFaceZoneID),
    slaveContactPointGap_
    (
        patch.boundaryMesh().mesh().boundaryMesh()[slavePatchID].nPoints(), 0.0
    )
{}


normalContactModel::normalContactModel(const normalContactModel& nm)
:
    name_(nm.name_),
    patch_(nm.patch_),
    masterPatchID_(nm.masterPatchID_),
    slavePatchID_(nm.slavePatchID_),
    masterFaceZoneID_(nm.masterFaceZoneID_),
    slaveFaceZoneID_(nm.slaveFaceZoneID_),
    slaveContactPointGap_(nm.slaveContactPointGap_)
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //


void normalContactModel::calculateSlaveFaceNormals
(
    vectorField& slaveFaceNormals,
    const nonLinearGeometry::nonLinearType nonLinear,
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch,
    const ggiZoneInterpolation& zoneToZone
)
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label slavePatchIndex = slavePatchID();

    // dirichletNeumannCalculateSlaveFaceNormals.H
    if (nonLinear == nonLinearGeometry::OFF)
    {
        // undeformed normals
        // remember that the mesh has not moved, only the global face zone
        // patches are moved to the deformed position

        // slaveFaceNormals =
        //-1*masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
        //        (
        //         mesh.boundaryMesh()[masterPatchIndex].faceNormals()
        //         );
        // slaveFaceNormals /= mag(slaveFaceNormals);

        // OPTION 2
        // Use deformed slave face normals
        // Get them from global face zone patches
        // const vectorField& globalSlaveFaceNormals =
        //slaveFaceZonePatch.faceNormals();
        // const label slavePatchStart
        //    = mesh.boundaryMesh()[slavePatchIndex].start();
        // forAll(slaveFaceNormals, facei)
        //    {
        //      slaveFaceNormals[facei] =
        //        globalSlaveFaceNormals
        //        [
        //    mesh.faceZones()[slaveFaceZoneID()].whichFace
        // (slavePatchStart + facei)
        //         ];
        //    }
        // // make sure they are unity
        // slaveFaceNormals /= mag(slaveFaceNormals);

        // OPTION 3
        // globally interpolate master normals to the slave
        // then get the local normals

        // undeformed master normals
        // const vectorField& actualSlaveFaceNormals =
        //mesh.boundaryMesh()[slavePatchIndex].faceNormals();
        // vectorField globalMasterFaceNormals
        //(masterFaceZonePatch.size(), vector::zero);
        // const label masterPatchStart
        //    = mesh.boundaryMesh()[masterPatchIndex].start();
        // // master undeformed face normals
        // vectorField masterFaceNormals =
        //mesh.boundaryMesh()[masterPatchIndex].faceNormals();
        // // put field into global
        // forAll(masterFaceNormals, i)
        //    {
        //      globalMasterFaceNormals[mesh.faceZones()
        //[masterFaceZoneID()].whichFace(masterPatchStart + i)] =
        //        masterFaceNormals[i];
        //    }

        // //- exchange parallel data
        // // sum because each face is only on one proc
        // reduce(globalMasterFaceNormals, sumOp<vectorField>());

        // use deformed master normals
        const vectorField& globalMasterFaceNormals =
            masterFaceZonePatch.faceNormals();

        // interpolate master normals to the slave and reverse the direction
        // we can use inverseDistance or GGI for the interpolation
        // GGI is much better is we want the actual master normals
        vectorField globalSlaveFaceNormals =
            -zoneToZone.masterToSlave(globalMasterFaceNormals);

        vectorField actualSlaveFaceNormals = slaveFaceZonePatch.faceNormals();

        // put global back into local
        const label slavePatchStart
            = mesh.boundaryMesh()[slavePatchIndex].start();
        forAll(slaveFaceNormals, facei)
        {
            slaveFaceNormals[facei] =
                globalSlaveFaceNormals
                [
                    mesh.faceZones()[slaveFaceZoneID()].whichFace
                    (
                        slavePatchStart + facei
                    )
                ];

            if (mag(slaveFaceNormals[facei]) < SMALL)
            {
                // interpolation sometimes does not work for far
                // away faces so we will use actual slave normals
                // but doesn't really matter as the face is probably
                // not in contact
                slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
            }
            else
            {
                // make sure they are unity
                slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
            }
        }
    }
    else if
    (
        nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN
        || nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN
    )
    {
        Info << "nonLinear UL/TL normals" << endl;

        // a few different ways to calculate deformed normals

        // OPTION 1
        // use deformed slave normals calculated using F
        // calculate deformed normals using deformation gradient
        // const volTensorField& gradField =
        // mesh.objectRegistry::lookupObject<volTensorField>
        //("grad("+fieldName+")");
        //        // deformation gradient
        //        tensorField F = gradField.boundaryField()[slavePatchIndex];
        //        if (fieldName == "DU")
        //    {
        //      const volTensorField& gradU =
        //        mesh.objectRegistry::lookupObject<volTensorField>("grad(U)");
        //      F += gradU.boundaryField()[slavePatchIndex];
        //    }
        //        const tensorField Finv = inv(I + F);
        // deformed normals
        //slaveFaceNormals =
        // J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();


        // OPTION 2
        // Use deformed slave face normals
        // Get them from global face zone patches
        const vectorField& globalSlaveFaceNormals =
            slaveFaceZonePatch.faceNormals();
        const label slavePatchStart
            = mesh.boundaryMesh()[slavePatchIndex].start();
        forAll(slaveFaceNormals, facei)
        {
            slaveFaceNormals[facei] =
                globalSlaveFaceNormals
                [
                    mesh.faceZones()[slaveFaceZoneID()].whichFace
                    (
                        slavePatchStart + facei
                    )
                ];
        }
        // make sure they are unity
        slaveFaceNormals /= mag(slaveFaceNormals);


        // OPTION 3
        // Use deformed master normals
        // globally interpolate deformed master normals to the slave
        // then get the local normals
        // globally interpolate master normals to the slave
        // then get the local normals

        // undeformed master normals
        //const vectorField& actualSlaveFaceNormals =
        //mesh.boundaryMesh()[slavePatchIndex].faceNormals();
        // vectorField globalMasterFaceNormals
        //(masterFaceZonePatch.size(), vector::zero);
        // const label masterPatchStart
        //    = mesh.boundaryMesh()[masterPatchIndex].start();
        // // master undeformed face normals
        // vectorField masterFaceNormals =
        //mesh.boundaryMesh()[masterPatchIndex].faceNormals();
        // // put field into global
        // forAll(masterFaceNormals, i)
        //    {
        //      globalMasterFaceNormals[mesh.faceZones()
        //[masterFaceZoneID()].whichFace(masterPatchStart + i)] =
        //        masterFaceNormals[i];
        //    }

        // //- exchange parallel data
        // // sum because each face is only on one proc
        // reduce(globalMasterFaceNormals, sumOp<vectorField>());

        // use deformed master normals
        // {
        // const vectorField& globalMasterFaceNormals =
        //masterFaceZonePatch.faceNormals();

        // // interpolate master normals to the slave and reverse the direction
        // // we can use inverseDistance or GGI for the interpolation
        // // GGI is much better is we want the actual master normals
        // vectorField globalSlaveFaceNormals
        //(slaveFaceZonePatch.size(), vector::zero);
        // if (ggiInterpolatorPtr)
        //    {
        //      //Info << "interpolating master normals with GGI" << endl;
        //      globalSlaveFaceNormals =
        //        -ggiInterpolatorPtr->masterToSlave(globalMasterFaceNormals);
        //    }
        // else // inverse distance
        //    {
        //      globalSlaveFaceNormals =
        //        -masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
        //        (
        //         globalMasterFaceNormals
        //         );
        //    }
        // vectorField actualSlaveFaceNormals =
        //    slaveFaceZonePatch.faceNormals();

        // // get local normals from global
        // const label slavePatchStart
        //    = mesh.boundaryMesh()[slavePatchIndex].start();
        // forAll(slaveFaceNormals, facei)
        //    {
        //      slaveFaceNormals[facei] =
        //        globalSlaveFaceNormals
        //        [
        //         mesh.faceZones()[slaveFaceZoneID()].whichFace
        //(slavePatchStart + facei)
        //         ];

        //      if (mag(slaveFaceNormals[facei]) < SMALL)
        //        {
        //          // interpolation sometimes does not work for far
        //          // away faces so we will use actual slave normals
        //          // but doesn't really matter as the face is probably
        //          // not in contact
        //          slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
        //        }
        //      else
        //        {
        //          // make sure they are unity
        //          slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
        //        }
        //    }
        // }
    }
    else if (nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF)
    {
        // OPTION 1
        // use deformed slave normals calculated using F
        // calculate deformed normals using deformation gradient

        // Lookup relative deformation gradient
        // const tensorField& Finv =
        //     mesh.objectRegistry::lookupObject<volTensorField>
        //     (
        //         "relFinv"
        //     ).boundaryField()[slavePatchIndex];
        // const scalarField& J =
        //     mesh.objectRegistry::lookupObject<volScalarField>
        //     (
        //         "J"
        //     ).boundaryField()[slavePatchIndex];

        // const label masterPatchIndex = masterPatchID();

        // Lookup relFinv
        const tensorField& Finv =
            mesh.objectRegistry::lookupObject<volTensorField>
            (
                "relFinv"
            ).boundaryField()[slavePatchIndex];
        // Lookup relFinv surface field
        // const tensorField& Finv =
        //     mesh.objectRegistry::lookupObject<surfaceTensorField>
        //     (
        //         "relFinvf"
        //     ).boundaryField()[slavePatchIndex];
        //).boundaryField()[masterPatchIndex];
        // const scalarField& J =
        //     mesh.objectRegistry::lookupObject<surfaceScalarField>
        //     (
        //         "Jf"
        //     ).boundaryField()[slavePatchIndex];

        // Calculate deformed normals by Nanson's formula

        // if (gMax(mag(slaveFaceNormals_)) < SMALL)
        // {
        //     slaveFaceNormals_ =
        //   //J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
        //       Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
        //     slaveFaceNormals_ /= mag(slaveFaceNormals_);
        // }
        // const vectorField slaveFaceNormalsPrevIter = slaveFaceNormals_;

        slaveFaceNormals =
            //J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
            Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();

        slaveFaceNormals /= mag(slaveFaceNormals);
    }
    else if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN_KIRCHHOFF)
    {
        // OPTION 1
        // use deformed slave normals calculated using F
        // calculate deformed normals using deformation gradient

        // Lookup relative deformation gradient and inverse from the solver
        const tensorField& Finv =
            mesh.objectRegistry::lookupObject<volTensorField>
            (
                "Finv"
            ).boundaryField()[slavePatchIndex];

        const scalarField& J =
            mesh.objectRegistry::lookupObject<volScalarField>
            (
                "J"
            ).boundaryField()[slavePatchIndex];


        // Calculate deformed normals by Nanson's formula
        slaveFaceNormals =
            J*Finv.T() & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
    }
    else if
    (
        nonLinear == nonLinearGeometry::DEFORMED_LAGRANGIAN
    )
    {
        //Info << "nonLinear DL normals" << endl;

        // The main mesh is always at the deformed position for
        // deformedLagrangian
        // So we will lookup up the slaveFaceZone from the main mesh as this was
        // moved using
        // least squares interpolation which is more accurate than
        //primitvePatchIntepolation

        // slave normals
        // Philipc: fix for parallel - use local patch not global face zone
        //slaveFaceNormals =
        //mesh.faceZones()[slaveFaceZoneID()]().faceNormals();
        slaveFaceNormals = mesh.boundaryMesh()[slavePatchID()].faceNormals();

        // or interpolate master normals
        //...
    }
    else
    {
        FatalError
            << "dirichletNeumann::correct()" << nl
                << "dirichletNeumannCalculateSlaveFaceNormals" << nl
                << "nonLinear option " << nonLinear
                << " is unknown"
                << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
