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

Description
    Mixing plane class dealing with transfer of data between two
    primitivePatches

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "MixingPlaneInterpolation.H"
#include "demandDrivenData.H"
#include "PrimitivePatch.H"
#include "IOmanip.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Compute master and slave patch addressing and weighting factors
//
//  Addressing:
//     Given the interpolation profile, find under which profile interval falls
//     any given  master and slave patch faces.
//
//  Weighting factors:
//     Find the amount of area intersection any given patch faces has with an
//     imaginary ribbon spanned from the profile interval that overlaps it.
//          Weighting factor == 1.0:  The face is completely overlapped
//           by its ribbon,
//          Weighting factor == 0.0:  The face is not overlapped by its ribbon.
//          GGI weighting factors are well suited for this.
//
//  For every interpolation profile intervals, we got:
//        'n' addresses and weighting factors for the master patch
//        'm' addresses and weighting factors for the slave patch
//

template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::calcAddressing() const
{
    if
    (
       masterPatchToProfileAddrPtr_
    || masterProfileToPatchAddrPtr_
    || masterPatchToProfileWeightsPtr_
    || masterProfileToPatchWeightsPtr_
    || slavePatchToProfileAddrPtr_
    || slaveProfileToPatchAddrPtr_
    || slavePatchToProfileWeightsPtr_
    || slaveProfileToPatchWeightsPtr_
    )
    {
        FatalErrorIn
        (
            "void MixingPlaneInterpolation<MasterPatch, "
            "SlavePatch>::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    Info<< "Creating internal GGIs: Large values for the master GGI "
        << "weighting factor corrections are expected."
        << endl;

    // Construct 2 GGIs in order to evaluate the interpolation weighting
    // factors and the addressing
    // - Since the profile cannot exactly match either master and slave
    // patches perfectly (discretization effects),
    // - it is best to use a GGI in order to compute those parameters

    // Basic GGI with no mesh transform nor separation
    const tensorField noTransform(0);
    const vectorField noTranslation(0);

    // NB: It is important here that the "ribbon patch" be the master patch
    // for each GGIs.
    //     We have to remember that when evaluating the GGI weighting factors,
    //     the slave patch faces points
    //     are projected onto a plane defined at every single master patch
    //     faces. We then compute a 2D polygon
    //     intersection using Sutherland-Hodgman, and so on...
    //     It is thus important that the "ribbon patch" be the same master
    //     patch for both GGI because the faces
    //     for both masterPatchCylCoord and slavePatchCylCoord will then be
    //     projected onto a common geometry.
    //     By doing so, we hope to minimize geometrical discretization errors
    //     introduced by the possibly different
    //     master and slave patch mesh resolution.

    GGIInterpolation<standAlonePatch, standAlonePatch>
        masterCircumAvgPatchToPatch
        (
            this->mixingPlanePatch(),
            this->transformedMasterPatch(),
            noTransform,
            noTransform,
            noTranslation,
            0.0,
            0.0,
            true    // Scale GGI weights
        );

    GGIInterpolation<standAlonePatch, standAlonePatch>
        slaveCircumAvgPatchToPatch
        (
            this->mixingPlanePatch(),
            this->transformedShadowPatch(),
            noTransform,
            noTransform,
            noTranslation,
            0.0,
            0.0,
            true    // Scale GGI weights
        );

    // Memorized the GGI addressing and weighting factors
    //
    // Basically, the master/slave weights and the master/slave addr values
    //  from both GGI gives us the information necessary to compute the
    // fields circumferential average on the master side, and to transfer
    // that averaged value properly back to the slave side.
    //
    // We get, for each master/slave patches: which face contribute to which
    // ribbon, and in which proportion.
    //
    // The GGI weighting factors will be use to compute the circumferential
    // weighted average.
    // Since the GGI weighting factors are in fact already factoring in the
    // ratio of intersected area versus the full facet area,
    // we do end up with a surface-area weighted average when using the GGI
    // weighting factors... Pretty neat... MB :)

    // Once we collect this information, we will can simply discard the GGIs
    // we no longer need them

    // The transfer from patch to profile requires the master
    //  weighting factors so we can do a proper facet surface weighted average
    masterPatchToProfileAddrPtr_ =
        new labelListList(masterCircumAvgPatchToPatch.masterAddr());

    masterPatchToProfileWeightsPtr_ =
        new scalarListList(masterCircumAvgPatchToPatch.masterWeights());

    slavePatchToProfileAddrPtr_ =
        new labelListList(slaveCircumAvgPatchToPatch.masterAddr());

    slavePatchToProfileWeightsPtr_ =
        new scalarListList(slaveCircumAvgPatchToPatch.masterWeights());


    // The transfer from profile to patch requires the slave weighting
    // factors so we can transfer the information back to the
    // slave patch using the standard GGI mechanism
    masterProfileToPatchAddrPtr_ =
        new labelListList(masterCircumAvgPatchToPatch.slaveAddr());

    masterProfileToPatchWeightsPtr_ =
        new scalarListList(masterCircumAvgPatchToPatch.slaveWeights());

    slaveProfileToPatchAddrPtr_ =
        new labelListList(slaveCircumAvgPatchToPatch.slaveAddr());

    slaveProfileToPatchWeightsPtr_ =
        new scalarListList(slaveCircumAvgPatchToPatch.slaveWeights());
}


//  Generate Cartesian/Cylindrical tranformation tensor fields for
//  master and slave patches
template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::calcTransforms() const
{
    // Collapse in spawise direction
    const direction spanDir = spanwiseSwitch();

    // Master side

    masterPatchToProfileTPtr_ = new tensorField(masterPatch_.size());
    tensorField& mPatchToProfileT = *masterPatchToProfileTPtr_;

    masterProfileToPatchTPtr_ = new tensorField(masterPatch_.size());
    tensorField& mProfileToPatchT = *masterProfileToPatchTPtr_;

    // Get master patch face normals
    const vectorField& globalMasterNormals = masterPatch_.faceNormals();

    // Move face normals into the local coordinate system
    vectorField localMasterNormals = cs_.localVector(globalMasterNormals);

    localMasterNormals.replace(spanDir, 0);

    // Transform back to global
    vectorField transformMasterNormals = cs_.globalVector(localMasterNormals);

    // Calculate transform tensors
    mPatchToProfileT =
        rotationTensor
        (
            globalMasterNormals,
            transformMasterNormals
        );

    mProfileToPatchT =
        rotationTensor
        (
            transformMasterNormals,
            globalMasterNormals
        );

    if (debug != 0)
    {
        // We can remove this later when agree this is ok.
        const vectorField& masterFaceCntr = masterPatch_.faceCentres();

        pointField transformedFaceCentre1 = transform
        (
            mPatchToProfileT,
            masterFaceCntr
        );

        pointField transformedFaceCentre2 = transform
        (
            mProfileToPatchT,
            transformedFaceCentre1
        );

        Info << "::calcTransforms(): globalMasterNormals: " << globalMasterNormals << endl;
        Info << "::calcTransforms(): localMasterNormals: " << localMasterNormals << endl;
        Info << "::calcTransforms(): transformMasterNormals: " << transformMasterNormals << endl;
        Info << "::calcTransforms(): masterFaceCntr: " << masterFaceCntr << endl;
        Info << "::calcTransforms(): transformedFaceCentre1: " << transformedFaceCentre1 << endl;
        Info << "::calcTransforms(): transformedFaceCentre2: " << transformedFaceCentre2 << endl;

        InfoIn
        (
            "MixingPlaneInterpolation"
            "<MasterPatch, SlavePatch>::calcTransforms()"
        )   << "master face centre transformation check: "
            << "(should be zero!) = "
            << " mag    : " << mag(transformedFaceCentre2 - masterFaceCntr) << endl
            << " sum mag: " << sum(mag(transformedFaceCentre2 - masterFaceCntr))
            << endl;
    }


    // Slave side

    slavePatchToProfileTPtr_ = new tensorField(slavePatch_.size());
    tensorField& sPatchToProfileT = *slavePatchToProfileTPtr_;

    slaveProfileToPatchTPtr_ = new tensorField(slavePatch_.size());
    tensorField& sProfileToPatchT = *slaveProfileToPatchTPtr_;

    // Get slave patch face normals
    const vectorField& globalSlaveNormals = slavePatch_.faceNormals();

    // Move face normals into the local coordinate system
    vectorField localSlaveNormals = cs_.localVector(globalSlaveNormals);

    localSlaveNormals.replace(spanDir, 0);

    // Transform back to global
    vectorField transformSlaveNormals = cs_.globalVector(localSlaveNormals);

    // Calculate transform tensors
    sPatchToProfileT =
        rotationTensor
        (
            globalSlaveNormals,
            transformSlaveNormals
        );

    sProfileToPatchT =
        rotationTensor
        (
            transformSlaveNormals,
            globalSlaveNormals
        );

    if (debug != 0)
    {
        // We can remove this later when agree this is ok.
        const vectorField& slaveFaceCntr = slavePatch_.faceCentres();

        pointField transformedFaceCentre1 = transform
        (
            sPatchToProfileT,
            slaveFaceCntr
        );

        pointField transformedFaceCentre2 = transform
        (
            sProfileToPatchT,
            transformedFaceCentre1
        );

        Info << "::calcTransforms(): globalSlaveNormals: " << globalSlaveNormals << endl;
        Info << "::calcTransforms(): localSlaveNormals: " << localSlaveNormals << endl;
        Info << "::calcTransforms(): transformSlaveNormals: " << transformSlaveNormals << endl;
        Info << "::calcTransforms(): slaveFaceCntr: " << slaveFaceCntr << endl;
        Info << "::calcTransforms(): transformedFaceCentre1: " << transformedFaceCentre1 << endl;
        Info << "::calcTransforms(): transformedFaceCentre2: " << transformedFaceCentre2 << endl;

        InfoIn
        (
            "MixingPlaneInterpolation"
            "<SlavePatch, SlavePatch>::calcTransforms()"
        )   << "slave face centre transformation check: "
            << "(should be zero!) = "
            << " mag    : " << mag(transformedFaceCentre2 - slaveFaceCntr) << endl
            << " sum mag: " << sum(mag(transformedFaceCentre2 - slaveFaceCntr))
            << endl;
    }
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearAddressing()
{
    // Master side
    deleteDemandDrivenData(masterPatchToProfileAddrPtr_);
    deleteDemandDrivenData(masterProfileToPatchAddrPtr_);
    deleteDemandDrivenData(masterPatchToProfileWeightsPtr_);
    deleteDemandDrivenData(masterProfileToPatchWeightsPtr_);

    // Slave side
    deleteDemandDrivenData(slavePatchToProfileAddrPtr_);
    deleteDemandDrivenData(slaveProfileToPatchAddrPtr_);
    deleteDemandDrivenData(slavePatchToProfileWeightsPtr_);
    deleteDemandDrivenData(slaveProfileToPatchWeightsPtr_);
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearTransforms()
{
    // Master side
    deleteDemandDrivenData(masterPatchToProfileTPtr_);
    deleteDemandDrivenData(masterProfileToPatchTPtr_);

    // Slave side
    deleteDemandDrivenData(slavePatchToProfileTPtr_);
    deleteDemandDrivenData(slaveProfileToPatchTPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
