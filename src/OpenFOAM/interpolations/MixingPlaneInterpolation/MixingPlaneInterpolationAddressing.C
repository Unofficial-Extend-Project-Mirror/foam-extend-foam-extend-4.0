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
#include "RodriguesRotation.H"

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
    // Collapse in spanwise direction
    const direction spanDir = spanwiseSwitch();

#define useRodrigues
#ifndef useRodrigues
    const direction directionDir = directionalSwitch();
#endif

    // Master side
    masterPatchToProfileTPtr_ = new tensorField(masterPatch_.size());
    tensorField& mPatchToProfileT = *masterPatchToProfileTPtr_;

    masterProfileToPatchTPtr_ = new tensorField(masterPatch_.size());
    tensorField& mProfileToPatchT = *masterProfileToPatchTPtr_;

    if (cs_.name() == "coordinateSystem")
    {
        // Identity tensor for cartesian space
        mPatchToProfileT = tensor(sphericalTensor::I);
        mProfileToPatchT = tensor(sphericalTensor::I);
    }
    else
    {
        // Get master patch face centers.
        // Warning: Face normals are not a good choice here. We need a vector
        // that will intersect the rotation axis. Furthermore, if the surface
        // normals are all parallel with the rotation axis, no valid rotation
        // tensors can be computed
        // Use face cell centers cause the fields are taken from cell centers
        vectorField globalMasterVectors =
            masterPatch_.faceCellCentres() - cs_.origin(); 

        // We need unit vectors for computing rotation tensor
        // We also need vectors lying into the plane normal to the rotation
        // axis. This is a major limitation of the current implementation of
        // rotationTensor() that is being used for computing rotation tensors:
        // we cannot specify the rotation axis. Instead, when using
        // rotationTensor(), the rotation axis will be aligned with the normal
        // of the plane spanned by the two vectors.
        // RodriguesRotation() was designed to overcome these limitations.
#ifndef useRodrigues
       globalMasterVectors.replace(directionDir, 0);

        Info << "Debug: " << masterPatch_.faceCellCentres()  - globalMasterVectors << endl;
        globalMasterVectors /= mag(globalMasterVectors);
#endif

        // Move face vector into the local coordinate system
        vectorField localMasterVectors = cs_.localVector(globalMasterVectors);

        // Translate everything to theta=0
        localMasterVectors.replace(spanDir, 0);

        // Transform back to global
        vectorField transformMasterVectors = cs_.globalVector(localMasterVectors);

        if( debug )
        {
            Info << "Debug: globalMasterVectors: " << globalMasterVectors << endl;
            Info << "Debug: transformMasterVectors: " << transformMasterVectors << endl;
        }

        // Calculate transform tensors. These are pure rotation tensors, aligned
        // with the reference frame rotation axis
        mPatchToProfileT =
#ifndef useRodrigues
            rotationTensor
            (
#else
            RodriguesRotation
            (
                cs_.axis(),
#endif
                globalMasterVectors,
                transformMasterVectors
            );

        mProfileToPatchT =
#ifndef useRodrigues
            rotationTensor
            (
#else
            RodriguesRotation
            (
                cs_.axis(),
#endif
                transformMasterVectors,
                globalMasterVectors
            );

#if 0
        {
            vectorField globalMasterVectors2 = masterPatch_.faceCellCentres() - cs_.origin(); 
            globalMasterVectors2 /= mag(globalMasterVectors);
            
            vectorField localMasterVectors2 = cs_.localVector(globalMasterVectors2);
            localMasterVectors2.replace(spanDir, 0);
            vectorField transformMasterVectors2 = cs_.globalVector(localMasterVectors2);

            Info << "master delta face transf: " << mag(mPatchToProfileT - RodriguesRotation(cs_.axis(), globalMasterVectors, transformMasterVectors)) << endl;  
            Info << "master delta face transf: " << mag(mProfileToPatchT - RodriguesRotation(cs_.axis(), transformMasterVectors, globalMasterVectors )) << endl;  

            Info << "master delta face: " << mag(mPatchToProfileT - RodriguesRotation(cs_.axis(), globalMasterVectors2, transformMasterVectors2)) << endl;  
            Info << "master delta face: " << mag(mProfileToPatchT - RodriguesRotation(cs_.axis(), transformMasterVectors2, globalMasterVectors2 )) << endl;  

        }
#endif

        if (debug != 0)
        {
            // We can remove this later when agree this is ok.
            const vectorField& masterFaceCntr = masterPatch_.faceCellCentres();

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

            Info << "::calcTransforms(): mPatchToProfileT: " << mPatchToProfileT << endl;
            Info << "::calcTransforms(): mProfileToPatchT: " << mProfileToPatchT << endl;

            Info << "::calcTransforms(): globalMasterVectors: " << globalMasterVectors << endl;
            Info << "::calcTransforms(): localMasterVectors: " << localMasterVectors << endl;
            Info << "::calcTransforms(): transformMasterVectors: " << transformMasterVectors << endl;
            Info << "::calcTransforms(): masterFaceCntr: " << masterFaceCntr << endl;
            Info << "::calcTransforms(): transformedFaceCentre1: " << transformedFaceCentre1 << endl;
            Info << "::calcTransforms(): transformedFaceCentre2: " << transformedFaceCentre2 << endl;
            Info << "::calcTransforms(): deltaMag direct: " << mag(transformedFaceCentre1) - mag(masterFaceCntr) << endl;
            Info << "::calcTransforms(): deltaMag reverse: " << mag(transformedFaceCentre2) - mag(masterFaceCntr) << endl;
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

    }

    // Slave side

    slavePatchToProfileTPtr_ = new tensorField(slavePatch_.size());
    tensorField& sPatchToProfileT = *slavePatchToProfileTPtr_;

    slaveProfileToPatchTPtr_ = new tensorField(slavePatch_.size());
    tensorField& sProfileToPatchT = *slaveProfileToPatchTPtr_;

    if (cs_.name() == "coordinateSystem")
    {
        // Identity tensor for cartesian space
        sPatchToProfileT = tensor(sphericalTensor::I);
        sProfileToPatchT = tensor(sphericalTensor::I);
    }
    else
    {
        // Get slave patch face centers.
        // Warning: Face normals are not a good choice here. We need a vector
        // that will intersect the rotation axis. Furthermore, if the surface
        // normals are all parallel with the rotation axis, no valid rotation
        // tensors can be computed 
        vectorField globalSlaveVectors =
            slavePatch_.faceCellCentres() - cs_.origin(); 

        // We need unit vectors for computing rotation tensor
        // We also need vectors lying into the plane normal to the rotation
        // axis. This is a major limitation of the current implementation of
        // rotationTensor() that is being used for computing rotation tensors:
        // we cannot specify the rotation axis. Instead, when using
        // rotationTensor(), the rotation axis will be aligned with the normal
        // of the plane spanned by the two vectors.
        // RodriguesRotation() was designed to overcome these limitations.
#ifndef useRodrigues
        globalSlaveVectors.replace(directionDir, 0);

        Info << "Debug: " << slavePatch_.faceCellCentres()  - globalSlaveVectors << endl;

        globalSlaveVectors /= mag(globalSlaveVectors);
#endif

        // Move face vector into the local coordinate system
        vectorField localSlaveVectors = cs_.localVector(globalSlaveVectors);

        // Translate everything to theta=0
        localSlaveVectors.replace( spanDir, 0 );

        // Transform back to global
        vectorField transformSlaveVectors = cs_.globalVector(localSlaveVectors);

        // Calculate transform tensors. These are pure rotation tensors, aligned
        // with the reference frame rotation axis
        sPatchToProfileT =
#ifndef useRodrigues
            rotationTensor
            (
#else
            RodriguesRotation
            (
                cs_.axis(),
#endif
                globalSlaveVectors,
                transformSlaveVectors
            );

        sProfileToPatchT =
#ifndef useRodrigues
            rotationTensor
            (
#else
            RodriguesRotation
            (
                cs_.axis(),
#endif
                transformSlaveVectors,
                globalSlaveVectors
            );

#if 0
        {
            vectorField globalSlaveVectors2 = slavePatch_.faceCellCentres() - cs_.origin(); 
            globalSlaveVectors2 /= mag(globalSlaveVectors);
            
            vectorField localSlaveVectors2 = cs_.localVector(globalSlaveVectors2);
            localSlaveVectors2.replace(spanDir, 0);
            vectorField transformSlaveVectors2 = cs_.globalVector(localSlaveVectors2);

            Info << "slave delta face transf: " << mag(sPatchToProfileT - RodriguesRotation(cs_.axis(), globalSlaveVectors, transformSlaveVectors)) << endl;  
            Info << "slave delta face transf: " << mag(sProfileToPatchT - RodriguesRotation(cs_.axis(), transformSlaveVectors, globalSlaveVectors )) << endl;  

            Info << "slave delta face: " << mag(sPatchToProfileT - RodriguesRotation(cs_.axis(), globalSlaveVectors2, transformSlaveVectors2)) << endl;  
            Info << "slave delta face: " << mag(sProfileToPatchT - RodriguesRotation(cs_.axis(), transformSlaveVectors2, globalSlaveVectors2 )) << endl;  

        }
#endif
        if (debug != 0)
        {
            // We can remove this later when agree this is ok.
            const vectorField& slaveFaceCntr = slavePatch_.faceCellCentres();

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

            Info << "::calcTransforms(): globalSlaveVectors: " << globalSlaveVectors << endl;
            Info << "::calcTransforms(): localSlaveVectors: " << localSlaveVectors << endl;
            Info << "::calcTransforms(): transformSlaveVectors: " << transformSlaveVectors << endl;
            Info << "::calcTransforms(): slaveFaceCntr: " << slaveFaceCntr << endl;
            Info << "::calcTransforms(): transformedFaceCentre1: " << transformedFaceCentre1 << endl;
            Info << "::calcTransforms(): transformedFaceCentre2: " << transformedFaceCentre2 << endl;
            Info << "::calcTransforms(): deltaMag direct: " << mag(transformedFaceCentre1) - mag(slaveFaceCntr) << endl;
            Info << "::calcTransforms(): deltaMag reverse: " << mag(transformedFaceCentre2) - mag(slaveFaceCntr) << endl;

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
