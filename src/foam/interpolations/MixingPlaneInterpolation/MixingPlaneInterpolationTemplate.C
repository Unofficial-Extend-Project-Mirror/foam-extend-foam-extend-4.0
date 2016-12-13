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

Description
    Mixing plane class dealing with transfer of data between two
    primitivePatches

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "MixingPlaneInterpolationTemplate.H"
#include "demandDrivenData.H"
#include "PrimitivePatchTemplate.H"
#include "IOmanip.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
Foam::direction
MixingPlaneInterpolation<MasterPatch, SlavePatch>::sweepAxisSwitch() const
{
    // SweepAxis switch
    switch (sweepAxisType_)
    {
        case SWEEP_X:
        case SWEEP_R:
        {
            return vector::X;
        }
        break;

        case SWEEP_Y:
        case SWEEP_THETA:
        {
            return vector::Y;
        }
        break;

        case SWEEP_Z:
        {
            return vector::Z;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "direction MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::sweepAxisSwitch() const"
            )   << "Bad sweepAxis type: "
                << MixingPlaneInterpolationName::sweepAxisNames_
                       [sweepAxisType_]
                << "Available types: "
                << MixingPlaneInterpolationName::sweepAxisNames_
                << abort(FatalError);

            // Dummy return
            return vector::X;
        }
    }
}


template<class MasterPatch, class SlavePatch>
Foam::direction
MixingPlaneInterpolation<MasterPatch, SlavePatch>::stackAxisSwitch() const
{
    // stackAxis switch
    switch (stackAxisType_)
    {
        case STACK_X:
        case STACK_R:
        {
            return vector::X;
        }
        break;

        case STACK_Y:
        case STACK_THETA:
        {
            return vector::Y;
        }
        break;

        case STACK_Z:
        {
            return vector::Z;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "direction MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::stackAxisSwitch() const"
            )   << "Bad stackAxis type: "
                << MixingPlaneInterpolationName::stackAxisNames_
                       [stackAxisType_]
                << "Available types: "
                << MixingPlaneInterpolationName::stackAxisNames_
                << abort(FatalError);

            // Dummy return
            return vector::X;
        }
    }
}


template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearOut()
{
    clearTransfomedPatches();
    clearMixingPlanePatch();

    clearAddressing();
    clearTransforms();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
MixingPlaneInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch& slavePatch,
    const coordinateSystem& cs,
    const MixingPlaneInterpolationName::discretisation& discretisationType,
    const MixingPlaneInterpolationName::sweepAxis& sweepAxisType,
    const MixingPlaneInterpolationName::stackAxis& stackAxisType,
    const pointField& interpolationProfile
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    cs_(cs),
    discretisationType_(discretisationType),
    sweepAxisType_(sweepAxisType),
    stackAxisType_(stackAxisType),
    interpolationProfile_(interpolationProfile),

    forwardT_(),
    reverseT_(),
    forwardSep_(),

    transformedMasterPatchPtr_(NULL),
    transformedShadowPatchPtr_(NULL),
    mixingPlanePatchPtr_(NULL),

    masterPatchToProfileTPtr_(NULL),
    masterProfileToPatchTPtr_(NULL),
    slavePatchToProfileTPtr_(NULL),
    slaveProfileToPatchTPtr_(NULL),

    masterPatchToProfileAddrPtr_(NULL),
    masterProfileToPatchAddrPtr_(NULL),
    masterPatchToProfileWeightsPtr_(NULL),
    masterProfileToPatchWeightsPtr_(NULL),

    slavePatchToProfileAddrPtr_(NULL),
    slaveProfileToPatchAddrPtr_(NULL),
    slavePatchToProfileWeightsPtr_(NULL),
    slaveProfileToPatchWeightsPtr_(NULL)

{
    // Check size of transform.  They should be equal to slave patch size
    // if the transform is not constant
    if (forwardT_.size() > 1 || reverseT_.size() > 1)
    {
        if
        (
            forwardT_.size() != slavePatch_.size()
         || reverseT_.size() != masterPatch_.size()
        )
        {
            FatalErrorIn
            (
                "MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::MixingPlaneInterpolation"
            )   << "Incorrectly defined transform: forwardT: "
                << forwardT_.size() << " patch: " << slavePatch_.size()
                << " reverseT: " << reverseT_.size()
                << " patch: " << masterPatch_.size()
                << abort(FatalError);
        }
    }

    masterPatchToProfileAddr();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::~MixingPlaneInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::mixingPlanePatch() const
{
    if (!mixingPlanePatchPtr_)
    {
        calcMixingPlanePatch();
    }

    return *mixingPlanePatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
transformedMasterPatch() const
{
    if (!transformedMasterPatchPtr_)
    {
        calcTransformedPatches();
    }

    return *transformedMasterPatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
transformedShadowPatch() const
{
    if (!transformedShadowPatchPtr_)
    {
        calcTransformedPatches();
    }

    return *transformedShadowPatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileAddr() const
{
    if (!masterPatchToProfileAddrPtr_)
    {
        calcAddressing();
    }

    return *masterPatchToProfileAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchAddr() const
{
    if (!masterProfileToPatchAddrPtr_)
    {
        calcAddressing();
    }

    return *masterProfileToPatchAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileWeights() const
{
    if (!masterPatchToProfileWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterPatchToProfileWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchWeights() const
{
    if (!masterProfileToPatchWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterProfileToPatchWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slavePatchToProfileAddr() const
{
    if (!slavePatchToProfileAddrPtr_)
    {
        calcAddressing();
    }

    return *slavePatchToProfileAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slaveProfileToPatchAddr() const
{
    if (!slaveProfileToPatchAddrPtr_)
    {
        calcAddressing();
    }

    return *slaveProfileToPatchAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slavePatchToProfileWeights() const
{
    if (!slavePatchToProfileWeightsPtr_)
    {
        calcAddressing();
    }

    return *slavePatchToProfileWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slaveProfileToPatchWeights() const
{
    if (!slaveProfileToPatchWeightsPtr_)
    {
        calcAddressing();
    }

    return *slaveProfileToPatchWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileT() const
{
    if (!masterPatchToProfileTPtr_)
    {
        calcTransforms();
    }

    return *masterPatchToProfileTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchT() const
{
    if (!masterProfileToPatchTPtr_)
    {
        calcTransforms();
    }

    return *masterProfileToPatchTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slavePatchToProfileT() const
{
    if (!slavePatchToProfileTPtr_)
    {
        calcTransforms();
    }

    return *slavePatchToProfileTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveProfileToPatchT() const
{
    if (!slaveProfileToPatchTPtr_)
    {
        calcTransforms();
    }

    return *slaveProfileToPatchTPtr_;
}


template<class MasterPatch, class SlavePatch>
bool MixingPlaneInterpolation<MasterPatch, SlavePatch>::movePoints()
{
    clearOut();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
