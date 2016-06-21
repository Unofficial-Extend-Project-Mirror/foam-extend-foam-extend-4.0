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
    MixingPlane interpolation functions

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::toProfile
(
    const Field<Type>& srcF,
    const labelListList& srcAddr,
    const scalarListList& srcWeights,
    Field<Type>& profileBandValues
)
{
    // The src to profile transfer is done using weighted averaging
    // evaluation of srcF.
    // See:  http://en.wikipedia.org/wiki/Weighted_mean
    //
    //                            sum(w*phi)
    //         average(phi)  ==  -----------
    //                             sum(w)
    //
    // Note however that in our case, the weighting factors utilized
    // for the computation of the weighted average are GGI weighting
    // factors. So for every interpolation band, we have the following
    // situation, thanks to the intrinsic conservativeness of the GGI
    // interface:
    //
    //         sum(w) == 1.0
    //
    //  which leads to the following simplification:
    //
    //         average(phi)  ==  sum(w*phi)

    forAll (srcAddr, bandI)
    {
        forAll (srcAddr[bandI], faceI)
        {
            profileBandValues[bandI] +=
                srcF[srcAddr[bandI][faceI]]*srcWeights[bandI][faceI];
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::fromProfile
(
    const Field<Type>& profileBandValues,
    const labelListList& dstAddr,
    const scalarListList& dstWeights,
    Field<Type>& dstResultF
)
{
    // The profile to dst transfer is done by distributing the
    // average value accordingly to the dst weighting factors
    forAll (dstAddr, faceI)
    {
        const labelList& curAddr = dstAddr[faceI];
        const scalarList& curW = dstWeights[faceI];

        dstResultF[faceI] = pTraits<Type>::zero;

        forAll (curAddr, bandI)
        {
            dstResultF[faceI] += profileBandValues[curAddr[bandI]]*curW[bandI];
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::maskedFromProfile
(
    const Field<Type>& profileBandValues,
    const labelListList& dstAddr,
    const scalarListList& dstWeights,
    Field<Type>& dstResultF,
    const labelList& mask
)
{
    // The profile to dst transfer is done by distributing the
    // average value accordingly to the dst weighting factors
    forAll (mask, maskI)
    {
        // Pick the masked face
        const label faceI = mask[maskI];

        const labelList& curAddr = dstAddr[faceI];
        const scalarList& curW = dstWeights[faceI];

        dstResultF[maskI] = pTraits<Type>::zero;

        forAll (curAddr, bandI)
        {
            dstResultF[maskI] += profileBandValues[curAddr[bandI]]*curW[bandI];
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::maskedTransform
(
    Field<Type>& transField,
    const tensorField& t,
    const Field<Type>& inField,
    const labelList& mask
)
{
    // The profile to dst transfer is done by distributing the
    // average value accordingly to the dst weighting factors
    forAll (mask, maskI)
    {
        // Pick the masked face
        const label faceI = mask[maskI];

        transField[maskI] = transform(t[faceI], inField[maskI]);
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::interpolate
(
    const Field<Type>& srcF,
    const labelListList& srcAddr,
    const scalarListList& srcWeights,
    const labelListList& dstAddr,
    const scalarListList& dstWeights,
    Field<Type>& dstResultF
) const
{
    Field<Type> profileBandValues(nProfileBands(), pTraits<Type>::zero);

    // Interpolate from patch to profile
    toProfile
    (
        srcF,
        srcAddr,
        srcWeights,
        profileBandValues
    );

    // profileBandValues are now the circumferentially averaged values

    // Collect from profile to patch
    fromProfile
    (
        profileBandValues,
        dstAddr,
        dstWeights,
        dstResultF
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToSlave
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::masterToSlave("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move master data from 'patch space' to 'profile space'
    // MB: We need this back
    Field<Type> profileFF = transform(masterPatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            slavePatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    interpolate
    (
        profileFF,                     // Master data in 'profile space'
        masterPatchToProfileAddr(),    // From master: compute the average
        masterPatchToProfileWeights(),
        slaveProfileToPatchAddr(),     // To slave distribute average from
        slaveProfileToPatchWeights(),  // profile to patch
        result
    );

    // Apply transform to bring the slave field back from 'profile space'
    // to 'patch space'
    transform(result, slaveProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToSlave
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = masterToSlave(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToMaster
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::slaveToMaster("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move slave data from 'patch space' to 'profile space'
     Field<Type> profileFF = transform(slavePatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            masterPatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    interpolate
    (
        profileFF,                     // Slave data in 'profile space'
        slavePatchToProfileAddr(),     // From slave: compute the average
        slavePatchToProfileWeights(),
        masterProfileToPatchAddr(),    // To master: distribute the average
        masterProfileToPatchWeights(),
        result
    );

    // Apply transform to bring the master field back from 'profile space'
    // to 'patch space'
    transform(result, masterProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToMaster
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = slaveToMaster(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToProfile
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::masterToProfile("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move master data from 'patch space' to 'profile space'
    Field<Type> profileFF = transform(masterPatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            nProfileBands(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    toProfile
    (
        profileFF,                      // Master data in 'profile space'
        masterPatchToProfileAddr(),     // From master: compute the average
        masterPatchToProfileWeights(),
        result
    );

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToProfile
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = masterToProfile(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToProfile
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::slaveToProfile("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move slave data from 'patch space' to 'profile space'
    Field<Type> profileFF = transform(slavePatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            nProfileBands(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    toProfile
    (
        profileFF,                     // Slave data in 'profile space'
        slavePatchToProfileAddr(),     // From slave: compute the average
        slavePatchToProfileWeights(),
        result
    );

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToProfile
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = slaveToProfile(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::profileToMaster
(
    const Field<Type>& profileFF
) const
{
    if (profileFF.size() != nProfileBands())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::profileToMaster("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to profile.  Profile size: "
            << nProfileBands() << " field size: " << profileFF.size()
            << abort(FatalError);
    }

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            masterPatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    fromProfile
    (
        profileFF,                      // Master data in 'profile space'
        masterProfileToPatchAddr(),     // To master: distribute the average
        masterProfileToPatchWeights(),
        result
    );

    // Apply transform to bring the master field back from 'profile space'
    // to 'patch space'
    transform(result, masterProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::profileToMaster
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = profileToMaster(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::maskedProfileToMaster
(
    const Field<Type>& profileFF,
    Field<Type>& result,
    const labelList& mask
) const
{
    if (profileFF.size() != nProfileBands() || result.size() != mask.size())
    {
        FatalErrorIn
        (
            "bvoid MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "maskedProfileToMaster\n"
            "(\n"
            "    const Field<Type>& profileFF,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to profile.  Profile size: "
            << nProfileBands() << " field size: " << profileFF.size()
            << " result size: " << result.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    // Do interpolation
    maskedFromProfile
    (
        profileFF,                      // Master data in 'profile space'
        masterProfileToPatchAddr(),     // To master: distribute the average
        masterProfileToPatchWeights(),
        result,
        mask
    );

    // Apply transform to bring the master field back from 'profile space'
    // to 'patch space'
    maskedTransform(result, masterProfileToPatchT(), result, mask);
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::profileToSlave
(
    const Field<Type>& profileFF
) const
{
    if (profileFF.size() != nProfileBands())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::profileToSlave("
            "const Field<Type> ff) const"
        )   << "given field does not correspond to profile.  Profile size: "
            << nProfileBands() << " field size: " << profileFF.size()
            << abort(FatalError);
    }

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            slavePatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    fromProfile
    (
        profileFF,                     // Slave data in 'profile space'
        slaveProfileToPatchAddr(),     // To slave distribute average from
        slaveProfileToPatchWeights(),  // profile to patch
        result
    );

    // Apply transform to bring the slave field back from 'profile space'
    // to 'patch space'
    transform(result, slaveProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::profileToSlave
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = profileToSlave(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::maskedProfileToSlave
(
    const Field<Type>& profileFF,
    Field<Type>& result,
    const labelList& mask
) const
{
    if (profileFF.size() != nProfileBands() || result.size() != mask.size())
    {
        FatalErrorIn
        (
            "void MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "maskedProfileToSlave\n"
            "(\n"
            "    const Field<Type>& profileFF,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to profile.  Profile size: "
            << nProfileBands() << " field size: " << profileFF.size()
            << " result size: " << result.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedFromProfile
    (
        profileFF,                     // Slave data in 'profile space'
        slaveProfileToPatchAddr(),     // To slave distribute average from
        slaveProfileToPatchWeights(),  // profile to patch
        result,
        mask
    );

    // Apply transform to bring the slave field back from 'profile space'
    // to 'patch space'
    maskedTransform(result, slaveProfileToPatchT(), result, mask);
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToMaster
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::masterToMaster("
            "const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move master data from 'patch space' to 'profile space'
    Field<Type> profileFF = transform(masterPatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            masterPatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    interpolate
    (
        profileFF,                      // Master data in 'profile space'
        masterPatchToProfileAddr(),     // From master: compute the average
        masterPatchToProfileWeights(),
        masterProfileToPatchAddr(),     // To master: distribute the average
        masterProfileToPatchWeights(),
        result
    );

    // Apply transform to bring the master field back from 'profile space'
    // to 'patch space'
    transform(result, masterProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToMaster
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = masterToMaster(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToSlave
(
    const Field<Type>& patchFF
) const
{
    if (patchFF.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "MixingPlaneInterpolation::slaveToSlave("
            "const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move slave data from 'patch space' to 'profile space'
    Field<Type> profileFF = transform(slavePatchToProfileT(), patchFF);

    // Do interpolation
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            slavePatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    interpolate
    (
        profileFF,                     // Slave data in 'profile space'
        slavePatchToProfileAddr(),     // From slave: compute the average
        slavePatchToProfileWeights(),
        slaveProfileToPatchAddr(),     // To slave: distribute the average
        slaveProfileToPatchWeights(),
        result
    );

    // Apply transform to bring the slave field back from 'profile space'
    // to 'patch space'
    transform(result, slaveProfileToPatchT(), result);

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToSlave
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = slaveToSlave(tff());
    tff.clear();
    return tint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
