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
    MixingPlane interpolation functions

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include <sys/time.h>
#include <sstream>
#include <iomanip>

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
template<class Type>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
interpolate
(
    const Field<Type>& srcF,
    const labelListList& srcAddr,
    const scalarListList& srcWeights,
    const labelListList& dstAddr,
    const scalarListList& dstWeights,
    Field<Type>& dstResultF
) const
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


    int nbrProfileBands = interpolationProfile_.size() - 1;

    List<Type> profileBandValues(nbrProfileBands, pTraits<Type>::zero);
    scalarField srcScalingValues(nbrProfileBands, 0.0);

    forAll (srcAddr, bandI)
    {
        forAll (srcAddr[bandI], faceI)
        {
            profileBandValues[bandI] +=
                srcF[srcAddr[bandI][faceI]]*srcWeights[bandI][faceI];

            // NB: The next operation should be computed only
            // once... and should sum up to 1.0 Let's keep this
            // operation for now, until the mixingPlane interface is
            // fully debugged (MB, 07/2010)
            srcScalingValues[bandI] += srcWeights[bandI][faceI];

            if (debug <= -200)
            {
                Info << "bande: " << bandI
                    << "  src valeur: " << srcF[srcAddr[bandI][faceI]]
                    << "  src weight: " << srcWeights[bandI][faceI] << endl;
            }
        }
    }

    // We don't need to divide the profileBandValues by the
    // srcScalingValues because the srcScalingValues are identically
    // equal to 1.0, thanks to the conservativeness of the GGI
    // weighting factors
    //profileBandValues = profileBandValues/srcScalingValues;

    // profileBandValues are now the circumferentially averaged values

    // The profile to dst transfer is done by simply distributing the
    // average value accordingly to the dst weighting factors
    forAll (dstAddr, faceI)
    {
        const labelList& curAddr = dstAddr[faceI];
        const scalarList& curW = dstWeights[faceI];

        forAll (curAddr, bandI)
        {
            dstResultF[faceI] += profileBandValues[curAddr[bandI]]*curW[bandI];

            if (debug <= -200)
            {
                Info<< "bande: " << dstAddr[faceI][bandI]
                    << "  dst valeur: " << dstResultF[faceI]
                    << "  dst weight: " << dstWeights[faceI][bandI] << endl;
            }
        }
    }

    if (debug <= -500)
    {
        error::printStack(Info);

        Info<< "srcF      :  " << srcF << nl
            << "srcAddr   :  " << srcAddr << nl
            << "srcWeights:  " << srcWeights << nl
            << "profileBandValues:  " << profileBandValues << nl
            << "srcScalingValues:  " << srcScalingValues << nl
            << "dstAddr   :  " << dstAddr << nl
            << "dstWeights:  " << dstWeights << nl
            << "dstResultF:  " << dstResultF << nl
            << "srcScalingValues: " << srcScalingValues << endl;
    }

    if (debug <= -999)
    {
        fileName traceFileDir("./mixingPlaneTraceFiles");

        if (!exists(traceFileDir))
        {
            mkDir(traceFileDir);
        }

        struct timeval tod;
        gettimeofday(&tod, NULL);

        //struct timespec tp;
        //clock_gettime(CLOCK_MONOTONIC,  &tp);
        std::ostringstream osBuffer;

        osBuffer
            << Foam::name(tod.tv_sec)
            << "." << std::setfill('0') << std::setw(6)
            <<  Foam::name(tod.tv_usec);

        fileName traceFileName(traceFileDir/"profileValues_" + osBuffer.str());
        OFstream dumpFileSrc(traceFileName + "_orig");
        OFstream dumpFileDst(traceFileName + "_interpolated");
        OFstream dumpFileProfile(traceFileName + "_profile");

        //Foam::error::printStack(Info);

        InfoIn
        (
            "MixingPlaneInterpolation::interpolate"
        )   << "Dumping src profiles to: " << traceFileName + "_orig" << nl
            << "Dumping interpolated profiles to : "
            << traceFileName + "_interpolated" << nl
            << "Dumping profile values to: "
            << traceFileName + "_profile" << endl;

        dumpFileDst << dstResultF << endl;
        dumpFileSrc << srcF << endl;
        dumpFileProfile << profileBandValues << endl;
    }
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
            "const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move master data from 'patch space' to 'profile space'
    // MB: We need this back
    Field<Type> profileFF = transform(masterPatchToProfileT(), patchFF);

    if (debug > 1)
    {
        Info << "MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToSlave: "
            << "patchFF: " << patchFF << endl
            << "profileFF: " << profileFF << endl
            << "masterPatchToProfileT(): " << masterPatchToProfileT() << endl
            << "slaveProfileToPatchT():  " << slaveProfileToPatchT() << endl
            << endl;
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

    interpolate
    (
        profileFF,                     // Master data in 'profile space'
        masterPatchToProfileAddr(),    // From master: compute the average
        masterPatchToProfileWeights(),
        slaveProfileToPatchAddr(),     // To slave we distribute the average from
        slaveProfileToPatchWeights(),  // profile to patch
        result
    );

    // Apply transform to bring the slave field back from 'profile space'
    // to 'patch space'
    transform(result, slaveProfileToPatchT(), result); // MB: We need this back

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
            "const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << patchFF.size()
            << abort(FatalError);
    }

    // Move slave data from 'patch space' to 'profile space'
     Field<Type> profileFF = transform(slavePatchToProfileT(), patchFF);

    if (debug > 1)
    {
        Info << "MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToMaster: "
            << "patchFF: " << patchFF << endl
            << "profileFF: " << profileFF << endl
            << "slavePatchToProfileT(): " << slavePatchToProfileT() << endl
            << "masterProfileToPatchT(): " << masterProfileToPatchT() << endl
            << endl;
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

    if (debug > 1)
    {
        Info << "MixingPlaneInterpolation<MasterPatch, SlavePatch>::masterToMaster: "
            << "patchFF: " << patchFF << endl
            << "profileFF: " << profileFF << endl
            << "masterPatchToProfileT(): " << masterPatchToProfileT() << endl
            << endl;
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
        
    interpolate
    (
        profileFF,                      // Master data in 'profile space'
        masterPatchToProfileAddr(),     // From master: compute the average
        masterPatchToProfileWeights(),
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

    if (debug > 1)
    {
        Info << "MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveToSlave: "
            << "patchFF: " << patchFF << endl
            << "profileFF: " << profileFF << endl
            << "slavePatchToProfileT(): " << slavePatchToProfileT() << endl
            << endl;
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
