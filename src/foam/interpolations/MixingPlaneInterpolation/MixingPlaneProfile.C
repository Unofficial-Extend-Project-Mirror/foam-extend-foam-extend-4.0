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
tmp<pointField>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::computeProfileFromHistograms
(
    const profileHistogram& masterHisto,
    const profileHistogram& slaveHisto,
    const scalar halfSizeBin
) const
{
    // Find min, max bounds
    scalar histoMinValue =
        Foam::min
        (
            masterHisto.begin()->first,
            slaveHisto.begin()->first
        );

    scalar histoMaxValue =
        Foam::max
        (
            (--masterHisto.end())->first,
            (--slaveHisto.end())->first
        );

    // Next, we compare both histograms, leap-frogging from one histogram
    // to the other, everytime jumping to the next largest value from
    // a given position
    std::list<point> leapFrogProfile(0);

    scalar curRvalue = histoMinValue;

    if (masterHisto.find(curRvalue) != masterHisto.end())
    {
        leapFrogProfile.push_back
        (
            *(masterHisto.find(curRvalue)->second.begin())
        );
    }
    else
    {
        leapFrogProfile.push_back
         (
             *(slaveHisto.find(curRvalue)->second.begin())
         );
    }

    do
    {
        profileHistogram::const_iterator nextMaster =
            masterHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        profileHistogram::const_iterator nextSlave  =
            slaveHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        if
        (
            nextMaster == masterHisto.end()
         || nextSlave  == slaveHisto.end()
        )
        {
            // We are done
            if (curRvalue != histoMaxValue)
            {
                if (masterHisto.find(histoMaxValue) != masterHisto.end())
                {
                    leapFrogProfile.push_back
                    (
                        *(masterHisto.find(histoMaxValue)->second.begin())
                    );
                }
                else
                {
                    leapFrogProfile.push_back
                    (
                        *(slaveHisto.find(histoMaxValue)->second.begin())
                    );
                }
            }
            break;
        }

        // Leap frog to the next largest delta
        // Please note that we are not taking into account the number of
        // values at the specific map index (the attribute called 'second'
        // of the map container)
        // So by leap-frogging this way, we might ran into the situation
        // where the next jump will bring us to a radius value shared by only
        // one patch face point, and we might oversee a slightly smaller
        // radius shared by 100s of points.
        // There is a definite advantage of sticking to existing radius values
        // for the profile because this will tend to keep the number of
        // GGI face neighbours low.
        // I am not sure if there is a nice solution to this problem.... MB
        //
        curRvalue = max(nextMaster->first, nextSlave->first);

        if (masterHisto.find(curRvalue) != masterHisto.end())
        {
            leapFrogProfile.push_back
            (
                *(masterHisto.find(curRvalue)->second.begin())
            );
        }
        else
        {
            leapFrogProfile.push_back
            (
                *(slaveHisto.find(curRvalue)->second.begin())
            );
        }
    } while (curRvalue != histoMaxValue);

    // Re-package the data into pointField
    tmp<pointField> tprofile(new pointField(leapFrogProfile.size()));

    pointField& profile = tprofile();

    label pI = 0;

    forAllIter (std::list<point>, leapFrogProfile, lI)
    {
        profile[pI++] = *lI;
    }

    return tprofile;
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::updateProfileHistogram
(
    profileHistogram& histo,
    const point& profileCoord,  // 3D point reference
    const direction dir,        // Sorting dimension 0: x, 1: y, 2: z
    const scalar halfSizeBin    // half size of min width for histogram bins
) const
{
    bool foundNewBin = true;

    scalar keyValue = profileCoord.component(dir);

    forAllIter (profileHistogram, histo, histoI)
    {
        if
        (
            keyValue >= histoI->first - halfSizeBin
         && keyValue < histoI->first + halfSizeBin
        )
        {
            foundNewBin = false;
            histoI->second.push_back(profileCoord);
            break;
        }
    }

    if (foundNewBin)
    {
        std::list<point> initValue;

        initValue.push_back(profileCoord);
        histo.insert
        (
            std::pair<scalar, std::list<point> >(keyValue, initValue)
        );
    }
}


template<class MasterPatch, class SlavePatch>
tmp<pointField>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::calcProfile() const
{
    // First, transform patch points over one single global profile
    // into local coordinates

    pointField masterGlobalProfile =
        cs_.localPosition(masterPatch_.localPoints());

    pointField slaveGlobalProfile =
        cs_.localPosition(slavePatch_.localPoints());

    // Collapse all points in the sweep direction

    const direction sweepDir = sweepAxisSwitch();

    masterGlobalProfile.replace(sweepDir, 0);
    slaveGlobalProfile.replace(sweepDir, 0);

    if (debug)
    {
        InfoIn
        (
            "tmp<pointField> MixingPlaneInterpolation::calcProfile()"
        )   << "masterGlobalProfile: " << masterGlobalProfile << nl
            << "slaveGlobalProfile: " << slaveGlobalProfile << endl;
    }

    // Find the smallest edge length from both patches.
    // This length will control the size of the bins for the histograms
    // we are about to build.
    scalar masterMinEdgeLength = GREAT;
    const edgeList& masterEdgeList = masterPatch_.edges();

    forAll (masterEdgeList, mEi)
    {
        masterMinEdgeLength =
            Foam::min
            (
                masterMinEdgeLength,
                masterEdgeList[mEi].mag(masterPatch_.localPoints())
            );
    }

    scalar slaveMinEdgeLength = GREAT;
    const edgeList& slaveEdgeList = slavePatch_.edges();

    forAll (slaveEdgeList, sEi)
    {
        slaveMinEdgeLength =
            Foam::min
            (
                slaveMinEdgeLength,
                slaveEdgeList[sEi].mag(slavePatch_.localPoints())
            );
    }

    // There is no point classifying data to a resolution smaller than the
    // largest of the two minimum edges found.
    // This will drive the size of the smallest bin needed to construct the
    // z and r histograms
    scalar halfMinSizeBin =
        Foam::max(masterMinEdgeLength, slaveMinEdgeLength)/2.0;

    if (debug)
    {
        InfoIn
        (
            "tmp<pointField> MixingPlaneInterpolation::calcProfile()"
        )   << "halfMinSizeBin: " << halfMinSizeBin << endl;
    }

    // Build master and slave histogram
    profileHistogram masterHistogram;
    profileHistogram slaveHistogram;

    // stackingDir switch
    const direction stackingDir = stackAxisSwitch();

    // Master side
    forAll (masterGlobalProfile, mI)
    {
        updateProfileHistogram
        (
            masterHistogram,
            masterGlobalProfile[mI],
            stackingDir,
            halfMinSizeBin
        );
    }

    // Shadow side
    forAll (slaveGlobalProfile, sI)
    {
        updateProfileHistogram
        (
            slaveHistogram,
            slaveGlobalProfile[sI],
            stackingDir,
            halfMinSizeBin
        );
    }

    if (debug > 1)
    {
        // Write histograms
        forAllIter (profileHistogram, masterHistogram, zHi)
        {
            Info<< "master histo (z, n): (" << zHi->first << " "
                << static_cast<int>(zHi->second.size()) << ")" << endl;
        }

        forAllIter (profileHistogram, slaveHistogram, zHi)
        {
            Info<< "slave histo (z, n): (" << zHi->first << " "
                << static_cast<int>(zHi->second.size()) << ")" << endl;
        }
    }

    pointField profileBeforeClip;

    // Select which type of profile we need
    switch (discretisationType_)
    {
        // To Do: add uniform mixing plane

        case MASTER_PATCH:
        {
            profileBeforeClip = computeProfileFromHistograms
            (
                masterHistogram,
                masterHistogram,
                halfMinSizeBin
            );
        }
        break;

        case SLAVE_PATCH:
        {
            profileBeforeClip = computeProfileFromHistograms
            (
                slaveHistogram,
                slaveHistogram,
                halfMinSizeBin
            );
        }
        break;

        case BOTH_PATCHES:
        {
            profileBeforeClip = computeProfileFromHistograms
            (
                masterHistogram,
                slaveHistogram,
                halfMinSizeBin
            );
        }
        break;

        case UNIFORM:
        case USER_DEFINED:
        default:
        {
            FatalErrorIn
            (
                "tmp<pointField> MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::calcProfile() const"
            )   << "Bad type of mixing plane discretisation: "
                << discretisationNames_[discretisationType_]
                << nl << "Available types are: " << nl
                << MixingPlaneInterpolationName::discretisationNames_
                << abort(FatalError);
        }
    }

    // Remove interpolationProfile_ points located outside of either
    // master/slave patch boundingBox,
    // with the exception of the first and last profile points

    // Martin to work here. HJ, 27/Jan/2011

    boundBox masterBB
    (
        masterGlobalProfile,
        false
    );

    boundBox slaveBB
    (
        slaveGlobalProfile,
        false
    );

    // Expand the bounding box in the sweepAxis-wise direction
    // Note: All points are collapsed to zero in sweepDir
    // It is sufficient to expand the box by 1 in this direction

    masterBB.min().replace(sweepDir, -1);
    masterBB.max().replace(sweepDir, 1);

    slaveBB.min().replace(sweepDir, -1);
    slaveBB.max().replace(sweepDir, 1);

    boundBox globSpanBB
    (
        point(Foam::min(masterBB.min(), slaveBB.min()))
      - point(SMALL, SMALL, SMALL),
        point(Foam::max(masterBB.max(), slaveBB.max()))
      + point(SMALL, SMALL, SMALL)
    );

    if (debug)
    {
        InfoIn
        (
            "tmp<pointField> MixingPlaneInterpolation<MasterPatch, "
            "SlavePatch>::calcProfile() const"
        )   << setprecision(12) << nl
            << "masterBB: " << masterBB << nl
            << "slaveBB: " <<  slaveBB << nl
            << "globSpanBB: " <<  globSpanBB << nl
            << "initial profile values: " << profileBeforeClip << endl;
    }

    // Iterate through profile, removing points located
    // outside of either master/slave BB

    tmp<pointField> tprofile(new pointField(profileBeforeClip.size()));
    pointField& profile = tprofile();
    label curIndex = 0;

    forAll (profileBeforeClip, pI)
    {
        if (globSpanBB.contains(profileBeforeClip[pI]))
        {
            // We keep that profile point
            profile[curIndex] = profileBeforeClip[pI];
            curIndex++;  // Next slot
        }
        else
        {
            if (debug)
            {
                InfoIn
                (
                    "MixingPlaneInterpolation"
                    "<MasterPatch, SlavePatch>::"
                    "removeNonOverlappedProfilePoints"
                )   << setprecision(12) << "   Removing point: "
                    << profileBeforeClip[pI] << endl;
            }
        }
    }

    profile.setSize(curIndex);

    if (profile.size() < 2)
    {
        FatalErrorIn
        (
            "tmp<pointField> MixingPlaneInterpolation<MasterPatch, "
            "SlavePatch>::calcProfile() const"
        )   << "Lost all points in profile: " << profile
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "removeNonOverlappedProfilePoints"
        )   << "cleaned-up profile values: " << profile << endl;
    }

    return tprofile;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::pointField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::interpolationProfile() const
{
    if (interpolationProfile_.size() == 0)
    {
        // Not a user-defined profile: calculate as per algorithm
        interpolationProfile_ = calcProfile();
    }

    return interpolationProfile_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
