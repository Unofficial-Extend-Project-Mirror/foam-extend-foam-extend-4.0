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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const polyMesh& moleculeCloud::mesh() const
{
    return mesh_;
}


// MD solution parameters

inline const moleculeCloud::integrationMethods&
    moleculeCloud::integrationMethod() const
{
    return integrationMethod_;
}


inline scalar moleculeCloud::potentialEnergyLimit() const
{
    return potentialEnergyLimit_;
}


inline label moleculeCloud::nPairPotentials() const
{
    return pairPotentials_.size();
}


inline const labelList& moleculeCloud::removalOrder() const
{
    return removalOrder_;
}


inline const labelListList& moleculeCloud::directInteractionList() const
{
    return directInteractionList_;
}


inline const referredCellList& moleculeCloud::referredInteractionList() const
{
    return referredInteractionList_;
}


inline const labelList&
    moleculeCloud::realCellsWithinRCutMaxOfAnyReferringPatch() const
{
    return realCellsWithinRCutMaxOfAnyReferringPatch_;
}


inline const labelList&
    moleculeCloud::realFacesWithinRCutMaxOfAnyReferringPatch() const
{
    return realFacesWithinRCutMaxOfAnyReferringPatch_;
}


inline const labelList&
    moleculeCloud::realEdgesWithinRCutMaxOfAnyReferringPatch() const
{
    return realEdgesWithinRCutMaxOfAnyReferringPatch_;
}


inline const labelList&
    moleculeCloud::realPointsWithinRCutMaxOfAnyReferringPatch() const
{
    return realPointsWithinRCutMaxOfAnyReferringPatch_;
}


inline const List<sendingReferralList>&
    moleculeCloud::cellSendingReferralLists() const
{
    return cellSendingReferralLists_;
}


inline const List<receivingReferralList>&
    moleculeCloud::cellReceivingReferralLists() const
{
    return cellReceivingReferralLists_;
}


inline label moleculeCloud::nInteractingProcs() const
{
    return cellReceivingReferralLists_.size();
}


inline const pairPotentialList& moleculeCloud::pairPotentials() const
{
    return pairPotentials_;
}


inline const tetherPotentialList& moleculeCloud::tetherPotentials() const
{
    return tetherPotentials_;
}


inline const vector& moleculeCloud::gravity() const
{
    return gravity_;
}


inline const List< DynamicList<molecule*> >&
    moleculeCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
