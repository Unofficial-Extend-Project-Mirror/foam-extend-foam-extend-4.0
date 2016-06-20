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
    frictionContactModel

\*---------------------------------------------------------------------------*/

#include "frictionContactModel.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(frictionContactModel, 0);
defineRunTimeSelectionTable(frictionContactModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frictionContactModel::frictionContactModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID
)
:
    name_(name),
    patch_(patch),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID),
    masterFaceZoneID_(masterFaceZoneID),
    slaveFaceZoneID_(slaveFaceZoneID),
    stickSlipFaces_(patch.boundaryMesh()[slavePatchID].size(), 0.0)
{}


frictionContactModel::frictionContactModel(const frictionContactModel& fm)
:
    name_(fm.name_),
    patch_(fm.patch_),
    masterPatchID_(fm.masterPatchID_),
    slavePatchID_(fm.slavePatchID_),
    masterFaceZoneID_(fm.masterFaceZoneID_),
    slaveFaceZoneID_(fm.slaveFaceZoneID_),
    stickSlipFaces_(fm.stickSlipFaces_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
