/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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
    const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
    const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
)
:
    name_(name),
    patch_(patch),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID),
    masterFaceZoneID_(masterFaceZoneID),
    slaveFaceZoneID_(slaveFaceZoneID),
    slaveContactPointGap_(patch.boundaryMesh().mesh().boundaryMesh()[slavePatchID].nPoints(), 0.0),
    masterToSlaveInterpolatorPtr_
    (
     // new zoneToZoneInterpolation
     new PatchToPatchInterpolation<PrimitivePatch<face, List, pointField>, PrimitivePatch<face, List, pointField> >
     (
      // patch().boundaryMesh().mesh().faceZones()[masterFaceZoneID](), // from zone
      // patch().boundaryMesh().mesh().faceZones()[slaveFaceZoneID](), // to zone
      masterFaceZonePatch,
      slaveFaceZonePatch,
      intersection::algorithmNames_.read(dict.lookup("projectionAlgo")),
      intersection::directionNames_.read(dict.lookup("projectionDir"))
      )
     )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
