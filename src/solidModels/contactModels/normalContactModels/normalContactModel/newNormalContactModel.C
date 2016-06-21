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
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<normalContactModel> normalContactModel::New
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
{
    Info<< "\tNormal contact model: " << name << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "normalContactModel::New(\n"
            "    const word& name,\n"
            "    const fvPatch& patch,\n"
            "    const dictionary& dict,\n"
        "    const label masterPatchID,\n"
        "    const label slavePatchID,\n"
        "    const label masterFaceZoneID,\n"
        "    const label slaveFaceZoneID,\n"
        "    const const PrimitivePatch<face, List, pointField>& "
            "masterFaceZonePatch,\n"
        "    const const PrimitivePatch<face, List, pointField>& "
            "slaveFaceZonePatch\n"
            ")",
            dict
        )   << "Unknown normalContactModel type "
            << name << endl << endl
            << "Valid  normalContactModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<normalContactModel>
      (
       cstrIter()
       (
    name,
    patch,
    dict,
    masterPatchID,
    slavePatchID,
    masterFaceZoneID,
    slaveFaceZoneID,
    masterFaceZonePatch,
    slaveFaceZonePatch
    )
       );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
