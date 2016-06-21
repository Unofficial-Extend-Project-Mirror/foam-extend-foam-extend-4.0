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

\*---------------------------------------------------------------------------*/

#include "ePatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<ePatch> ePatch::New
(
    const word& patchType,
    const word& name,
    const label size,
    const label start,
    const label index,
    const eBoundaryMesh& bm
)
{
    if (debug)
    {
        Info<< "ePatch::New(const word&, const word&, const label, "
               "const label, const label, const eBoundaryMesh&) : "
               "constructing ePatch"
            << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(patchType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ePatch::New(const word&, const word&, const label, "
            "const label, const label, const eBoundaryMesh&) "
        )   << "Unknown ePatch type " << patchType << " for patch " << name
            << endl << endl
            << "Valid ePatch types are :" << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<ePatch>(cstrIter()(name, size, start, index, bm));
}


autoPtr<ePatch> ePatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const eBoundaryMesh& bm
)
{
    if (debug)
    {
        Info<< "ePatch::New(const word&, const dictionary&, const label, "
               "const eBoundaryMesh&) : constructing ePatch"
            << endl;
    }

    word patchType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(patchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "ePatch::New(const word&, const dictionary&, "
            "const label, const eBoundaryMesh&)",
            dict
        )   << "Unknown ePatch type " << patchType << endl << endl
            << "Valid ePatch types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<ePatch>(cstrIter()(name, dict, index, bm));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
