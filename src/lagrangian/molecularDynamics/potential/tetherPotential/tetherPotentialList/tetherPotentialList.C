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

\*----------------------------------------------------------------------------*/

#include "tetherPotentialList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetherPotentialList::readTetherPotentialDict
(
    const dictionary& tetherPotentialDict,
    const List<word>& idList,
    const List<label>& tetherIds
)
{
    if (!tetherIds.size())
    {
        Info<< nl << "No tethered molecules found." << endl;

        idMap_.setSize(0);
    }
    else
    {
        Info<< nl << "Building tether potentials." << endl;

        idMap_ = List<label>(idList.size(), -1);

        label tetherMapIndex = 0;

        forAll(tetherIds, t)
        {
            const label tetherId = tetherIds[t];

            word tetherPotentialName = idList[tetherId];

            if (!tetherPotentialDict.found(tetherPotentialName))
            {
                FatalErrorIn("tetherPotentialList::readTetherPotentialDict")
                    << "tether potential specification subDict "
                    << tetherPotentialName << " not found" << nl
                    << abort(FatalError);
            }

            this->set
            (
                tetherMapIndex,
                tetherPotential::New
                (
                    tetherPotentialName,
                    tetherPotentialDict.subDict(tetherPotentialName)
                )
            );

            idMap_[tetherId] = tetherMapIndex;

            tetherMapIndex++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetherPotentialList::tetherPotentialList()
:
    PtrList<tetherPotential>(),
    idMap_()
{}


Foam::tetherPotentialList::tetherPotentialList
(
    const dictionary& idListDict,
    const dictionary& tetherPotentialDict,
    const List<label>& tetherIds
)
:
    PtrList<tetherPotential>(),
    idMap_()
{
    buildPotentials(idListDict, tetherPotentialDict, tetherIds);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetherPotentialList::~tetherPotentialList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetherPotentialList::buildPotentials
(
    const dictionary& idListDict,
    const dictionary& tetherPotentialDict,
    const List<label>& tetherIds
)
{
    setSize(tetherIds.size());

    List<word> idList(idListDict.lookup("idList"));

    readTetherPotentialDict(tetherPotentialDict, idList, tetherIds);
}


const Foam::tetherPotential& Foam::tetherPotentialList::tetherPotentialFunction
(
    const label a
) const
{
    return (*this)[tetherPotentialIndex(a)];
}


Foam::scalar Foam::tetherPotentialList::force
(
    const label a,
    const scalar rITMag
) const
{
    scalar f = (*this)[tetherPotentialIndex(a)].force(rITMag);

    return f;
}


Foam::scalar Foam::tetherPotentialList::energy
(
    const label a,
    const scalar rITMag
) const
{
    scalar e = (*this)[tetherPotentialIndex(a)].energy(rITMag);

    return e;
}


// ************************************************************************* //
