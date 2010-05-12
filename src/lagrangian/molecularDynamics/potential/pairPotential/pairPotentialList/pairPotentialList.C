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

#include "pairPotentialList.H"
#include "OFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairPotentialList::readPairPotentialDict
(
    const dictionary& pairPotentialDict,
    const polyMesh& mesh
)
{
    Info<< nl << "Building pair potentials." << endl;

    rCutMax_ = 0.0;

    for (label a = 0; a < nIds(); ++a)
    {
        word idA = idList_[a];

        for (label b = a; b < nIds(); ++b)
        {
            word idB = idList_[b];

            word pairPotentialName;

            if (a == b)
            {
                if (pairPotentialDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }
                else
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                            << "Pair pairPotential specification subDict "
                            << idA << "-" << idB << " not found"
                            << nl << abort(FatalError);
                }
            }
            else
            {
                if (pairPotentialDict.found(idA + "-" + idB))
                {
                    pairPotentialName = idA + "-" + idB;
                }

                else if (pairPotentialDict.found(idB + "-" + idA))
                {
                    pairPotentialName = idB + "-" + idA;
                }

                else
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " or "
                        << idB << "-" << idA << " not found"
                        << nl << abort(FatalError);
                }

                if
                (
                    pairPotentialDict.found(idA + "-" + idB)
                 && pairPotentialDict.found(idB + "-" + idA)
                )
                {
                    FatalErrorIn("pairPotentialList::buildPotentials") << nl
                        << "Pair pairPotential specification subDict "
                        << idA << "-" << idB << " and "
                        << idB << "-" << idA << " found multiple definition"
                        << nl << abort(FatalError);
                }
            }

            (*this).set
            (
                pairPotentialIndex(a, b),
                pairPotential::New
                (
                    pairPotentialName,
                    pairPotentialDict.subDict(pairPotentialName)
                )
            );

            if ((*this)[pairPotentialIndex(a, b)].rCut() > rCutMax_)
            {
                rCutMax_ = (*this)[pairPotentialIndex(a, b)].rCut();
            }

            if ((*this)[pairPotentialIndex(a, b)].writeTables())
            {
                OFstream ppTabFile(mesh.time().path()/pairPotentialName);

                if
                (
                    !(*this)[pairPotentialIndex(a, b)].writeEnergyAndForceTables
                    (
                        ppTabFile
                    )
                )
                {
                    FatalErrorIn("pairPotentialList::readPairPotentialDict")
                        << "Failed writing to "
                        << ppTabFile.name() << nl
                        << abort(FatalError);
                }
            }
        }
    }

    rCutMaxSqr_ = rCutMax_*rCutMax_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPotentialList::pairPotentialList()
:
    PtrList<pairPotential>(),
    idList_()
{}


Foam::pairPotentialList::pairPotentialList
(
    const dictionary& idListDict,
    const dictionary& pairPotentialDict,
    const polyMesh& mesh
)
:
    PtrList<pairPotential>(),
    idList_()
{
    buildPotentials(idListDict, pairPotentialDict, mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairPotentialList::~pairPotentialList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pairPotentialList::buildPotentials
(
    const dictionary& idListDict,
    const dictionary& pairPotentialDict,
    const polyMesh& mesh
)
{
    idList_ = List<word>(idListDict.lookup("idList"));

    setSize(((idList_.size()*(idList_.size() + 1))/2));

    readPairPotentialDict(pairPotentialDict, mesh);
}


const Foam::pairPotential& Foam::pairPotentialList::pairPotentialFunction
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex(a, b)];
}


bool Foam::pairPotentialList::rCutSqr
(
    const label a,
    const label b,
    const scalar rIJMagSqr
) const
{
    if (rIJMagSqr < rCutSqr(a, b))
    {
        return true;
    }
    else
    {
        return false;
    }
}


Foam::scalar Foam::pairPotentialList::rMin
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex(a, b)].rMin();
}


Foam::scalar Foam::pairPotentialList::dr
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex(a, b)].dr();
}


Foam::scalar Foam::pairPotentialList::rCutSqr
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex(a, b)].rCutSqr();
}


Foam::scalar Foam::pairPotentialList::rCut
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex(a, b)].rCut();
}


Foam::scalar Foam::pairPotentialList::force
(
    const label a,
    const label b,
    const scalar rIJMag
) const
{
    scalar f = (*this)[pairPotentialIndex(a, b)].forceLookup(rIJMag);

    return f;
}


Foam::scalar Foam::pairPotentialList::energy
(
    const label a,
    const label b,
    const scalar rIJMag
) const
{
    scalar e = (*this)[pairPotentialIndex(a, b)].energyLookup(rIJMag);

    return e;
}


// ************************************************************************* //
