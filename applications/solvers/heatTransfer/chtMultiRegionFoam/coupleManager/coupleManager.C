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

#include "coupleManager.H"
#include "OFstream.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupleManager::coupleManager(const fvPatch& patch)
:
    patch_(patch),
    neighbourRegionName_("undefined-neighbourRegionName"),
    neighbourPatchName_("undefined-neighbourPatchName"),
    neighbourFieldName_("undefined-neighbourFieldName"),
    localRegion_(patch_.boundaryMesh().mesh())
{}


Foam::coupleManager::coupleManager
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    neighbourRegionName_(dict.lookup("neighbourRegionName")),
    neighbourPatchName_(dict.lookup("neighbourPatchName")),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    localRegion_(patch_.boundaryMesh().mesh())
{}


Foam::coupleManager::coupleManager
(
    const coupleManager& cm
)
:
    patch_(cm.patch()),
    neighbourRegionName_(cm.neighbourRegionName()),
    neighbourPatchName_(cm.neighbourPatchName()),
    neighbourFieldName_(cm.neighbourFieldName()),
    localRegion_(patch_.boundaryMesh().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupleManager::~coupleManager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coupleManager::regionOwner() const
{
    const fvMesh& nbrRegion = neighbourRegion();

    const regionProperties& props =
        localRegion_.objectRegistry::parent().lookupObject<regionProperties>
        (
            "regionProperties"
        );

    label myIndex = findIndex(props.fluidRegionNames(), localRegion_.name());
    if (myIndex == -1)
    {
        label i = findIndex(props.solidRegionNames(), localRegion_.name());

        if (i == -1)
        {
            FatalErrorIn("coupleManager::regionOwner() const")
                << "Cannot find region " << localRegion_.name()
                << " neither in fluids " << props.fluidRegionNames()
                << " nor in solids " << props.solidRegionNames()
                << exit(FatalError);
        }
        myIndex = props.fluidRegionNames().size() + i;
    }
    label nbrIndex = findIndex(props.fluidRegionNames(), nbrRegion.name());
    if (nbrIndex == -1)
    {
        label i = findIndex(props.solidRegionNames(), nbrRegion.name());

        if (i == -1)
        {
            FatalErrorIn("coupleManager::regionOwner() const")
                << "Cannot find region " << nbrRegion.name()
                << " neither in fluids " << props.fluidRegionNames()
                << " nor in solids " << props.solidRegionNames()
                << exit(FatalError);
        }
        nbrIndex = props.fluidRegionNames().size() + i;
    }

    return myIndex < nbrIndex;
}


void Foam::coupleManager::checkCouple() const
{
    Info<< "neighbourRegionName_ = " << neighbourRegionName_ << endl;
    Info<< "neighbourPatchName_ = " << neighbourPatchName_ << endl;
    Info<< "neighbourFieldName_ = " << neighbourFieldName_ << endl;

    const fvPatch& nPatch = neighbourPatch();

    if (patch_.size() != nPatch.size())
    {
        FatalErrorIn("Foam::coupleManager::checkCouple()")
            << "Unequal patch sizes:" << nl
            << "    patch name (size) = " << patch_.name()
            << "(" << patch_.size() << ")" << nl
            << "    neighbour patch name (size) = "
            << nPatch.name() << "(" << patch_.size() << ")" << nl
            << abort(FatalError);
    }
}


void Foam::coupleManager::coupleToObj() const
{
    const fvPatch& nPatch = neighbourPatch();

    OFstream obj
    (
         patch_.name() + "_to_" + nPatch.name() + "_couple.obj"
    );
    const vectorField& c1 = patch_.Cf();
    const vectorField& c2 = neighbourPatch().Cf();

    if (c1.size() != c2.size())
    {
        FatalErrorIn("coupleManager::coupleToObj() const")
            << "Coupled patches are of unequal size:" << nl
            << "    patch0 = " << patch_.name()
            << "(" << patch_.size() <<  ")" << nl
            << "    patch1 = " << nPatch.name()
            << "(" << nPatch.size() <<  ")" << nl
            << abort(FatalError);
    }

    forAll(c1, i)
    {
        obj << "v " << c1[i].x() << " " << c1[i].y() << " " << c1[i].z() << nl
            << "v " << c2[i].x() << " " << c2[i].y() << " " << c2[i].z() << nl
            << "l " << (2*i + 1) << " " << (2*i + 2) << endl;
    }
}


void Foam::coupleManager::writeEntries(Ostream& os) const
{
    os.writeKeyword("neighbourRegionName");
    os << neighbourRegionName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatchName");
    os << neighbourPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldName");
    os << neighbourFieldName_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
