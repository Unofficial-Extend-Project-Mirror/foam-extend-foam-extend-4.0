/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "forcesFiltered.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forcesFiltered, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forcesFiltered::forcesFiltered
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    forces(name, obr, dict, loadFromFiles),
    filterFieldName_("")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forcesFiltered::~forcesFiltered()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forcesFiltered::read(const dictionary& dict)
{
    if (active_)
    {
        forces::read(dict);

        dict.lookup("filterFieldName") >> filterFieldName_;

        // Check whether filterFieldName_ exists, if not deactivate forces
        if (!obr_.foundObject<volScalarField>(filterFieldName_))
        {
            active_ = false;
            WarningIn("void forces::read(const dictionary& dict)")
                << "Could not find " << filterFieldName_
                << " in database." << nl
                << "    De-activating forces."
               << endl;
        }
    }
}


Foam::forces::forcesMoments Foam::forcesFiltered::calcForces() const
{
    const volScalarField& filterField =
    obr_.lookupObject<volScalarField>(filterFieldName_);

    return forces::calcForces(filterField*devRhoReff());
}


// ************************************************************************* //
