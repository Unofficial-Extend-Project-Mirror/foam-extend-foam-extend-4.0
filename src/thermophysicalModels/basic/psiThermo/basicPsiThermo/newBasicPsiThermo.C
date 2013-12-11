/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "basicPsiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicPsiThermo> Foam::basicPsiThermo::New
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    word thermoTypeName;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                obj,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("thermoType") >> thermoTypeName;
    }

    Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(thermoTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("basicPsiThermo::New(const fvMesh&, const objectRegistry&)")
            << "Unknown basicPsiThermo type " << thermoTypeName << nl << nl
            << "Valid basicPsiThermo types are:" << nl
            << fvMeshConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicPsiThermo>(cstrIter()(mesh, obj));
}


// ************************************************************************* //
