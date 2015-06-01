/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoChemistryModel> Foam::rhoChemistryModel::New
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
{
    word rhoChemistryModelType;
    word thermoTypeName;
    word userModel;

    // Enclose the creation of the chemistrtyProperties to ensure it is
    // deleted before the chemistrtyProperties is created otherwise the
    // dictionary is entered in the database twice
    {
        IOdictionary chemistryPropertiesDict
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                obj,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        chemistryPropertiesDict.lookup("rhoChemistryModel") >> userModel;

        // construct chemistry model type name by inserting first template
        // argument
        label tempOpen = userModel.find('<');
        label tempClose = userModel.find('>');

        word className = userModel(0, tempOpen);
        thermoTypeName = userModel(tempOpen + 1, tempClose - tempOpen - 1);

        rhoChemistryModelType =
            className + '<' + typeName + ',' + thermoTypeName + '>';
    }

    if (debug)
    {
        Info<< "Selecting rhoChemistryModel " << rhoChemistryModelType << endl;
    }
    else
    {
        Info<< "Selecting rhoChemistryModel " << userModel << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(rhoChemistryModelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        if (debug)
        {
            FatalErrorIn("rhoChemistryModelBase::New(const mesh&, const objectRegistry&)")
                << "Unknown rhoChemistryModel type " << rhoChemistryModelType
                << nl << nl << "Valid rhoChemistryModel types are:" << nl
                << fvMeshConstructorTablePtr_->sortedToc() << nl << exit(FatalError);
        }
        else
        {
            wordList models = fvMeshConstructorTablePtr_->sortedToc();
            forAll(models, i)
            {
                models[i] = models[i].replace(typeName + ',', "");
            }

            FatalErrorIn("rhoChemistryModelBase::New(const mesh&, const objectRegistry&)")
                << "Unknown rhoChemistryModel type " << userModel
                << nl << nl << "Valid rhoChemistryModel types are:" << nl
                << models << nl << exit(FatalError);
        }
    }

    return autoPtr<rhoChemistryModel>
        (cstrIter()(mesh, obj, typeName, thermoTypeName));
}


// ************************************************************************* //
