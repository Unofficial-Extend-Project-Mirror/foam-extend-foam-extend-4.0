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

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistrySolver<CompType, ThermoType> >
Foam::chemistrySolver<CompType, ThermoType>::New
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& compTypeName,
    const word& thermoTypeName
)
{
    word modelName(model.lookup("chemistrySolver"));

    word chemistrySolverType =
        modelName + '<' + compTypeName + ',' + thermoTypeName + '>';

    Info<< "Selecting chemistrySolver " << modelName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistrySolverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        wordList models = dictionaryConstructorTablePtr_->sortedToc();
        forAll(models, i)
        {
            models[i] = models[i].replace
            (
                '<' + compTypeName + ',' + thermoTypeName + '>',
                ""
            );
        }

        FatalErrorIn
        (
            "chemistrySolver::New"
            "("
                "const ODEChemistryModel&, "
                "const word&, "
                "const word&"
            ")"
        )   << "Unknown chemistrySolver type " << modelName
            << nl << nl << "Valid chemistrySolver types are:" << nl
            << models << nl << exit(FatalError);
    }

    return autoPtr<chemistrySolver<CompType, ThermoType> >
        (cstrIter()(model, modelName));
}


// ************************************************************************* //
