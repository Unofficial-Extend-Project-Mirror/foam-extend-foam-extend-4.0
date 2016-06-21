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

#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SRF
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<SRFModel> SRFModel::New
(
    const volVectorField& Urel
)
{
    word SRFModelTypeName;

    // Enclose the creation of the SRFPropertiesDict to ensure it is
    // deleted before the SRFModel is created - otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary SRFPropertiesDict
        (
            IOobject
            (
                "SRFProperties",
                Urel.time().constant(),
                Urel.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        SRFPropertiesDict.lookup("SRFModel") >> SRFModelTypeName;
    }

    Info<< "Selecting SRFModel " << SRFModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(SRFModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "SRFModel::New(const fvMesh&)"
        )   << "Unknown SRFModel type " << SRFModelTypeName
            << nl << nl
            << "Valid SRFModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<SRFModel>(cstrIter()(Urel));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SRF
} // End namespace Foam

// ************************************************************************* //
