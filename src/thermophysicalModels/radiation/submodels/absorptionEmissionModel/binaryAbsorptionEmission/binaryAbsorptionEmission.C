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

#include "binaryAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(binaryAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            binaryAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::binaryAbsorptionEmission::binaryAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    model1_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model1"), mesh)
    ),
    model2_
    (
        absorptionEmissionModel::New(coeffsDict_.subDict("model2"), mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::binaryAbsorptionEmission::~binaryAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::aCont() const
{
    return model1_->aCont() + model2_->aCont();
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::aDisp() const
{
    return model1_->aDisp() + model2_->aDisp();
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::eCont() const
{
    return model1_->eCont() + model2_->eCont();
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::eDisp() const
{
    return model1_->eDisp() + model2_->eDisp();
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::ECont() const
{
    return model1_->ECont() + model2_->ECont();
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::binaryAbsorptionEmission::EDisp() const
{
    return model1_->EDisp() + model2_->EDisp();
}


// ************************************************************************* //
