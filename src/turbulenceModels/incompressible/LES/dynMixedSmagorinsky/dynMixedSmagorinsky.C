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

#include "dynMixedSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynMixedSmagorinsky, 0);
addToRunTimeSelectionTable(LESModel, dynMixedSmagorinsky, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynMixedSmagorinsky::dynMixedSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),
    scaleSimilarity(U, phi, transport),
    dynSmagorinsky(U, phi, transport)
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> dynMixedSmagorinsky::k() const
{
    return
    (
        scaleSimilarity::k()
      + dynSmagorinsky::k()
    );
}


tmp<volScalarField> dynMixedSmagorinsky::epsilon() const
{
    return
    (
        scaleSimilarity::epsilon()
      + dynSmagorinsky::epsilon()
    );
}


tmp<volSymmTensorField> dynMixedSmagorinsky::B() const
{
    return
    (
        scaleSimilarity::B()
      + dynSmagorinsky::B()
    );
}


tmp<volSymmTensorField> dynMixedSmagorinsky::devBeff() const
{
    return
    (
        scaleSimilarity::devBeff()
      + dynSmagorinsky::devBeff()
    );
}


tmp<fvVectorMatrix> dynMixedSmagorinsky::divDevBeff(volVectorField& U) const
{
    return
    (
        scaleSimilarity::divDevBeff(U)
      + dynSmagorinsky::divDevBeff(U)
    );
}


void dynMixedSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    scaleSimilarity::correct(gradU);
    dynSmagorinsky::correct(gradU());
}


bool dynMixedSmagorinsky::read()
{
    if (LESModel::read())
    {
        scaleSimilarity::read();
        dynSmagorinsky::read();

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // namespace incompressible
} // End namespace Foam

// ************************************************************************* //
