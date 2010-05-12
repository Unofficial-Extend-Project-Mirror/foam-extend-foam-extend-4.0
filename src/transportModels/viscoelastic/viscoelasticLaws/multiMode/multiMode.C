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

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(multiMode, 0);
addToRunTimeSelectionTable(viscoelasticLaw, multiMode, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
multiMode::multiMode
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    viscoelasticLaw(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    nulo_
    (
        IOobject
        (
            "nulo",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    models_()
{
    PtrList<entry> modelEntries(dict.lookup("models"));
    models_.setSize(modelEntries.size());

    forAll (models_, modelI)
    {
        models_.set
        (
            modelI,
            viscoelasticLaw::New
            (
                modelEntries[modelI].keyword(),
                U,
                phi,
                modelEntries[modelI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> multiMode::divTau(volVectorField& U) const
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTau(U);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix() += models_[i].divTau(U);
    }

    return divMatrix;
}


tmp<volSymmTensorField> multiMode::tau() const
{

    tau_ = nulo_; 

    for (label i = 0; i < models_.size(); i++)
    {
        tau_ += models_[i].tau();
    }

    return tau_;

}


void multiMode::correct()
{
    forAll (models_, i)
    {
        Info<< "Model mode "  << i+1 << endl;
        models_[i].correct();
    }

    tau();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
