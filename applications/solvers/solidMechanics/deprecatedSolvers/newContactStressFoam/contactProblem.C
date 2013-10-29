/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Class describes a multiple body contact problem.  Each individual contact
    is described by a contactPatchPair.  contactProblem handles
    multiple contact updates and sets the boundary conditions on the
    displacement field.

\*---------------------------------------------------------------------------*/

#include "contactProblem.H"
#include "fvMesh.H"
#include "FieldFields.H"
#include "directionMixedFvPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactProblem, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject
contactProblem::contactProblem
(
    volVectorField& U,
    const volTensorField& gradU
)
:
    IOdictionary
    (
        IOobject
        (
            "contactProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    contactPatchPairList(),
    U_(U),
    gradU_(gradU),
    urfValue_(readScalar(lookup("urfValue"))),
    urfTraction_(readScalar(lookup("urfTraction"))),
    urfFraction_(readScalar(lookup("urfFraction")))
{
    // Read contactPatchPairList
    Istream& is = lookup("contacts");

    PtrList<entry> contactEntries(is);

    contactPatchPairList& contacts = *this;

    contacts.setSize(contactEntries.size());

    forAll(contacts, contactI)
    {
        contacts.set
        (
            contactI,
            new contactPatchPair
            (
                contactEntries[contactI].keyword(),
                *this,
                contactEntries[contactI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void contactProblem::correct()
{
    contactPatchPairList& contacts = *this;

    // Create fields for accumulation
    volVectorField::GeometricBoundaryField& Upatches = U().boundaryField();

    FieldField<Field, vector> curTraction(Upatches.size());
    FieldField<Field, vector> newTraction(Upatches.size());
    FieldField<Field, vector> refValue(Upatches.size());
    FieldField<Field, scalar> valueFraction(Upatches.size());

    forAll (Upatches, patchI)
    {
        curTraction.set
        (
            patchI,
            new vectorField(Upatches[patchI].size(), vector::zero)
        );

        newTraction.set
        (
            patchI,
            new vectorField(Upatches[patchI].size(), vector::zero)
        );

        refValue.set
        (
            patchI,
            new vectorField(Upatches[patchI].size(), vector::zero)
        );

        valueFraction.set
        (
            patchI,
            new scalarField(Upatches[patchI].size(), 0)
        );
    }

    // Collect patches involved in contact
    boolList contactPatches(Upatches.size(), false);

    forAll (contacts, contactI)
    {
        contactPatches[contacts[contactI].masterPatch().index()] = true;
        contactPatches[contacts[contactI].slavePatch().index()] = true;
    }

    // Calculate the traction for all involved patches

    // Collect fields
    const volTensorField::GeometricBoundaryField& gradUpatches =
        gradU().boundaryField();

    const surfaceVectorField::GeometricBoundaryField& Apatches =
        mesh().Sf().boundaryField();
    const surfaceScalarField::GeometricBoundaryField& magApatches =
        mesh().magSf().boundaryField();

    // Lookup mu and lambda form object registry
    const volScalarField& mu =
        mesh().objectRegistry::lookupObject<volScalarField>("mu");

    const volScalarField::GeometricBoundaryField& muPatches =
        mu.boundaryField();

    const volScalarField& lambda =
        mesh().objectRegistry::lookupObject<volScalarField>("lambda");

    const volScalarField::GeometricBoundaryField& lambdaPatches =
        lambda.boundaryField();

    forAll (Upatches, patchI)
    {
        if (contactPatches[patchI])
        {
            vectorField nPatch = Apatches[patchI]/magApatches[patchI];

            curTraction[patchI] =
                nPatch &
                (
                    muPatches[patchI]*
                    (
                        gradUpatches[patchI]
                      + gradUpatches[patchI].T()
                    )
                  + I*(lambdaPatches[patchI]*tr(gradUpatches[patchI]))
                );
        }
    }

    // Accumulate contact data and active patches
    forAll (contacts, contactI)
    {
        contacts[contactI].correct
        (
            curTraction,
            newTraction,
            refValue,
            valueFraction
        );
    }

    // Enforce accumulated contact onto the patches
    forAll (Upatches, patchI)
    {
        if (contactPatches[patchI])
        {
            // Cast the patch into direction mixed type
            directionMixedFvPatchVectorField& curUPatch =
                refCast<directionMixedFvPatchVectorField>(Upatches[patchI]);

            // Set the values using under-relaxation
            curUPatch.refValue() =
                (1.0 - urfValue_)*curUPatch.refValue()
              + urfValue_*refValue[patchI];

            // Calculate the gradient from under-relaxad accumulated traction
            vectorField nPatch = Apatches[patchI]/magApatches[patchI];

            curUPatch.refGrad() =
            (
                (1.0 - urfTraction_)*curTraction[patchI]
              + urfTraction_*newTraction[patchI]

              - (nPatch &
                (
                    muPatches[patchI]*gradUpatches[patchI].T()
                  - (
                        muPatches[patchI]
                      + lambdaPatches[patchI]
                    )*gradUpatches[patchI]
                )
                )

              - nPatch*
                (
                    lambdaPatches[patchI]*tr(gradUpatches[patchI])
                )
        
            )/(2.0*muPatches[patchI] + lambdaPatches[patchI]);

            // Set the value fractions
            curUPatch.valueFraction() =
                (1.0 - urfFraction_)*curUPatch.valueFraction()
              + I*urfFraction_*valueFraction[patchI];
        }
    }
}


tmp<volScalarField> contactProblem::contactArea() const
{
    tmp<volScalarField> tca
    (
        new volScalarField
        (
            IOobject
            (
                "contactArea",
                U().time().timeName(),
                U().db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar(0)
        )
    );

    volScalarField& ca = tca();

    // Set contact area boundary
    const contactPatchPairList& contacts = *this;

    forAll (contacts, contactI)
    {
        // Get master contact
        ca.boundaryField()[contacts[contactI].masterPatch().index()] +=
            contacts[contactI].masterTouchFraction();

        // Get slave contact
        ca.boundaryField()[contacts[contactI].slavePatch().index()] +=
            contacts[contactI].slaveTouchFraction();
    }

    return tca;
}


// Return a list of contactPatchPair names
wordList contactProblem::names() const
{
    const contactPatchPairList& contacts = *this;

    wordList t(contacts.size());

    forAll (contacts, contactI)
    {
        t[contactI] = contacts[contactI].name();
    }

    return t;
}


bool contactProblem::read()
{
    if (regIOobject::read())
    {
        urfValue_ = readScalar(lookup("urfValue"));
        urfTraction_ = readScalar(lookup("urfTraction"));
        urfFraction_ = readScalar(lookup("urfFraction"));

        // Decided not to re-read contactPatchPairList.  HJ, 10/Jul/2004
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
