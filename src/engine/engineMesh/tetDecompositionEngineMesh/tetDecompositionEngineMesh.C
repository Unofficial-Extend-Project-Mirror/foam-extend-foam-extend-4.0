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

#include "tetDecompositionEngineMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"

#include "tetFem.H"

#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(tetDecompositionEngineMesh, 0);

addToRunTimeSelectionTable(engineMesh, tetDecompositionEngineMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from objectRegistry, and read/write options
tetDecompositionEngineMesh::tetDecompositionEngineMesh(const IOobject& io)
:
    engineMesh(io),
    pistonLayers_("pistonLayers", dimLength, 0.0)
{
    if (engineDB_.engineDict().found("pistonLayers"))
    {
        engineDB_.engineDict().lookup("pistonLayers") >> pistonLayers_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetDecompositionEngineMesh::~tetDecompositionEngineMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tetDecompositionEngineMesh::move()
{
    scalar deltaZ = engineDB_.pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;

    // Position of the top of the static mesh layers above the piston
    scalar pistonPlusLayers = pistonPosition_.value() + pistonLayers_.value();

    pointField newPoints = points();

    tetPolyMesh tetMesh(*this);

    // Select the set of boundary condition types.  For symmetry planes
    // and wedge boundaries, the slip condition should be used;
    // otherwise, use the fixedValue
    wordList boundaryTypes
    (
        boundary().size(),
        fixedValueTetPolyPatchScalarField::typeName
    );

    forAll (boundary(), patchI)
    {
        if
        (
            isType<symmetryFvPatch>(boundary()[patchI])
         || isType<wedgeFvPatch>(boundary()[patchI])
         || isType<emptyFvPatch>(boundary()[patchI])
        )
        {
            boundaryTypes[patchI] = slipTetPolyPatchScalarField::typeName;
        }
    }

    tetPointScalarField motionUz
    (
        IOobject
        (
            "motionUz",
            engineDB_.timeName(),
            engineDB_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh,
        dimensionedScalar("0", dimLength, 0),
        boundaryTypes
    );

    motionUz.boundaryField()[pistonIndex_] == deltaZ;

    {
        scalarField linerPoints =
            motionUz.boundaryField()[linerIndex_].patch()
            .localPoints().component(vector::Z);

        motionUz.boundaryField()[linerIndex_] ==
            deltaZ*pos(deckHeight_.value() - linerPoints)
            *(deckHeight_.value() - linerPoints)
            /(deckHeight_.value() - pistonPlusLayers);
    }

    elementScalarField diffusion
    (
        IOobject
        (
            "motionDiffusion",
            engineDB_.timeName(),
            engineDB_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh,
        dimensionedScalar("d", dimless, 1.0)
    );

    const fvPatchList& patches = boundary();

    forAll (patches, patchI)
    {
        const unallocLabelList& fc = patches[patchI].faceCells();

        forAll (fc, fcI)
        {
            diffusion[fc[fcI]] = 2;
        }
    }

    solve(tetFem::laplacian(diffusion, motionUz));

    newPoints.replace
    (
        vector::Z,
        newPoints.component(vector::Z)
      + scalarField::subField
        (
            motionUz.internalField(),
            newPoints.size()
        )
    );

    if (engineDB_.foundObject<surfaceScalarField>("phi"))
    {
        surfaceScalarField& phi = 
            const_cast<surfaceScalarField&>
            (engineDB_.lookupObject<surfaceScalarField>("phi"));

        const volScalarField& rho =
            engineDB_.lookupObject<volScalarField>("rho");

        const volVectorField& U = 
            engineDB_.lookupObject<volVectorField>("U");

        bool absolutePhi = false;
        if (moving())
        {
            phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
            absolutePhi = true;
        }

        movePoints(newPoints);

        if (absolutePhi)
        {
            phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        }
    }

    pistonPosition_.value() += deltaZ;
    scalar pistonSpeed = deltaZ/engineDB_.deltaT().value();

    Info<< "clearance: " << deckHeight_.value() - pistonPosition_.value() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
