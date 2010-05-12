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

Application
    yPlusRAS

Description
    Calculates and reports yPlus for all wall patches, for the specified times.
    Extended version for being able to handle two phase flows using the
    -twoPhase option.

    Frank Albina, 16/Nov/2009

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseMixture.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RASModel/RASModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate single phase Y+
void calcSinglePhaseYPlus
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    volScalarField& yPlus
)
{
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> RASModel
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const fvPatchList& patches = U.mesh().boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            yPlus.boundaryField()[patchi] = RASModel->yPlus(patchi);
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << currPatch.name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }
}

// Calculate two phase Y+
void calcTwoPhaseYPlus
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    volScalarField& yPlus
)
{
    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "gamma");

    autoPtr<incompressible::RASModel> RASModel
    (
        incompressible::RASModel::New(U, phi, twoPhaseProperties)
    );

    const fvPatchList& patches = U.mesh().boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            yPlus.boundaryField()[patchi] = RASModel->yPlus(patchi);
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << currPatch.name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    //Additional option for twoPhase transport model
    argList::validOptions.insert("twoPhase","");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    // Check if two phase model was selected
    bool twoPhase = args.options().found("twoPhase");

    forAll (timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        Info << "Reading field U\n" << endl;
        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info << "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

#           include "createPhi.H"

            if (twoPhase)
            {
                Info<< "Reading field gamma\n" << endl;
                IOobject gammaHeader
                (
                    "gamma",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (gammaHeader.headerOk())
                {
                    volScalarField gamma(gammaHeader, mesh);
                    calcTwoPhaseYPlus(U, phi, yPlus);
                }
                else
                {
                    Info << "    no gamma field!" << endl;
                }
            }
            else
            {
                calcSinglePhaseYPlus(U, phi, yPlus);
            }
        }
        else
        {
            Info  << "    no U field" << endl;
        }

        Info<< "Writing yPlus to field " << yPlus.name() << nl << endl;
        yPlus.write();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
