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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "patchAvgTractionHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchAvgTractionHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchAvgTractionHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::patchAvgTractionHistory::writeData()
{
    Info << "Writing average traction history" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volSymmTensorField& sigma =
        mesh.lookupObject<volSymmTensorField>("sigma");

    const symmTensorField& patchSigma =
        sigma.boundaryField()[patchIndex_];

//     const vectorField& patchS = mesh.Sf().boundaryField()[patchIndex_];

    const scalarField& patchMagS = mesh.magSf().boundaryField()[patchIndex_];

//     vector totForce =
//         gSum(patchS & patchSigma);

//     const volTensorField& gradD =
//         mesh.lookupObject<volTensorField>("grad(D)");

    vector totForce =
        gSum
        (
            mesh.Sf().boundaryField()[patchIndex_] & patchSigma
        );

//     vector totForce =
//         gSum
//         (
//             mesh.Sf().boundaryField()[patchIndex_]
//           & (patchSigma & (I+gradD.boundaryField()[patchIndex_]))
//         );

    vector avgTraction = totForce/(gSum(patchMagS) - SMALL);

    Info << "Average traction at patch " << patchName_
        << ": " << avgTraction << "(" << totForce
        << ")" << "\n" << endl;

    if (Pstream::master())
    {
        historyFilePtr_()
            << mesh.time().value()
                << tab << avgTraction.x()
                << tab << avgTraction.y()
                << tab << avgTraction.z()
                << tab << totForce.x()
                << tab << totForce.y()
                << tab << totForce.z() << endl;
    }

    if (time_.outputTime())
    {
        OFstream file
        (
            time_.timePath()/"tensile-sigmayy-sigmaxy-radius.dat"
        );

        file.precision(6);

        const vectorField Cf = mesh.boundary()[patchIndex_].Cf();

        file << "x" << tab << "sigmayy" << tab << "sigmaxy" << endl;
        forAll(Cf, faceI)
        {
            file << Cf[faceI].x() << tab
//                 << Cf[faceI].y() << tab
//                 << Cf[faceI].z() << tab
                << patchSigma[faceI].yy() << tab
                << patchSigma[faceI].xy() << endl;
        }
    }


    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchAvgTractionHistory::patchAvgTractionHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(NULL),
    patchName_(dict.lookup("patchName")),
    patchIndex_(-1)
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    polyPatchID patch(patchName_, mesh.boundaryMesh());

    if (!patch.active())
    {
        FatalErrorIn("patchAvgTractionHistory::patchAvgTractionHistory()")
            << "Patch name " << patchName_ << " not found."
                << abort(FatalError);
    }

    patchIndex_ = patch.index();

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            fileName file("avgTraction_" + patchName_ + ".dat");

            // Open new file at start up
            historyFilePtr_.reset(new OFstream(historyDir/file));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                        << tab << "avgTractionX"
                        << tab << "avgTractionY"
                        << tab << "avgTractionZ"
                        << tab << "totForceX"
                        << tab << "totForceY"
                        << tab << "totForceZ" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchAvgTractionHistory::start()
{
    return writeData();
}


bool Foam::patchAvgTractionHistory::execute()
{
    return writeData();
}


bool Foam::patchAvgTractionHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
