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

#include "writeFsiData.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeFsiData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeFsiData,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::writeFsiData::writeFsiData
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
    totalSubIterations_(0)
//     aMeshPtr_(NULL)
{
    if (Pstream::parRun())
    {
        FatalErrorIn("writeFsiData::writeFsiData(...)")
            << "writeFsiData objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

//     const fvMesh& mesh =
//         time_.lookupObject<fvMesh>(regionName_);

//     if (mesh.foundObject<areaScalarField>("h"))
//     {
//         const areaScalarField& h =
//             mesh.lookupObject<areaScalarField>("h");

//         aMeshPtr_ = &h.mesh();
//     }
//     else
//     {
//         FatalErrorIn("writeFsiData::writeFsiData(...)")
//             << " Can not find film thickness filed h"
//                 << abort(FatalError);
//     }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::writeFsiData::start()
{
    return false;
}


bool Foam::writeFsiData::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const fluidSolidInterface& fsi =
        mesh.lookupObject<fluidSolidInterface>("fsiProperties");

    const volScalarField& p =
        mesh.lookupObject<volScalarField>("p");

    if (time_.outputTime())
    {
        OFstream file
        (
            time_.timePath()/"interfacePressure.dat"
        );

        file.precision(12);

        vectorField Cf =
            mesh.boundary()[fsi.fluidPatchIndex()].Cf();

        const scalarField pf = p.boundaryField()[fsi.fluidPatchIndex()];

        file << "x" << tab << "y" << tab << "p" << endl;

        forAll(Cf, faceI)
        {
            file << Cf[faceI].x() << tab << Cf[faceI].y() << tab
                << pf[faceI] << endl;
        }
    }

    if (mesh.foundObject<fluidSolidInterface>("fsiProperties"))
    {
        const fluidSolidInterface& fsi =
            mesh.lookupObject<fluidSolidInterface>("fsiProperties");

        totalSubIterations_ += fsi.outerCorr();
    }

    Info << "Avg num of sub-iterations: "
        << totalSubIterations_/mesh.time().timeIndex() << endl;;

    return false;
}


bool Foam::writeFsiData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
